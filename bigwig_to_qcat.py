#!/usr/bin/env python3
"""
bigwig_to_qcat.py
=================
Full pipeline: BigWig files (one per category) -> 200 bp binned KL-divergence
scores -> pyGenomeTracks-compatible qcat.bgz + .tbi index.

USAGE
-----
    python3 bigwig_to_qcat.py \
        --bw  atac.bw rnaseq_plus.bw rnaseq_minus.bw ctcf.bw cohesin.bw \
               h3k27ac.bw \
        --out  my_experiment.qcat.bgz \
        --genome mm10

    # or point at a directory:
    python3 bigwig_to_qcat.py --bw-dir ./bigwigs/ --out my_experiment.qcat.bgz

    # testing mode -- only process intervals in a BED file:
    python3 bigwig_to_qcat.py --bw *.bw --out test.qcat.bgz --regions my_regions.bed

    # sheet mode -- run one or more BigWig->qcat jobs from a TSV:
    python3 bigwig_to_qcat.py --sheet samples.tsv --out placeholder.qcat.bgz

Pass between 5 and 15 BigWig files. The first N categories from the hardcoded
list will be used automatically -- no configuration needed.
In --sheet mode, each row can define its own BigWig list and output name.
The CLI --out value is required by argparse but ignored when --sheet is used.

DEPENDENCIES
------------
    Python 3.8+
    pip3 install pyBigWig pysam numpy scipy

    Optional (for --stats-plots):
    pip3 install matplotlib

WHAT THIS SCRIPT DOES
---------------------
1. Reads the mm10 chromosome sizes (built-in, or from --chrom-sizes).
2. Tiles every chromosome into 200 bp bins (or only bins overlapping --regions).
3. For each bin, reads the mean signal from each BigWig -> raw signal vector.
4. Normalises the signal vector to a probability distribution P (clip negatives,
   add pseudocount, divide by row sum).
5. Computes the background distribution Q as the mean of P across all bins.
6. Computes per-state S1 epilogos KL scores for each bin:
       score_i = P_i * log2( P_i / Q_i )
   Negative values are clamped to 0.
7. Writes a tab-separated qcat file:
       chr  start  end  id:N,qcat:[[score1,state1],...],raw:[s1,s2,...,sK]
   where raw values are pre-normalization signal (abs, clipped, before
   pseudocount), one value per track in input order.
8. Bgzip-compresses and tabix-indexes the output.
9. Writes a companion categories JSON and example .ini for pyGenomeTracks.

NOTES ON KL SCORING
-------------------
The epilogos S1 metric is designed for discrete chromatin state calls. Here it
is adapted for continuous BigWig signal by treating the normalised signal vector
as an observed probability distribution P and computing per-state KL
contributions against the background Q. Bins where a state's signal is
unusually high relative to its background average receive large positive scores.
"""

import argparse
import bisect
import json
import math
import time
import multiprocessing
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

_json = json

# ---------------------------------------------------------------------------
# Full 15-state category list.
# The script uses the first N entries, where N = number of BigWig files
# supplied (must be between MIN_STATES and MAX_STATES).
# Edit names or colors here as needed.
# ---------------------------------------------------------------------------
ALL_CATEGORIES = [
    ("ATAC",                       "#be92e0"),
    ("RNAseq +",                   "#6495ed"),
    ("RNAseq -",                   "#1a3a8f"),
    ("CTCF",                       "#ff2200"),
    ("Cohesin",                    "#8b0000"),
    ("H3K27ac",                    "#00e676"),
    ("Bivalent/Poised TSS",        "#cd5c5c"),
    ("Flanking Bivalent TSS/Enh",  "#e9967a"),
    ("Bivalent Enhancer",          "#bdb76b"),
    ("Repressed PolyComb",         "#808080"),
    ("Weak Repressed PolyComb",    "#c0c0c0"),
    ("Quiescent/Low",              "#ffffff"),
]
MIN_STATES = 5
MAX_STATES = len(ALL_CATEGORIES)

# Placeholder built-in feature sets for mm10 around the TCRb locus.
# TODO: replace with genome-wide coordinates from ENCODE/Ensembl BED files.
BUILTIN_FEATURE_SETS = {
    "promoters": [
        ("chr6", 41500000, 41502000),
        ("chr6", 41800000, 41802000),
    ],
    "gene_bodies": [
        ("chr6", 41470000, 41560000),
        ("chr6", 41760000, 41870000),
    ],
    "all_genes": [
        ("chr6", 41470000, 41560000),
        ("chr6", 41760000, 41870000),
        ("chr6", 41500000, 41502000),
        ("chr6", 41800000, 41802000),
    ],
    "tss_2kb": [
        ("chr6", 41498000, 41502000),
        ("chr6", 41798000, 41802000),
    ],
    "ctcf_sites": [
        ("chr6", 41254000, 41254200),
        ("chr6", 41618000, 41618200),
        ("chr6", 42192000, 42192400),
    ],
}

# 1-based indices of tracks whose BigWig signal is negative (negative strand).
# For these tracks, abs() is applied before normalization so the signal
# magnitude is preserved. The INI writer also flips their display axes.
NEGATIVE_STRAND_STATES = {3}  # RNAseq -


def _detect_unmappable_for_chrom(task):
    """Worker for per-chromosome unmappable detection."""
    (
        chrom,
        chrom_len,
        intervals,
        interval_bin_counts,
        n_bins,
        all_bw_paths,
        n_samples,
        n_tracks_per_sample,
        zero_frac,
        min_samples,
        bin_size,
        signal_threshold,
    ) = task

    try:
        import pyBigWig
    except ImportError:
        raise RuntimeError("pyBigWig is required for --detect-unmappable")

    if n_bins <= 0:
        return chrom, np.zeros(0, dtype=bool), 0

    combo_zero_counts = np.zeros(n_bins, dtype=np.uint16)
    sample_track_zero_counts = None
    if min_samples is not None:
        sample_track_zero_counts = np.zeros((n_samples, n_bins), dtype=np.uint16)

    for si, bw_path in all_bw_paths:
        try:
            bw = pyBigWig.open(bw_path)
        except Exception as exc:
            raise RuntimeError(f"Could not open BigWig: {bw_path}") from exc
        if bw is None:
            raise RuntimeError(f"Could not open BigWig (open returned None): {bw_path}")

        try:
            if intervals is not None:
                raw = []
                for (s, e), n_subbins in zip(intervals, interval_bin_counts):
                    if n_subbins <= 0:
                        continue
                    vals = bw.stats(chrom, s, e, type="mean", nBins=n_subbins)
                    if vals is None:
                        vals = [None] * n_subbins
                    raw.extend(vals)
            else:
                # Fast path: one bulk values() read + numpy binning is ~10x
                # faster than stats(nBins=...), which issues a summary query per
                # bin. Only zero/nonzero detection matters here, so binned means
                # are sufficient. Falls back to stats() on any error.
                try:
                    _v = bw.values(chrom, 0, chrom_len, numpy=True)
                    _v = np.nan_to_num(_v, nan=0.0).astype(np.float32)
                    _pad = (-len(_v)) % bin_size
                    if _pad:
                        _v = np.concatenate([_v, np.zeros(_pad, dtype=np.float32)])
                    _binned = _v.reshape(-1, bin_size).mean(axis=1)
                    raw = _binned[:n_bins].tolist()
                    if len(raw) < n_bins:
                        raw = raw + [None] * (n_bins - len(raw))
                except Exception:
                    raw = bw.stats(chrom, 0, chrom_len, type="mean", nBins=n_bins)
                    if raw is None:
                        raw = [None] * n_bins
        except Exception:
            # Missing chrom or inaccessible values -> treat as all-zero.
            zero_flags = np.ones(n_bins, dtype=bool)
        else:
            zero_flags = np.zeros(n_bins, dtype=bool)
            for bi, v in enumerate(raw):
                if v is None:
                    zero_flags[bi] = True
                    continue
                fv = float(v)
                if math.isnan(fv) or abs(fv) < signal_threshold:
                    zero_flags[bi] = True
        finally:
            bw.close()

        combo_zero_counts += zero_flags.astype(np.uint16)
        if sample_track_zero_counts is not None:
            sample_track_zero_counts[si] += zero_flags.astype(np.uint16)

    frac_zero = combo_zero_counts.astype(np.float32) / float(len(all_bw_paths))
    flagged = frac_zero >= float(zero_frac)

    if sample_track_zero_counts is not None:
        n_samples_zero = np.sum(
            sample_track_zero_counts == n_tracks_per_sample,
            axis=0,
        )
        flagged = flagged & (n_samples_zero >= int(min_samples))

    n_flagged = int(flagged.sum())
    return chrom, flagged, n_flagged


def load_categories_yaml(path):
    """
    Load category definitions from a YAML or JSON file.

    Returns (categories_list, negative_strand_states) where:
      categories_list        : list of (name, color) tuples
      negative_strand_states : set of 1-based indices where negative_strand=true
    
    Supports two formats:
    
    1. Standard YAML/JSON (list of dicts):
       categories:
         - name: "ATAC"
           color: "#be92e0"
         - name: "CTCF"
           color: "#ff2200"
    
    2. Numeric-keyed JSON (from bearmon):
       {
         "categories": {
           "1": ["ATAC", "#be92e0"],
           "2": ["CTCF", "#ff2200"]
         }
       }
    """
    try:
        import yaml
    except ImportError:
        sys.exit(
            "ERROR: PyYAML is required for --categories.  "
            "Install with:  pip install pyyaml"
        )
    with open(path) as fh:
        doc = yaml.safe_load(fh)
    if not isinstance(doc, dict) or "categories" not in doc:
        sys.exit(f"ERROR: {path} must contain a top-level 'categories' key.")
    
    entries = doc["categories"]
    categories_list = []
    negative_strand_states = set()
    
    # Handle two different structures for the categories list
    
    if isinstance(entries, list):
        # Format 1: list of dicts
        for i, entry in enumerate(entries, start=1):
            if not isinstance(entry, dict):
                sys.exit(
                    f"ERROR: category {i} in {path} is not a dict (got {type(entry).__name__}: {repr(entry)})"
                )
            name = entry.get("name")
            color = entry.get("color")
            if not name or not color:
                sys.exit(
                    f"ERROR: category {i} in {path} is missing 'name' or 'color'."
                )
            categories_list.append((name, color))
            if entry.get("negative_strand", False):
                negative_strand_states.add(i)
    
    elif isinstance(entries, dict):
        # Format 2: numeric-key dict with [name, color] values (from bearmon JSON)
        # Sort by numeric key to ensure consistent ordering
        try:
            sorted_keys = sorted(entries.keys(), key=lambda x: int(x))
        except ValueError:
            # Keys are not numeric; just sort alphabetically
            sorted_keys = sorted(entries.keys())
        
        for idx, key in enumerate(sorted_keys, start=1):
            entry = entries[key]
            if isinstance(entry, (list, tuple)) and len(entry) >= 2:
                name = entry[0]
                color = entry[1]
                if not isinstance(name, str) or not isinstance(color, str):
                    sys.exit(
                        f"ERROR: category '{key}' in {path} has invalid types."
                    )
                categories_list.append((name, color))
            elif isinstance(entry, dict):
                name = entry.get("name")
                color = entry.get("color")
                if not name or not color:
                    sys.exit(
                        f"ERROR: category '{key}' in {path} is missing 'name' or 'color'."
                    )
                categories_list.append((name, color))
                if entry.get("negative_strand", False):
                    negative_strand_states.add(idx)
            else:
                sys.exit(
                    f"ERROR: category '{key}' in {path} has unexpected format."
                )
    else:
        sys.exit(
            f"ERROR: 'categories' in {path} must be a list or dict, "
            f"got {type(entries).__name__}."
        )
    
    if not categories_list:
        sys.exit(f"ERROR: no categories found in {path}")
    
    return categories_list, negative_strand_states


def build_categories(n, categories_list=None):
    """Return a categories dict for the first n states."""
    src = categories_list if categories_list is not None else ALL_CATEGORIES
    return {str(i + 1): list(src[i]) for i in range(n)}


# ---------------------------------------------------------------------------
# Built-in mm10 chromosome sizes (GRCm38)
# ---------------------------------------------------------------------------
MM10_CHROM_SIZES = {
    "chr1":  195471971,
    "chr2":  182113224,
    "chr3":  160039680,
    "chr4":  156508116,
    "chr5":  151834684,
    "chr6":  149736546,
    "chr7":  145441459,
    "chr8":  129401213,
    "chr9":  124595110,
    "chr10": 130694993,
    "chr11": 122082543,
    "chr12": 120129022,
    "chr13": 120421639,
    "chr14": 124902244,
    "chr15": 104043685,
    "chr16":  98207768,
    "chr17":  94987271,
    "chr18":  90702639,
    "chr19":  61431566,
    "chrX":  171031299,
    "chrY":   91744698,
}

PSEUDOCOUNT = 1e-6
BIN_SIZE    = 200    # bp
MIN_SIGNAL  = 0.01   # bins with total raw signal below this are zeroed out

# Per-track noise-floor defaults (in CPM). Tracks whose per-bin signal falls
# below their per-track floor are zeroed out BEFORE the bin's signal vector
# is converted to a probability distribution. This prevents sparse/noisy
# tracks (e.g. RNA-seq) from contributing noise-floor signal to the
# compositional Q vector in bins where other denser tracks (ATAC, ChIP-seq)
# carry the bin past --min-signal. Use --min-signal-per-track to override.
MIN_SIGNAL_PER_TRACK_DEFAULTS = {
    "RNAseq +": 0.1,
    "RNAseq -": 0.1,
    "ATAC":     0.05,
    "CTCF":     0.05,
    "Cohesin":  0.05,
    "H3K27ac":  0.05,
}


def parse_min_signal_per_track(spec, all_categories):
    """Parse a --min-signal-per-track spec string.

    Accepts:
      * the literal "default" (case-insensitive) -> returns the built-in
        MIN_SIGNAL_PER_TRACK_DEFAULTS dict, filtered to track names present
        in all_categories.
      * a comma-separated list of "TrackName=value" pairs. Track names with
        spaces are supported when quoted at the shell level. Track names
        not present in all_categories raise an error to catch typos.

    Returns a dict {track_name: float_floor}.
    """
    if spec is None:
        return {}
    if isinstance(spec, dict):
        return dict(spec)
    spec = str(spec).strip()
    if spec == "" or spec.lower() == "none":
        return {}
    valid_names = {name for name, _ in all_categories}
    if spec.lower() == "default":
        return {n: v for n, v in MIN_SIGNAL_PER_TRACK_DEFAULTS.items()
                if n in valid_names}
    out = {}
    for entry in spec.split(","):
        entry = entry.strip()
        if not entry:
            continue
        if "=" not in entry:
            sys.exit(
                f"ERROR: --min-signal-per-track entry '{entry}' is not in "
                f"NAME=VALUE form. Example: 'RNAseq +=0.1,ATAC=0.05'."
            )
        name, val = entry.split("=", 1)
        name = name.strip()
        try:
            fval = float(val.strip())
        except ValueError:
            sys.exit(
                f"ERROR: --min-signal-per-track value for '{name}' is not a "
                f"number: '{val}'."
            )
        if fval < 0:
            sys.exit(
                f"ERROR: --min-signal-per-track value for '{name}' is "
                f"negative: {fval}."
            )
        if name not in valid_names:
            sys.exit(
                f"ERROR: --min-signal-per-track track name '{name}' is not "
                f"in the active category list. Valid names: "
                f"{sorted(valid_names)}."
            )
        out[name] = fval
    return out


def load_floors_from_tsv(tsv_path, sample_name, all_categories):
    """Load per-track floors from a TSV produced by an earlier percentile run.

    The TSV format is the same as the one written by --write-floors-tsv:
        sample <TAB> track <TAB> floor [<TAB> percentile <TAB> source]

    The TSV may include rows for many samples. Only rows whose 'sample'
    column matches sample_name (or the literal '*' / 'ALL') are loaded.

    Rationale: percentile-of-nonzero floors are invariant under circular
    shift, so floors computed once on the real BigWigs can be reused for
    all permutation rounds without recomputation, saving the percentile
    sampling pass (~4 min per scoring run).

    Parameters
    ----------
    tsv_path : path-like
        Path to a TSV file with columns: sample, track, floor [, ...].
        Header row is required (any column order, case-insensitive header
        names).
    sample_name : str
        Which sample's floors to load. Rows where sample equals this name
        (or '*' / 'ALL') are included. Sample-specific rows take precedence
        over wildcard rows when both are present for the same track.
    all_categories : list of (name, color) tuples
        Active track names. Floors for tracks not in this list are silently
        ignored (so a multi-sample TSV with extra tracks won't error).

    Returns
    -------
    floors : dict {track_name: float_floor}
        Per-track floors for sample_name. Empty if the TSV has no matching
        rows.
    """
    import csv
    valid_names = {name for name, _ in all_categories}
    sample_floors = {}
    wildcard_floors = {}
    sample_lc = str(sample_name).strip()
    try:
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter="\t", lineterminator="\n")
            if reader.fieldnames is None:
                sys.exit(f"ERROR: --floors-tsv '{tsv_path}' has no header row.")
            # Resolve column names case-insensitively
            colmap = {h.strip().lower(): h for h in reader.fieldnames}
            if "sample" not in colmap or "track" not in colmap \
                    or "floor" not in colmap:
                sys.exit(
                    f"ERROR: --floors-tsv '{tsv_path}' must have columns: "
                    f"sample, track, floor. Found: "
                    f"{list(reader.fieldnames)}"
                )
            s_col = colmap["sample"]
            t_col = colmap["track"]
            f_col = colmap["floor"]
            for row in reader:
                s = (row.get(s_col) or "").strip()
                t = (row.get(t_col) or "").strip()
                v = (row.get(f_col) or "").strip()
                if not s or not t or not v:
                    continue
                if t not in valid_names:
                    continue
                try:
                    fv = float(v)
                except ValueError:
                    print(f"  WARNING: --floors-tsv: non-numeric floor for "
                          f"{s}/{t}: '{v}' (skipped)", file=sys.stderr)
                    continue
                if fv < 0:
                    print(f"  WARNING: --floors-tsv: negative floor for "
                          f"{s}/{t}: {fv} (skipped)", file=sys.stderr)
                    continue
                if s == sample_lc:
                    sample_floors[t] = fv
                elif s in ("*", "ALL", "all"):
                    wildcard_floors[t] = fv
    except FileNotFoundError:
        sys.exit(f"ERROR: --floors-tsv: file not found: {tsv_path}")
    except Exception as ex:
        sys.exit(f"ERROR: --floors-tsv: failed to read '{tsv_path}': {ex}")

    # Merge: sample-specific wins over wildcard
    out = dict(wildcard_floors)
    out.update(sample_floors)
    return out


def write_floors_tsv_file(floors_dict, sample_name, out_path, percentile=None):
    """Write per-track floors to a TSV for reuse in subsequent runs."""
    import csv
    pct_str = f"{percentile:g}" if percentile is not None else ""
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "track", "floor", "percentile", "source"])
        for track, floor in sorted(floors_dict.items()):
            src = f"p{pct_str}_of_nonzero" if pct_str else "manual"
            writer.writerow([sample_name, track, f"{floor:.6g}",
                             pct_str, src])


def apply_per_track_noise_floor(signal_matrix, categories, floors_dict):
    """Zero out per-bin entries below each track's noise floor.

    Modifies signal_matrix in-place and returns the count of (bin, track)
    cells that were zeroed.

    Parameters
    ----------
    signal_matrix : (num_bins, num_states) np.ndarray
        Raw signal matrix. Negative-strand columns should already be
        abs()-converted (this function assumes the matrix represents
        non-negative magnitudes).
    categories : list of (name, color) tuples
        Track names in column order (column i -> categories[i][0]).
    floors_dict : dict {track_name: float_floor}
        Per-track minimum-signal thresholds. Tracks not present in
        floors_dict are not masked.

    Returns
    -------
    n_zeroed : int
        Total number of (bin, track) cells set to zero by this masking pass.
    """
    if not floors_dict:
        return 0
    n_zeroed = 0
    num_cols = signal_matrix.shape[1]
    for col, (name, _color) in enumerate(categories):
        if col >= num_cols:
            break
        floor = floors_dict.get(name)
        if floor is None or floor <= 0:
            continue
        below = signal_matrix[:, col] < floor
        n = int(below.sum())
        if n:
            signal_matrix[below, col] = 0.0
            n_zeroed += n
    return n_zeroed


def compute_percentile_floors_per_track(bw_paths, categories, percentile,
                                        chrom_sizes, chroms,
                                        regions=None,
                                        neg_strand_states=None,
                                        blacklist=None,
                                        max_sample_bins=200_000,
                                        collect_stats=False):
    """Compute per-track noise floors as the Nth percentile of non-zero values.

    Uses random spatial sampling rather than a whole-genome walk. For
    a representative percentile estimate, we sample up to max_sample_bins
    random non-overlapping 200 bp bin positions distributed across the
    processed chromosomes (or across the user-supplied regions), then read
    those bin means from each bigwig. This is ~50-200x faster than walking
    every bin on every chromosome and is statistically indistinguishable
    for percentile estimation (the standard error of the Nth percentile
    from 100k samples is well below 1% of the underlying value).

    For each track:
      * read the sampled bin values from the bigwig
      * apply negative-strand abs() and blacklist masking (matching the
        scoring pass)
      * collect non-zero magnitudes
      * return the Nth percentile of that pooled distribution as the floor

    Parameters
    ----------
    bw_paths : list of paths
        BigWig files, one per track, in column order.
    categories : list of (name, color) tuples
        Track names in column order.
    percentile : float
        Percentile of non-zero values to use as the per-track floor.
        E.g. 5.0 -> 5th percentile of non-zero values. Set to 0.0 to
        disable (returns an empty dict).
    chrom_sizes : dict {chrom: length}
    chroms : list of chrom names to process
    regions : dict {chrom: [(s,e), ...]} or None
        If provided, only bins overlapping these regions are sampled.
    neg_strand_states : list of 1-based state indices
        State indices whose values should be abs()-converted before
        percentile computation (RNA-seq minus strand etc.).
    blacklist : dict {chrom: [(s,e), ...]} or None
        Blacklisted intervals; sampled bins inside are excluded.
    max_sample_bins : int
        Cap on the number of random bin positions used per track for
        percentile computation. Default 200,000. Increase for more
        precise percentile estimates; decrease for faster runs.
    collect_stats : bool
        If True, also compute and return a per-track signal distribution
        summary (n_nonzero, frac_nonzero, p1, p5, p10, p25, p50, p75, p90,
        p95, p99, mean, max).

    Returns
    -------
    floors : dict {track_name: float_floor}
        Per-track Nth-percentile-of-nonzero floor. Tracks with no non-zero
        signal (e.g., a missing bigwig) are absent from the dict.
    signal_stats : dict {track_name: dict} (only when collect_stats=True)
        Per-track raw-signal distribution summary computed on the same
        pooled non-zero values used for the floor (after blacklist+abs
        masking, after sub-sampling). Returned alongside floors as a tuple
        when collect_stats is True.
    """
    if percentile is None or float(percentile) <= 0:
        if collect_stats:
            return {}, {}
        return {}
    import pyBigWig

    pct = float(percentile)
    neg_set = set(neg_strand_states) if neg_strand_states else set()
    rng = np.random.default_rng(seed=0)

    # Step 1: build the master sample of (chrom, start, end) bin positions
    # once, used for all tracks. This guarantees comparable percentile
    # estimates across tracks (same sampling support).
    if regions is not None:
        # Build candidate-bin list from regions (typically small)
        sampleable_bins = []
        for chrom in chroms:
            chrom_regions = regions.get(chrom, [])
            if not chrom_regions:
                continue
            for (rs, re_) in chrom_regions:
                rs = (rs // BIN_SIZE) * BIN_SIZE
                for b_start in range(rs, re_, BIN_SIZE):
                    b_end = b_start + BIN_SIZE
                    if b_end <= re_:
                        sampleable_bins.append((chrom, b_start, b_end))
        total_candidate_bins = len(sampleable_bins)
        if total_candidate_bins == 0:
            print("  WARNING: percentile-floor: no candidate bins in regions",
                  file=sys.stderr)
            if collect_stats:
                return {}, {}
            return {}
        if total_candidate_bins <= max_sample_bins:
            sampled = sampleable_bins
        else:
            idx = rng.choice(total_candidate_bins, size=max_sample_bins,
                             replace=False)
            sampled = [sampleable_bins[i] for i in idx]
    else:
        # Genome-wide: sample by drawing random bin starts proportional to
        # chromosome length (without materializing the full bin list).
        chrom_lens = []
        for chrom in chroms:
            chrom_len = chrom_sizes.get(chrom)
            if chrom_len is None:
                continue
            n_bins = chrom_len // BIN_SIZE
            if n_bins > 0:
                chrom_lens.append((chrom, n_bins))
        total_bins_genome = sum(n for _, n in chrom_lens)
        if total_bins_genome == 0:
            if collect_stats:
                return {}, {}
            return {}
        if total_bins_genome <= max_sample_bins:
            # Small genome / few chroms -- materialize everything
            sampled = []
            for chrom, n_bins in chrom_lens:
                for i in range(n_bins):
                    sampled.append((chrom, i * BIN_SIZE, (i + 1) * BIN_SIZE))
        else:
            # Stratified random sampling per chromosome
            sampled = []
            for chrom, n_bins in chrom_lens:
                n_take = int(round(max_sample_bins * n_bins / total_bins_genome))
                if n_take <= 0:
                    continue
                n_take = min(n_take, n_bins)
                idxs = rng.choice(n_bins, size=n_take, replace=False)
                for i in idxs:
                    sampled.append((chrom, int(i) * BIN_SIZE,
                                    (int(i) + 1) * BIN_SIZE))

    n_sampled_bins = len(sampled)
    print(f"  sampling {n_sampled_bins:,} random bins across "
          f"{len(chroms)} chromosomes for percentile estimation...",
          flush=True)

    # Pre-compute blacklist mask for sampled bins (one pass over all bins).
    if blacklist is not None:
        bl_mask = np.zeros(n_sampled_bins, dtype=bool)
        # Group sampled bins by chrom to batch blacklist queries
        by_chrom = {}
        for i, (chrom, s, e) in enumerate(sampled):
            by_chrom.setdefault(chrom, []).append((i, s, e))
        for chrom, items in by_chrom.items():
            chrom_bl = blacklist.get(chrom, [])
            if not chrom_bl:
                continue
            for (bi, bs, be) in items:
                for (xs, xe) in chrom_bl:
                    if xs < be and xe > bs:
                        bl_mask[bi] = True
                        break
    else:
        bl_mask = None

    # Group sampled bins by chrom once for efficient per-track reads.
    bins_by_chrom = {}
    for i, (chrom, s, e) in enumerate(sampled):
        bins_by_chrom.setdefault(chrom, []).append((i, s, e))

    floors = {}
    stats = {}
    for col, (name, _color) in enumerate(categories):
        if col >= len(bw_paths):
            break
        bw_path = bw_paths[col]
        try:
            bw = pyBigWig.open(str(bw_path))
        except Exception as ex:
            print(f"  WARNING: percentile-floor: could not open "
                  f"{bw_path}: {ex}", file=sys.stderr)
            continue

        vals_arr = np.full(n_sampled_bins, np.nan, dtype=np.float32)
        try:
            for chrom, items in bins_by_chrom.items():
                # Read each sampled bin individually. pyBigWig's stats()
                # with nBins=1 is fast (one IO per call) and handles
                # bigwig zooms internally.
                for (bi, bs, be) in items:
                    try:
                        v = bw.stats(chrom, bs, be, type="mean", nBins=1)[0]
                    except Exception:
                        continue
                    if v is None:
                        continue
                    vals_arr[bi] = float(v)
        finally:
            bw.close()

        # Negative-strand: abs
        if (col + 1) in neg_set:
            vals_arr = np.abs(vals_arr)
        else:
            vals_arr = np.where(np.isfinite(vals_arr),
                                np.clip(vals_arr, 0.0, None), np.nan)

        # Apply blacklist mask
        if bl_mask is not None:
            vals_arr = vals_arr.copy()
            vals_arr[bl_mask] = np.nan

        # Restrict to non-zero, finite values
        finite_mask = np.isfinite(vals_arr)
        nonzero_mask = finite_mask & (vals_arr > 0)
        pooled = vals_arr[nonzero_mask]
        if pooled.size == 0:
            print(f"  WARNING: percentile-floor: no non-zero signal for "
                  f"track '{name}' ({bw_path}); skipping", file=sys.stderr)
            if collect_stats:
                stats[name] = {
                    "n_bins_total": int(n_sampled_bins),
                    "n_nonzero": 0,
                    "n_sampled": int(n_sampled_bins),
                    "frac_nonzero": 0.0,
                }
            continue
        floor = float(np.percentile(pooled, pct))
        floors[name] = floor
        print(f"  percentile floor: {name:<20} "
              f"p{pct:g}-of-nonzero = {floor:.4g}  "
              f"(from {pooled.size:,} non-zero of {int(finite_mask.sum()):,} "
              f"sampled bins; frac_nonzero={pooled.size/max(int(finite_mask.sum()), 1):.3f})",
              flush=True)
        if collect_stats:
            pct_values = [1, 5, 10, 25, 50, 75, 90, 95, 99]
            qs = np.percentile(pooled, pct_values)
            stats[name] = {
                "n_bins_total": int(n_sampled_bins),
                "n_nonzero": int(pooled.size),
                "n_sampled": int(finite_mask.sum()),
                "frac_nonzero": (
                    pooled.size / max(int(finite_mask.sum()), 1)
                ),
                "mean": float(np.mean(pooled)),
                "max": float(np.max(pooled)),
                **{f"p{p}": float(v) for p, v in zip(pct_values, qs)},
            }
    if collect_stats:
        return floors, stats
    return floors


def write_signal_stats_tsv(stats_dict, sample_name, out_path, percentile):
    """Write per-track raw-signal distribution stats to TSV.

    One row per track. Columns include n_bins_total, n_nonzero,
    frac_nonzero, mean, max, and percentiles p1, p5, p10, p25, p50, p75,
    p90, p95, p99. The percentile actually used as the noise-floor (the
    `percentile` argument) is highlighted in the floor_used column.
    """
    import csv
    pct_keys = ["p1", "p5", "p10", "p25", "p50", "p75", "p90", "p95", "p99"]
    floor_col = f"p{int(percentile)}" if (
        percentile is not None and float(percentile).is_integer()
    ) else None
    with open(out_path, "w", newline="") as f:
        fieldnames = (["sample", "track", "n_bins_total", "n_nonzero",
                       "n_sampled", "frac_nonzero", "mean", "max"]
                      + pct_keys + ["floor_percentile", "floor_value"])
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for track in sorted(stats_dict.keys()):
            row = {"sample": sample_name, "track": track,
                   "floor_percentile": (
                       f"p{percentile:g}" if percentile is not None else ""
                   )}
            row.update(stats_dict[track])
            if floor_col is not None and floor_col in stats_dict[track]:
                row["floor_value"] = stats_dict[track][floor_col]
            for k in pct_keys:
                if k not in row:
                    row[k] = ""
            writer.writerow(row)


def plot_signal_stats(stats_dict, out_path, sample_name="", percentile=None):
    """Render a per-track signal-distribution diagnostic figure.

    Two panels:
      * left: percentile bars per track (p5/p25/p50/p75/p95/p99, log scale)
      * right: frac_nonzero per track (bar chart)
    Plus a horizontal marker at the per-track floor (== p[percentile]).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed; skipping signal-stats plot)")
        return
    if not stats_dict:
        return
    tracks = sorted(stats_dict.keys())
    pcts = ["p5", "p25", "p50", "p75", "p95", "p99"]

    fig, axes = plt.subplots(1, 2, figsize=(max(10, len(tracks) * 1.4), 5),
                             gridspec_kw={"width_ratios": [2.2, 1.0]})

    # Panel 1: percentile points per track on log scale
    ax = axes[0]
    cmap = plt.get_cmap("viridis")
    colors_pcts = [cmap(i / (len(pcts) - 1)) for i in range(len(pcts))]
    for i, p in enumerate(pcts):
        vals = [stats_dict[t].get(p, np.nan) for t in tracks]
        ax.plot(range(len(tracks)), vals, marker="o", linestyle="-",
                color=colors_pcts[i], label=p, linewidth=1.2)
    # Highlight the floor percentile (if it lines up with a plotted percentile)
    floor_key = (f"p{int(percentile)}"
                 if (percentile is not None and float(percentile).is_integer())
                 else None)
    if floor_key is not None:
        floor_vals = [stats_dict[t].get(floor_key, np.nan) for t in tracks]
        if any(np.isfinite(v) for v in floor_vals):
            ax.plot(range(len(tracks)), floor_vals, marker="X", linestyle="",
                    markersize=14, markerfacecolor="red",
                    markeredgecolor="black", markeredgewidth=1.0,
                    label=f"floor ({floor_key})", zorder=5)
    ax.set_yscale("log")
    ax.set_xticks(range(len(tracks)))
    ax.set_xticklabels(tracks, rotation=45, ha="right")
    ax.set_ylabel("Signal (track-specific units; log scale)")
    ax.set_title(f"Per-track non-zero signal distribution"
                 + (f"  --  {sample_name}" if sample_name else ""))
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3, which="both", axis="y")

    # Panel 2: frac_nonzero per track
    ax = axes[1]
    fracs = [stats_dict[t].get("frac_nonzero", 0.0) for t in tracks]
    ax.bar(range(len(tracks)), fracs, color="steelblue")
    ax.set_xticks(range(len(tracks)))
    ax.set_xticklabels(tracks, rotation=45, ha="right")
    ax.set_ylabel("Fraction of bins with non-zero signal")
    ax.set_ylim(0, 1.0)
    ax.set_title("Track sparsity")
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            sizes[parts[0]] = int(parts[1])
    return sizes


def bins_for_chrom(chrom_len, bin_size=BIN_SIZE):
    """Return (start, end) pairs for all full+partial bins on a chromosome."""
    starts = range(0, chrom_len, bin_size)
    return [(s, min(s + bin_size, chrom_len)) for s in starts]


def load_bins_bed(path):
    """
    Load a consensus segmentation BED into {chrom: [(start, end), ...]} with each
    chromosome's intervals sorted by start. Used for opt-in adaptive (variable-
    width) binning: when supplied, these intervals replace the fixed grid and the
    background Q is computed width-weighted. Extra columns are ignored.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 3:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if end > start:
                out.setdefault(chrom, []).append((start, end))
    for chrom in out:
        out[chrom].sort(key=lambda se: se[0])
    return out


def load_regions_bed(path):
    """
    Parse regions from either:
      1) standard BED-like lines: chrom  start  end [...]
      2) regions TSV lines:       name   chr:start-end [...]

    Returns a dict: chrom -> sorted list of (start, end).
    Intervals are snapped to 200 bp bin boundaries (start rounded down,
    end rounded up). Overlapping intervals are merged.
    """
    def _parse_region_string(region_s):
        region_s = region_s.strip()
        if ":" not in region_s or "-" not in region_s:
            raise ValueError("not a region string")
        chrom, coords = region_s.split(":", 1)
        start_s, end_s = coords.replace(",", "").split("-", 1)
        return chrom, int(start_s), int(end_s)

    regions = defaultdict(list)
    with open(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith("#") or low.startswith("track") or low.startswith("browser"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 2:
                continue

            # Skip common headers for either BED or regions-template TSV.
            if (parts[0].lower() in {"chrom", "chr", "name"} or
                    (len(parts) > 1 and parts[1].lower() == "region")):
                continue

            chrom = None
            start_i = None
            end_i = None

            # Mode A: regions-template TSV (name, region, ...)
            # e.g. "v1\tchr6:40872100-40908238\t3000\tV1-region"
            if len(parts) >= 2 and ":" in parts[1] and "-" in parts[1]:
                try:
                    chrom, start_i, end_i = _parse_region_string(parts[1])
                except ValueError:
                    pass

            # Mode B: standard BED-like (chrom, start, end, ...)
            if chrom is None and len(parts) >= 3:
                try:
                    chrom = parts[0]
                    start_i = int(parts[1])
                    end_i = int(parts[2])
                except ValueError:
                    chrom = None

            if chrom is None:
                print(
                    f"  WARNING: skipping unrecognized region line {line_no} in {path}: {line}"
                )
                continue

            start = (start_i // BIN_SIZE) * BIN_SIZE
            end   = math.ceil(end_i / BIN_SIZE) * BIN_SIZE
            regions[chrom].append((start, end))

    # Merge overlapping / adjacent intervals per chrom
    merged = {}
    for chrom, ivs in regions.items():
        ivs.sort()
        out = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def bins_for_regions(region_list, bin_size=BIN_SIZE):
    """
    Given a list of (start, end) intervals (already bin-snapped),
    return all 200 bp bins that fall within them.
    """
    bins = []
    for (rs, re) in region_list:
        s = rs
        while s < re:
            bins.append((s, s + bin_size))
            s += bin_size
    return bins


def load_blacklist(bed_path):
    """
    Load a BED blacklist and return {chrom: [(start, end), ...]}.
    Accepts 3+ column BED; only chrom, start, end are used.
    """
    blacklist = defaultdict(list)
    with open(bed_path) as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith("#") or low.startswith("track") or low.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                print(
                    f"  WARNING: skipping malformed blacklist line {line_no} in {bed_path}",
                    file=sys.stderr,
                )
                continue
            if end <= start:
                continue
            blacklist[chrom].append((start, end))

    out = {}
    for chrom, ivs in blacklist.items():
        ivs.sort()
        out[chrom] = ivs
    return out


def bins_overlapping_blacklist(bins, blacklist, chrom):
    """Return boolean mask where True marks bins overlapping blacklist intervals."""
    mask = np.zeros(len(bins), dtype=bool)
    intervals = blacklist.get(chrom, []) if blacklist is not None else []
    if not intervals:
        return mask

    starts = [iv[0] for iv in intervals]
    ends = [iv[1] for iv in intervals]

    for i, (bs, be) in enumerate(bins):
        idx = bisect.bisect_left(starts, int(be))
        for j in range(idx - 1, -1, -1):
            if ends[j] <= int(bs):
                break
            mask[i] = True
            break
    return mask


def detect_unmappable_bins(
    sheet_jobs,
    chrom_sizes,
    regions=None,
    zero_frac=0.90,
    min_samples=None,
    bin_size=BIN_SIZE,
    signal_threshold=0.01,
    workers=1,
):
    """
    Identify bins with near-zero signal across all tracks and samples.

    Returns dict: {chrom: bool ndarray(len=n_bins)}.
    """
    try:
        import pyBigWig
    except ImportError:
        sys.exit(
            "ERROR: pyBigWig is required for --detect-unmappable. "
            "Install with: pip3 install pyBigWig"
        )

    print("Detecting unmappable bins from BigWig data...", file=sys.stderr)

    all_bw_paths = []
    n_samples = len(sheet_jobs)
    for si, job in enumerate(sheet_jobs):
        for bw_path in job["bw_paths"]:
            all_bw_paths.append((si, str(bw_path)))

    # Fail fast before any per-chromosome work so missing/unreadable BigWigs
    # do not get interpreted as all-zero signal.
    failed = check_bigwig_paths([p for _, p in all_bw_paths])
    if failed:
        sys.exit(
            f"\nAborted: fix the {len(failed)} missing/unreadable BigWig file(s) above "
            "before running --detect-unmappable."
        )

    n_combos = len(all_bw_paths)
    if n_combos == 0:
        return {}

    print(
        f"  {n_samples} samples, {n_combos} total (track x sample) combinations",
        file=sys.stderr,
    )

    tracks_per_sample = [len(job["bw_paths"]) for job in sheet_jobs]
    if len(set(tracks_per_sample)) != 1:
        sys.exit(
            "ERROR: --detect-unmappable currently requires each sheet row "
            "to have the same number of BigWig tracks."
        )
    n_tracks_per_sample = tracks_per_sample[0]

    if regions:
        print(
            "  WARNING: --detect-unmappable is running in region-restricted smoke-test mode.",
            file=sys.stderr,
        )
        print(
            "  The resulting BED only covers the requested regions and is not a genome-wide blacklist.",
            file=sys.stderr,
        )

    workers = max(1, int(workers))

    # Precompute per-chromosome bin geometry once to avoid repeated work.
    chrom_bins = {}
    if regions:
        for chrom in regions.keys():
            ivs = regions.get(chrom, [])
            if not ivs:
                continue
            interval_bin_counts = [max(0, int((e - s) // bin_size)) for s, e in ivs]
            n_bins = int(sum(interval_bin_counts))
            chrom_bins[chrom] = {
                "intervals": ivs,
                "interval_bin_counts": interval_bin_counts,
                "n_bins": n_bins,
            }

    chrom_iter = list(regions.keys()) if regions else list(chrom_sizes.keys())
    tasks = []
    for chrom in chrom_iter:
        chrom_len = chrom_sizes.get(chrom)
        if chrom_len is None:
            continue
        if regions and chrom not in regions:
            continue

        if regions:
            meta = chrom_bins.get(chrom)
            if meta is None:
                continue
            n_bins = int(meta["n_bins"])
            intervals = meta["intervals"]
            interval_bin_counts = meta["interval_bin_counts"]
        else:
            n_bins = math.ceil(chrom_len / bin_size)
        if n_bins == 0:
            continue

        tasks.append(
            (
                chrom,
                chrom_len,
                intervals if regions else None,
                interval_bin_counts if regions else None,
                n_bins,
                all_bw_paths,
                n_samples,
                n_tracks_per_sample,
                zero_frac,
                min_samples,
                bin_size,
                signal_threshold,
            )
        )

    result = {}
    if workers > 1 and len(tasks) > 1:
        n_workers = min(workers, len(tasks))
        print(f"  Running detection in parallel with {n_workers} workers", file=sys.stderr)
        with multiprocessing.Pool(processes=n_workers) as pool:
            for chrom, flagged, n_flagged in pool.imap_unordered(
                _detect_unmappable_for_chrom,
                tasks,
            ):
                if n_flagged > 0:
                    print(
                        f"  {chrom}: {n_flagged:,} bins flagged as unmappable",
                        file=sys.stderr,
                    )
                result[chrom] = flagged
    else:
        for task in tasks:
            chrom, flagged, n_flagged = _detect_unmappable_for_chrom(task)
            if n_flagged > 0:
                print(f"  {chrom}: {n_flagged:,} bins flagged as unmappable", file=sys.stderr)
            result[chrom] = flagged

    total_flagged = int(sum(int(v.sum()) for v in result.values()))
    print(f"  Total: {total_flagged:,} bins flagged genome-wide", file=sys.stderr)
    return result


def write_blacklist_bed(flagged_dict, chrom_sizes, out_path, bin_size=BIN_SIZE,
                        regions=None):
    """Write detected unmappable bins as a sorted BED file."""
    with open(out_path, "w") as fh:
        fh.write("# Auto-detected unmappable regions from Bearing\n")
        fh.write("# Generated by bigwig_to_qcat.py --detect-unmappable\n")
        chrom_iter = sorted(regions.keys()) if regions else sorted(chrom_sizes.keys())
        for chrom in chrom_iter:
            mask = flagged_dict.get(chrom)
            if mask is None:
                continue
            if regions and chrom in regions:
                chrom_bins = bins_for_regions(regions[chrom], bin_size=bin_size)
            else:
                chrom_len = chrom_sizes[chrom]
                chrom_bins = bins_for_chrom(chrom_len, bin_size=bin_size)
            for bi, flagged in enumerate(mask[:len(chrom_bins)]):
                if bool(flagged):
                    s, e = chrom_bins[bi]
                    fh.write(f"{chrom}\t{s}\t{e}\n")
    print(f"  Detected blacklist written: {out_path}", file=sys.stderr)


def _compute_chrom_cache(chrom, chrom_len, bw_paths, normalize_tracks,
                         normalize_method, neg_strand_states,
                         regions_for_chrom, temp_dir, blacklist=None,
                         categories=None, min_signal_per_track=None,
                         cohort_ref=None, bins_override=None):
    """Compute per-chromosome P matrix and cache to a temp .npz file.

    Returns:
      (chrom, n_bins, P_sum, q_denom, npz_path, n_pertrack_zeroed)

    q_denom is n_bins for fixed/region bins (so Q stays the unweighted mean) and
    the sum of bin widths when bins_override (adaptive grid) is given (so Q is
    width-weighted). The cached file contains arrays: bins (n_bins x 2 int32), P,
    raw_clipped, raw_abs.
    """
    import pyBigWig

    if bins_override is not None:
        bins = bins_override
    elif regions_for_chrom is not None:
        bins = bins_for_regions(regions_for_chrom)
    else:
        bins = bins_for_chrom(chrom_len)
    n = len(bins)

    num_states = len(bw_paths)
    signal_matrix = np.zeros((n, num_states), dtype=np.float32)

    for si, bw_path in enumerate(bw_paths):
        try:
            with pyBigWig.open(str(bw_path)) as bw:
                signal_matrix[:, si] = mean_signal_in_bins(bw, chrom, bins)
        except Exception as e:
            print(f"  WARNING: could not read {chrom} from BigWig {si+1}: {e}")

    for neg_idx in neg_strand_states:
        col = neg_idx - 1
        if col < num_states:
            signal_matrix[:, col] = np.abs(signal_matrix[:, col])

    if blacklist is not None:
        bl_mask = bins_overlapping_blacklist(bins, blacklist, chrom)
        n_blacklisted = int(bl_mask.sum())
        if n_blacklisted > 0:
            signal_matrix[bl_mask, :] = 0.0
            print(
                f"  {chrom}: zeroed {n_blacklisted:,} blacklisted bins",
                file=sys.stderr,
            )

    raw_abs = np.clip(signal_matrix, 0.0, None).astype(np.float32)

    # Per-track noise-floor masking: zero out cells where a single track's
    # signal is below its per-track floor. Applied BEFORE normalization and
    # BEFORE signals_to_prob so the masked-out track does not contribute to
    # the bin's compositional Q vector.
    n_pertrack_zeroed = 0
    if min_signal_per_track and categories is not None:
        n_pertrack_zeroed = apply_per_track_noise_floor(
            signal_matrix, categories, min_signal_per_track
        )

    if normalize_tracks:
        signal_matrix = normalize_signal_matrix(
            signal_matrix,
            method=normalize_method,
            cohort_ref=cohort_ref,
        )

    P = signals_to_prob(signal_matrix)
    raw_clipped = np.clip(signal_matrix, 0.0, None).astype(np.float32)

    if bins_override is not None:
        widths = np.array([e - s for (s, e) in bins], dtype=np.float64)
        P_sum = (P * widths[:, None]).sum(axis=0)
        q_denom = float(widths.sum())
    else:
        P_sum = P.sum(axis=0)
        q_denom = float(n)

    npz_path = str(Path(temp_dir) / f"{chrom}.npz")
    np.savez_compressed(npz_path,
                        bins=np.array(bins, dtype=np.int32),
                        P=P.astype(np.float32),
                        raw_clipped=raw_clipped,
                        raw_abs=raw_abs)

    return chrom, n, P_sum, q_denom, npz_path, n_pertrack_zeroed


def _score_chrom_from_cache(chrom, q, npz_path, min_signal, normalize_score,
                            start_id, out_tmp_path, score_method="kl"):
    """Score a chromosome (from cached npz) and write qcat rows."""
    import numpy as np

    data = np.load(npz_path, allow_pickle=False)
    bins = data["bins"]
    P = data["P"]
    raw_clipped = data["raw_clipped"]
    raw_abs = data["raw_abs"]

    scores, n_masked = kl_scores_per_bin(
        P,
        q,
        # Floor on the PRE-normalization raw signal (raw_abs), not the
        # post-normalization raw_clipped, so min_signal gates true data
        # presence and the scorable-bin set is identical whether or not
        # normalization is applied. In an un-normalized run raw_abs ==
        # raw_clipped, so this does not change default behavior.
        raw_signal_matrix=raw_abs,
        min_signal=min_signal,
        normalize_score=normalize_score,
        score_method=score_method,
    )

    chrom_sum = float(scores.sum())
    state_totals = scores.sum(axis=0)

    # Write qcat file for this chromosome
    with open(out_tmp_path, "w") as fh:
        bin_id = start_id
        for i in range(len(bins)):
            s, e = int(bins[i][0]), int(bins[i][1])
            row_scores = scores[i]
            pairs = sorted(
                ([float(f"{sc:.6g}"), state_idx]
                 for state_idx, sc in enumerate(row_scores, start=1)),
                key=lambda x: -x[0],
            )
            raw_vec = raw_clipped[i].tolist()
            qcat_col = "id:" + str(bin_id) + ",qcat:" + _json.dumps(pairs, separators=(",", ":"))
            fh.write(f"{chrom}\t{s}\t{e}\t{qcat_col}\n")
            bin_id += 1

    return chrom, int(n_masked), chrom_sum, state_totals, len(bins), out_tmp_path


def mean_signal_in_bins(bw, chrom, bins):
    """
    Return a float array of mean BigWig signal for each bin.
    Missing / NaN values become 0.

    Fast path: read the spanned range once with values() and average per bin
    with numpy. Calling bw.stats() once per bin (the previous approach) issues
    millions of summary queries genome-wide and dominated scoring runtime.
    Falls back to per-bin stats() on any error for robustness.
    """
    n = len(bins)
    vals = np.zeros(n, dtype=np.float32)
    if n == 0:
        return vals
    span_start = int(bins[0][0])
    span_end = int(bins[-1][1])
    try:
        chrom_len = bw.chroms(chrom)
        if chrom_len is not None:
            span_end = min(span_end, int(chrom_len))
    except Exception:
        pass
    # Bulk read with a few retries: under heavy concurrent load on a shared
    # filesystem the pyBigWig reader can fail transiently; a brief retry
    # avoids killing the whole scoring job over a momentary read hiccup.
    raw = None
    for attempt in range(3):
        try:
            raw = bw.values(chrom, span_start, span_end, numpy=True)
            break
        except Exception:
            if attempt < 2:
                time.sleep(0.5 * (attempt + 1))
            else:
                raw = None
    if raw is not None:
        try:
            raw = np.nan_to_num(raw, nan=0.0).astype(np.float32)
            for i, (s, e) in enumerate(bins):
                a = int(s) - span_start
                b = int(e) - span_start
                if a < 0:
                    a = 0
                if b > len(raw):
                    b = len(raw)
                if b > a:
                    vals[i] = raw[a:b].mean()
            return vals
        except Exception:
            pass
    # Fallback: per-bin stats (slower but robust)
    for i, (s, e) in enumerate(bins):
        try:
            v = bw.stats(chrom, s, e, type="mean")[0]
            vals[i] = 0.0 if (v is None or math.isnan(v)) else float(v)
        except Exception:
            vals[i] = 0.0
    return vals


def quantile_normalize_columns(signal_matrix):
    """
    Quantile-normalize each column (track) of a (num_bins x num_states) matrix
    so that all tracks share the same global signal distribution.

    This removes the dynamic range imbalance between broad tracks (e.g.
    H3K27ac, RNAseq) and focal tracks (e.g. CTCF) before the
    joint probability normalization, ensuring no single track dominates Q
    simply because it is broadly active genome-wide.

     Strategy (rank-based quantile normalization):
        1. For each column, replace values with their rank-based quantiles
            derived from the sorted mean distribution across all columns,
            including zeros.
        2. Zero values stay zero -- only non-zero signal is redistributed.

     Because zeros are included in the reference distribution, weak non-zero
     values can be mapped down to zero. This is kept as a legacy option.
    """
    mat = np.clip(signal_matrix, 0.0, None).astype(np.float64)
    n_bins, n_states = mat.shape

    # Compute the reference distribution: mean of sorted values across all columns
    sorted_cols = np.sort(mat, axis=0)           # sort each column independently
    ref_distrib = sorted_cols.mean(axis=1)        # mean at each rank position

    # For each column, map values to the reference distribution by rank
    normalized = np.zeros_like(mat)
    for j in range(n_states):
        col = mat[:, j]
        # Only normalize non-zero entries; leave zeros as zero
        nonzero_mask = col > 0
        if nonzero_mask.sum() == 0:
            continue
        # Rank non-zero values (ties get average rank)
        from scipy.stats import rankdata
        ranks = rankdata(col[nonzero_mask], method="average")
        # Map ranks to reference distribution quantiles
        # ranks are 1-based; scale to [0, n_bins-1] index
        idx = ((ranks - 1) / (nonzero_mask.sum() - 1 + 1e-300) * (n_bins - 1)).astype(int)
        idx = np.clip(idx, 0, n_bins - 1)
        normalized[nonzero_mask, j] = ref_distrib[idx]

    return normalized


def nonzero_quantile_normalize_columns(signal_matrix):
    """
    Quantile-normalize each column using only non-zero values to build the
    reference distribution.

    This preserves low non-zero values better than the legacy method because
    the smallest non-zero values map to the smallest non-zero reference
    quantiles instead of being mixed with the zero mass.
    """
    mat = np.clip(signal_matrix, 0.0, None).astype(np.float64)
    n_bins, n_states = mat.shape
    normalized = np.zeros_like(mat)

    nonzero_cols = [mat[mat[:, j] > 0, j] for j in range(n_states)]
    nonzero_cols = [np.sort(col) for col in nonzero_cols if col.size > 0]
    if not nonzero_cols:
        return normalized

    ref_len = max(len(col) for col in nonzero_cols)
    ref_probs = np.linspace(0.0, 1.0, ref_len)
    ref_distrib = np.mean(
        [np.quantile(col, ref_probs) for col in nonzero_cols],
        axis=0,
    )

    for j in range(n_states):
        col = mat[:, j]
        nonzero_mask = col > 0
        n_nonzero = int(nonzero_mask.sum())
        if n_nonzero == 0:
            continue
        if n_nonzero == 1:
            normalized[nonzero_mask, j] = ref_distrib[0]
            continue

        from scipy.stats import rankdata

        ranks = rankdata(col[nonzero_mask], method="average")
        idx = ((ranks - 1) / (n_nonzero - 1) * (ref_len - 1)).astype(int)
        idx = np.clip(idx, 0, ref_len - 1)
        normalized[nonzero_mask, j] = ref_distrib[idx]

    return normalized


def cohort_quantile_normalize_columns(signal_matrix, cohort_ref):
    """Map each track column onto a COHORT-WIDE per-track reference distribution.

    Unlike the within-sample methods (which build the reference from this
    sample's own tracks), the reference here is supplied per track from a
    pooled-across-samples quantile grid (see build_cohort_quantile_reference.py).
    Mapping a sample's track onto the shared per-track reference equalizes
    between-sample dynamic-range differences for the same assay.

    cohort_ref : (n_states, L) float array -- per-track reference quantile grid,
                 column order matching signal_matrix columns.
    """
    from scipy.stats import rankdata
    mat = np.clip(signal_matrix, 0.0, None).astype(np.float64)
    n_bins, n_states = mat.shape
    normalized = np.zeros_like(mat)
    for j in range(n_states):
        if j >= cohort_ref.shape[0]:
            normalized[:, j] = mat[:, j]
            continue
        ref = np.asarray(cohort_ref[j], dtype=np.float64)
        ref_len = ref.size
        if ref_len == 0 or not ref.any():
            # no cohort reference for this track -> leave values unchanged
            normalized[:, j] = mat[:, j]
            continue
        col = mat[:, j]
        nonzero_mask = col > 0
        n_nonzero = int(nonzero_mask.sum())
        if n_nonzero == 0:
            continue
        if n_nonzero == 1:
            normalized[nonzero_mask, j] = ref[ref_len // 2]
            continue
        ranks = rankdata(col[nonzero_mask], method="average")
        idx = ((ranks - 1) / (n_nonzero - 1) * (ref_len - 1)).astype(int)
        idx = np.clip(idx, 0, ref_len - 1)
        normalized[nonzero_mask, j] = ref[idx]
    return normalized


def normalize_signal_matrix(signal_matrix, method="nonzero-quantile",
                            cohort_ref=None):
    """Apply the requested normalization method.

    method "cohort-quantile" requires cohort_ref (n_states, L); the other
    methods normalize within the sample across its tracks.
    """
    if method == "cohort-quantile":
        if cohort_ref is None:
            raise ValueError("cohort-quantile requires a cohort reference")
        return cohort_quantile_normalize_columns(signal_matrix, cohort_ref)
    if method == "nonzero-quantile":
        return nonzero_quantile_normalize_columns(signal_matrix)
    if method == "quantile":
        return quantile_normalize_columns(signal_matrix)
    raise ValueError(f"Unknown normalization method: {method}")


def signals_to_prob(signal_matrix):
    """
    Convert a (num_bins x num_states) raw-signal matrix to a probability
    distribution per bin (rows sum to 1).
    Clip negatives -> add pseudocount -> divide by row sum.
    """
    mat = np.clip(signal_matrix, 0.0, None) + PSEUDOCOUNT
    row_sums = mat.sum(axis=1, keepdims=True)
    return mat / row_sums


def kl_scores_per_bin(P, Q, raw_signal_matrix=None, min_signal=MIN_SIGNAL,
                      normalize_score=False, score_method="kl"):
    """
    Compute per-state per-bin BEARING scores.

    score_method="kl" (default): clamped per-state KL contribution
        score_{b,i} = max( P_{b,i} * log2(P_{b,i} / Q_i), 0 )
    Unbounded above; negative (depleted) states clamped to 0.

    score_method="jsd": bounded per-state Jensen-Shannon contribution
        m_i      = (P_{b,i} + Q_i) / 2
        jsd_{b,i}= 0.5 * (P_{b,i} log2(P_{b,i}/m_i) + Q_i log2(Q_i/m_i))
    kept only where P_{b,i} > Q_i (same enrichment-only convention as KL).
    Each per-state term is in [0, 1] and the per-bin total is <= 1, so JSD is a
    bounded, directly comparable analog of the KL score. base-2 logs.

    Negative contributions (suppressed/depleted states) are clamped to 0.

    Bins where the total raw signal across all states is below min_signal are
    explicitly zeroed out before scoring. This prevents spurious positive scores
    in low-signal or repeat-masked regions caused by the pseudocount making the
    uniform distribution P = [1/N, ..., 1/N] appear to diverge from Q for focal
    tracks with Q_i < 1/N.

    Parameters
    ----------
    P                : (num_bins, num_states) float array  -- observed distribution
    Q                : (num_states,) float array           -- background distribution
    raw_signal_matrix: (num_bins, num_states) float array  -- pre-normalization signal
                       (after abs() for negative strands, before pseudocount).
                       If None, no low-signal masking is applied.
    min_signal       : float -- bins with total raw signal below this are zeroed
    normalize_score  : bool  -- divide KL by log2(num_states) (KL only; JSD is
                       already bounded, so this is ignored for JSD).
    score_method     : "kl" | "jsd"

    Returns
    -------
    scores : (num_bins, num_states) float array
    """
    Qb = Q[np.newaxis, :]
    if score_method == "jsd":
        M = 0.5 * (P + Qb)
        termP = P * np.log2((P + 1e-300) / (M + 1e-300))
        termQ = Qb * np.log2((Qb + 1e-300) / (M + 1e-300))
        scores = 0.5 * (termP + termQ)          # per-state JSD term (>= 0)
        scores = np.where(P > Qb, scores, 0.0)  # enrichment-only (match KL)
        scores = np.clip(scores, 0.0, None)     # guard tiny numerical negatives
    else:
        ratio = P / (Qb + 1e-300)
        scores = P * np.log2(ratio + 1e-300)
        scores = np.clip(scores, 0.0, None)

    # JSD is already bounded in [0,1]; the log2(N) rescale only applies to KL.
    do_norm = normalize_score and score_method != "jsd"

    # Zero out bins with insufficient total signal
    if raw_signal_matrix is not None:
        total_signal = raw_signal_matrix.sum(axis=1)  # (num_bins,)
        low_signal_mask = total_signal < min_signal
        scores[low_signal_mask, :] = 0.0
        n_masked = low_signal_mask.sum()
        if do_norm:
            num_states = scores.shape[1]
            norm_denom = math.log2(max(num_states, 2))
            scores = scores / norm_denom
        if n_masked > 0:
            return scores, int(n_masked)

    if do_norm:
        num_states = scores.shape[1]
        norm_denom = math.log2(max(num_states, 2))
        scores = scores / norm_denom

    return scores, 0


def compute_score_statistics(prob_cache, categories, metrics=["mean", "median", "p90"],
                             normalize_score=False, score_method="kl"):
    """
    Compute per-category quality statistics from the scored data in prob_cache.
    
    Parameters
    ----------
    prob_cache : dict { chrom: (bins, P, raw_clipped) }
                 from Pass 1; we'll re-score to match Pass 2 output
    categories : dict { state_key : [name, color] }
    metrics    : list of str (e.g., ["mean", "median", "p90"])
    
    Returns
    -------
    stats_dict : dict { category_name : { metric : value } }
    """
    num_states = len(categories)
    
    # Reconstruct Q from prob_cache
    Q_accumulator = np.zeros(num_states, dtype=np.float64)
    total_bins = 0
    for chrom in prob_cache:
        bins, P, raw_clipped = prob_cache[chrom]
        Q_accumulator += P.sum(axis=0)
        total_bins += len(bins)
    Q = Q_accumulator / total_bins
    
    # Collect all scores per category across all bins
    scores_by_category = {categories[str(i + 1)][0]: [] for i in range(num_states)}
    
    for chrom in prob_cache:
        bins, P, raw_clipped = prob_cache[chrom]
        scores, _ = kl_scores_per_bin(
            P,
            Q,
            raw_signal_matrix=raw_clipped,
            normalize_score=normalize_score,
            score_method=score_method,
        )
        
        # Accumulate scores for each category
        for state_idx in range(num_states):
            category_name = categories[str(state_idx + 1)][0]
            scores_by_category[category_name].extend(scores[:, state_idx].tolist())
    
    # Compute metrics for each category; also attach Q_i as 'q_background'
    stats_dict = {}
    name_to_q = {
        categories[str(i + 1)][0]: float(Q[i])
        for i in range(num_states)
    }
    for category_name, score_list in scores_by_category.items():
        if not score_list:
            stats_dict[category_name] = {"q_background": name_to_q.get(category_name, 0.0),
                                         **{m: 0.0 for m in metrics}}
            continue
        
        arr = np.array(score_list)
        category_stats = {"q_background": name_to_q.get(category_name, 0.0)}
        
        for metric in metrics:
            if metric == "mean":
                category_stats[metric] = float(np.mean(arr))
            elif metric == "median":
                category_stats[metric] = float(np.median(arr))
            elif metric == "std":
                category_stats[metric] = float(np.std(arr))
            elif metric == "min":
                category_stats[metric] = float(np.min(arr))
            elif metric == "max":
                category_stats[metric] = float(np.max(arr))
            elif metric.startswith("p") and metric[1:].isdigit():
                percentile = int(metric[1:])
                category_stats[metric] = float(np.percentile(arr, percentile))
            else:
                category_stats[metric] = np.nan
        
        stats_dict[category_name] = category_stats
    
    return stats_dict


def write_statistics_tsv(stats_dict, sample_name, metrics, out_path):
    """
    Write per-category statistics to a TSV file.
    
    Format: category, q_background, metric1, metric2, ...
    q_background is the genome-wide mean probability for each assay (sums to ~1).
    KL metrics are per-bin enrichment scores (mean/percentiles of score distribution).
    """
    import csv
    
    with open(out_path, "w", newline="") as f:
        fieldnames = ["category", "q_background"] + metrics
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        
        for category in sorted(stats_dict.keys()):
            row = {"category": category}
            row.update(stats_dict[category])
            writer.writerow(row)
    
    print(f"Statistics TSV:  {out_path}")


def plot_statistics(stats_dict, metrics, out_path):
    """
    Generate bar plots of quality statistics across categories.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib not available, skipping plots")
        return
    
    categories_sorted = sorted(stats_dict.keys())
    
    for metric in metrics:
        values = [stats_dict[cat][metric] for cat in categories_sorted]
        
        fig, ax = plt.subplots(figsize=(max(8, len(categories_sorted) * 0.6), 5))
        ax.bar(range(len(categories_sorted)), values, color="steelblue")
        ax.set_xticks(range(len(categories_sorted)))
        ax.set_xticklabels(categories_sorted, rotation=45, ha="right")
        ax.set_ylabel(f"Quality Score ({metric})")
        ax.set_title(f"Per-Assay Quality Statistics ({metric})")
        plt.tight_layout()
        
        plot_path = str(out_path).replace(".pdf", f"_{metric}.pdf")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        print(f"  Plot: {plot_path}")
        plt.close()


def plot_signal_distributions(raw_cache, norm_cache, categories, out_path,
                               sample_name="", max_bins=500_000):
    """
    Plot per-assay signal value distributions before and after quantile
    normalization, as a stacked grid of violin + box plots.

    One row per assay. Left column = raw (abs) signal; right column =
    normalized signal. Both columns share the same x-axis scale so shifts are
    immediately visible.

    Parameters
    ----------
    raw_cache  : dict {chrom: np.ndarray (n_bins, n_states)} raw abs signal
    norm_cache : dict {chrom: np.ndarray (n_bins, n_states)} normalized signal
                 If None, only the raw column is drawn.
    categories : list of (name, color) tuples
    out_path   : output PDF/PNG path
    sample_name: label shown in the figure title
    max_bins   : max bins sampled per assay to keep plotting fast
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("  WARNING: matplotlib not available, skipping signal distribution plots")
        return

    n_states = len(categories)
    has_norm = norm_cache is not None
    n_cols = 2 if has_norm else 1
    col_labels = ["Raw signal (abs)", "After normalization"]

    # Gather all per-state values across chroms, sub-sampled to max_bins
    rng = np.random.default_rng(42)

    def _collect(cache):
        cols = [[] for _ in range(n_states)]
        for arr in cache.values():
            for si in range(n_states):
                vals = arr[:, si].astype(np.float64)
                cols[si].append(vals)
        # Concatenate and sub-sample
        result = []
        for si in range(n_states):
            v = np.concatenate(cols[si]) if cols[si] else np.array([])
            if v.size > max_bins:
                v = rng.choice(v, size=max_bins, replace=False)
            result.append(v)
        return result

    raw_vals  = _collect(raw_cache)
    norm_vals = _collect(norm_cache) if has_norm else None

    # Determine global x-axis range (99th percentile of raw, ignoring zeros)
    all_raw_nonzero = np.concatenate(
        [v[v > 0] for v in raw_vals if v.size > 0]
    ) if any(v.size > 0 for v in raw_vals) else np.array([1.0])
    x_max_raw  = float(np.percentile(all_raw_nonzero, 99)) if all_raw_nonzero.size > 0 else 1.0
    x_max_norm = x_max_raw  # start with same; may adjust below
    if has_norm:
        all_norm_nonzero = np.concatenate(
            [v[v > 0] for v in norm_vals if v.size > 0]
        ) if any(v.size > 0 for v in norm_vals) else np.array([1.0])
        x_max_norm = float(np.percentile(all_norm_nonzero, 99)) if all_norm_nonzero.size > 0 else 1.0
    x_global = max(x_max_raw, x_max_norm)  # shared scale

    row_h  = 0.75
    fig_h  = max(4, n_states * row_h + 1.5)
    fig_w  = 5.0 * n_cols + 1.5
    fig, axes = plt.subplots(
        n_states, n_cols,
        figsize=(fig_w, fig_h),
        squeeze=False,
    )
    fig.subplots_adjust(hspace=0.12, wspace=0.08,
                        left=0.22, right=0.97, top=0.94, bottom=0.06)

    title = "Signal distributions"
    if sample_name:
        title += f" - {sample_name}"
    if has_norm:
        title += "  (raw vs. normalized)"
    fig.suptitle(title, fontsize=9, y=0.98)

    for si in range(n_states):
        name, color = categories[si]
        if color.lower() == "#ffffff":
            color = "#cccccc"

        for ci in range(n_cols):
            ax = axes[si, ci]
            vals = raw_vals[si] if ci == 0 else norm_vals[si]

            # Separate zero vs. non-zero for display
            nonzero = vals[vals > 0] if vals.size > 0 else np.array([])
            frac_nonzero = nonzero.size / max(vals.size, 1)

            if nonzero.size >= 5:
                # Violin plot of non-zero values
                vp = ax.violinplot([nonzero], positions=[0.5], vert=False,
                                   widths=0.7, showmedians=False,
                                   showextrema=False)
                for body in vp["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor("#555555")
                    body.set_alpha(0.6)
                    body.set_linewidth(0.5)

                # Overlay thin box showing IQR + median
                q25, q50, q75 = np.percentile(nonzero, [25, 50, 75])
                ax.plot([q25, q75], [0.5, 0.5], color="#222222",
                        lw=2.0, solid_capstyle="butt", zorder=3)
                ax.plot([q50], [0.5], "|", color="white",
                        ms=6, mew=1.5, zorder=4)
            else:
                ax.text(0.5, 0.5, "no signal",
                        transform=ax.transAxes,
                        fontsize=6, ha="center", va="center",
                        color="#999999", style="italic")

            # Shared x scale
            ax.set_xlim(0, x_global * 1.05)
            ax.set_ylim(0, 1)
            ax.set_yticks([])

            # Add % non-zero annotation
            ax.text(0.98, 0.85, f"{frac_nonzero*100:.0f}% > 0",
                    transform=ax.transAxes, fontsize=5,
                    ha="right", va="top", color="#444444")

            # Row label on leftmost column
            if ci == 0:
                ax.set_ylabel(name, fontsize=6.5, labelpad=3,
                              rotation=0, ha="right", va="center")

            # Column header on top row
            if si == 0:
                ax.set_title(col_labels[ci], fontsize=7.5, pad=4)

            # x-axis label on bottom row only
            if si == n_states - 1:
                ax.set_xlabel("Signal", fontsize=6)
                ax.tick_params(axis="x", labelsize=5)
            else:
                ax.set_xticks([])

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Signal distributions: {out_path}")


def print_normalization_qc(raw_cache, norm_cache, categories,
                           max_bins=500_000):
    """
    Print a compact per-assay normalization QC summary using sampled values.

    Reports, for each assay:
      - percent of bins with non-zero signal (%>0)
      - non-zero median
      - non-zero p90
    for both raw and normalized signal matrices.
    """
    if raw_cache is None or norm_cache is None:
        return

    n_states = len(categories)
    rng = np.random.default_rng(42)

    def _collect(cache):
        cols = [[] for _ in range(n_states)]
        total_counts = np.zeros(n_states, dtype=np.int64)
        nonzero_counts = np.zeros(n_states, dtype=np.int64)

        for arr in cache.values():
            for si in range(n_states):
                vals = arr[:, si].astype(np.float64)
                total_counts[si] += vals.size
                nz = vals[vals > 0]
                nonzero_counts[si] += nz.size
                if nz.size > 0:
                    cols[si].append(nz)

        sampled = []
        for si in range(n_states):
            v = np.concatenate(cols[si]) if cols[si] else np.array([], dtype=np.float64)
            if v.size > max_bins:
                v = rng.choice(v, size=max_bins, replace=False)
            sampled.append(v)

        return sampled, nonzero_counts, total_counts

    raw_vals, raw_nz, raw_tot = _collect(raw_cache)
    norm_vals, norm_nz, norm_tot = _collect(norm_cache)

    print("\nNormalization QC (sampled non-zero values):")
    print("  " +
          "Assay".ljust(18) +
          "Raw %>0".rjust(9) +
          "Raw med".rjust(10) +
          "Raw p90".rjust(10) +
          "Norm %>0".rjust(10) +
          "Norm med".rjust(10) +
          "Norm p90".rjust(10))

    for si in range(n_states):
        name = categories[si][0]
        raw_pct = (100.0 * raw_nz[si] / raw_tot[si]) if raw_tot[si] > 0 else 0.0
        norm_pct = (100.0 * norm_nz[si] / norm_tot[si]) if norm_tot[si] > 0 else 0.0

        raw_med = float(np.median(raw_vals[si])) if raw_vals[si].size > 0 else 0.0
        raw_p90 = float(np.percentile(raw_vals[si], 90)) if raw_vals[si].size > 0 else 0.0
        norm_med = float(np.median(norm_vals[si])) if norm_vals[si].size > 0 else 0.0
        norm_p90 = float(np.percentile(norm_vals[si], 90)) if norm_vals[si].size > 0 else 0.0

        print("  " +
              f"{name[:18]:<18}" +
              f"{raw_pct:>9.1f}" +
              f"{raw_med:>10.4g}" +
              f"{raw_p90:>10.4g}" +
              f"{norm_pct:>10.1f}" +
              f"{norm_med:>10.4g}" +
              f"{norm_p90:>10.4g}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def check_bigwig_paths(bw_paths):
    """
    Test that every BigWig path (local file or remote URL) can be opened by
    pyBigWig before the main run starts.  Prints a summary and returns the
    list of failed paths.  Does NOT call sys.exit() so the caller can decide
    whether to abort or continue.
    """
    try:
        import pyBigWig
    except ImportError:
        return []  # can't check without pyBigWig; run() will handle the error

    failed = []
    for p in bw_paths:
        p_str = str(p)
        try:
            bw = pyBigWig.open(p_str)
            if bw is None:
                raise RuntimeError("open returned None")
            bw.close()
        except Exception as exc:
            failed.append((p_str, str(exc)))

    if failed:
        print(f"\nERROR: {len(failed)} BigWig file(s) could not be opened:")
        for path, reason in failed:
            print(f"  MISSING/UNREADABLE: {path}")
            if reason:
                print(f"                      ({reason})")
    else:
        print(f"  Pre-flight check: all {len(bw_paths)} BigWig file(s) accessible.")
    return failed


def run(bw_paths, out_path, chrom_sizes, chroms=None, regions=None,
    normalize_tracks=False, normalize_method="nonzero-quantile",
    min_signal=MIN_SIGNAL, stats_metrics=None, stats_plots=False,
    signal_plots=False, categories=None, negative_strand_states=None,
    jobs=1, summary_chrom=None, skip_preflight=False,
    normalize_score=False, blacklist=None, min_signal_per_track=None,
    min_signal_percentile=None, floors_tsv=None, write_floors_tsv=None,
    sample_name=None, cohort_ref=None, score_method="kl", bins_bed=None):
    try:
        import pyBigWig
    except ImportError:
        sys.exit("ERROR: pyBigWig is required.  Install with:  pip3 install pyBigWig")
    try:
        import pysam
    except ImportError:
        sys.exit("ERROR: pysam is required.  Install with:  pip3 install pysam")

    _cats = categories if categories is not None else ALL_CATEGORIES
    _neg_strand = negative_strand_states if negative_strand_states is not None else NEGATIVE_STRAND_STATES

    # Resolve sample-name for floors-tsv operations
    _sample_name = sample_name
    if _sample_name is None and out_path is not None:
        _sample_name = Path(str(out_path)).name
        for suf in (".qcat.bgz", ".qcat.gz", ".qcat", ".bgz", ".gz"):
            if _sample_name.endswith(suf):
                _sample_name = _sample_name[:-len(suf)]
                break

    # Resolve per-track noise floors. Three modes (in precedence order):
    #   1. --floors-tsv (load precomputed floors, skip percentile pass)
    #   2. --min-signal-percentile (data-driven, sampling pass)
    #   3. --min-signal-per-track (absolute floors, manual)
    _min_signal_per_track = parse_min_signal_per_track(min_signal_per_track, _cats)
    if floors_tsv is not None:
        _tsv_floors = load_floors_from_tsv(
            floors_tsv, _sample_name or "", _cats
        )
        if _tsv_floors:
            if _min_signal_per_track:
                print(f"  (floors-tsv overrides --min-signal-per-track for "
                      f"sample '{_sample_name}')", flush=True)
            _min_signal_per_track = _tsv_floors
            print(f"Loaded per-track floors for sample '{_sample_name}' from "
                  f"{floors_tsv}:", flush=True)
            for k, v in sorted(_tsv_floors.items()):
                print(f"  {k:<20} = {v:.4g}", flush=True)
        else:
            print(f"  WARNING: --floors-tsv '{floors_tsv}' has no matching "
                  f"rows for sample '{_sample_name}'. Falling back to other "
                  f"floor sources.", file=sys.stderr)
    if _min_signal_per_track:
        active = ", ".join(f"{k}={v}" for k, v in _min_signal_per_track.items())
        print(f"Per-track minimum-signal floors (absolute): {active}")
    total_pertrack_zeroed = 0

    # Pre-flight: verify all BigWig files are readable before starting
    if not skip_preflight:
        print("Pre-flight: checking BigWig accessibility...")
        failed = check_bigwig_paths(bw_paths)
        if failed:
            sys.exit(
                f"\nAborted: fix the {len(failed)} missing/unreadable BigWig file(s) above "
                f"before re-running."
            )

    num_states = len(bw_paths)
    if not (MIN_STATES <= num_states <= MAX_STATES):
        sys.exit(
            f"ERROR: Expected between {MIN_STATES} and {MAX_STATES} BigWig files, "
            f"got {num_states}."
        )

    # Determine chromosomes to process
    if chroms is None:
        chroms = sorted(
            chrom_sizes.keys(),
            key=lambda c: (
                0 if c.startswith("chr") else 1,
                int(c.replace("chr", "").replace("X", "23").replace("Y", "24"))
                if c.replace("chr", "").replace("X", "23").replace("Y", "24").isdigit()
                else 99,
                c,
            ),
        )
    if len(chroms) == 0:
        sys.exit(
            "ERROR: no chromosomes selected for processing. "
            "Check --regions content and chromosome names."
        )

    # An explicit bin grid (adaptive segmentation) defines the chromosome
    # universe: restrict processing to chroms present in the BED (e.g. main
    # chromosomes only) so contigs absent from the grid are skipped instead of
    # raising a KeyError on chrom_bin_counts downstream.
    if bins_bed is not None:
        bed_chroms = set(bins_bed.keys())
        dropped = [c for c in chroms if c not in bed_chroms]
        chroms = [c for c in chroms if c in bed_chroms]
        if dropped:
            _ex = ", ".join(dropped[:5]) + ("..." if len(dropped) > 5 else "")
            print(f"  --bins-bed covers {len(bed_chroms)} chroms; skipping "
                  f"{len(dropped)} chrom(s) absent from the grid ({_ex}).")
        if len(chroms) == 0:
            sys.exit(
                "ERROR: none of the requested chromosomes are present in "
                "--bins-bed. Check chromosome naming between the grid and "
                "--chrom-sizes/--regions."
            )

    print(f"Processing {len(chroms)} chromosomes: {', '.join(chroms)}")

    # Data-driven per-track floors: compute Nth percentile of non-zero
    # values per track, once chromosome list is finalized. Overrides any
    # absolute floors from --min-signal-per-track. Skipped if floors were
    # already loaded via --floors-tsv (those have higher precedence).
    _signal_stats = {}
    _percentile_used = None
    _floors_from_percentile = False
    if (min_signal_percentile is not None and float(min_signal_percentile) > 0
            and not (floors_tsv is not None and _min_signal_per_track)):
        print(f"\nComputing per-track noise floors as the "
              f"{float(min_signal_percentile):g}th percentile of non-zero "
              f"values per track...", flush=True)
        _percentile_floors, _signal_stats = compute_percentile_floors_per_track(
            bw_paths=bw_paths,
            categories=_cats[:num_states],
            percentile=min_signal_percentile,
            chrom_sizes=chrom_sizes,
            chroms=chroms,
            regions=regions,
            neg_strand_states=_neg_strand,
            blacklist=blacklist,
            collect_stats=True,
        )
        _percentile_used = float(min_signal_percentile)
        if _percentile_floors:
            if _min_signal_per_track:
                print("  (percentile floors override absolute floors)")
            _min_signal_per_track = _percentile_floors
            _floors_from_percentile = True

    # Write out the computed floors for reuse in later runs (e.g. permutations).
    if write_floors_tsv is not None and _min_signal_per_track:
        _out_floors_path = str(write_floors_tsv)
        write_floors_tsv_file(
            _min_signal_per_track,
            _sample_name or "",
            _out_floors_path,
            percentile=_percentile_used,
        )
        print(f"Wrote per-track floors -> {_out_floors_path}",
              flush=True)

    # Determine parallelism (CPU cores) for per-chromosome work
    jobs = int(jobs) if jobs is not None else 1
    if jobs < 1:
        jobs = 1
    jobs = min(jobs, len(chroms))
    print(f"Using {jobs} parallel worker(s)")

    if jobs == 1:
        print(f"Opening {num_states} BigWig files...")
        bw_handles = [pyBigWig.open(str(p)) for p in bw_paths]

    # ---- Pass 1: build P matrices and accumulate Q --------------------------
    prob_cache         = {}
    _raw_signal_cache  = {}  # populated only when signal_plots=True
    _norm_signal_cache = {}  # populated only when signal_plots=True
    Q_accumulator = np.zeros(num_states, dtype=np.float64)
    total_bins    = 0
    Q_denom       = 0.0   # bins (no-bed) or sum of widths (bed) -> width-weighted Q

    tmp_path = str(out_path) + ".tmp.tsv"

    if jobs == 1:
        for chrom in chroms:
            chrom_len = chrom_sizes[chrom]
            if bins_bed is not None:
                bins = bins_bed.get(chrom, [])
                if not bins:
                    continue
            elif regions is not None:
                chrom_regions = regions.get(chrom, [])
                if not chrom_regions:
                    continue
                bins = bins_for_regions(chrom_regions)
            else:
                bins = bins_for_chrom(chrom_len)
            n = len(bins)

            signal_matrix = np.zeros((n, num_states), dtype=np.float32)
            for si, bw in enumerate(bw_handles):
                try:
                    signal_matrix[:, si] = mean_signal_in_bins(bw, chrom, bins)
                except Exception as e:
                    print(f"  WARNING: could not read {chrom} from BigWig {si+1}: {e}")

            for neg_idx in _neg_strand:
                col = neg_idx - 1  # convert to 0-based
                if col < num_states:
                    signal_matrix[:, col] = np.abs(signal_matrix[:, col])

            if blacklist is not None:
                bl_mask = bins_overlapping_blacklist(bins, blacklist, chrom)
                n_blacklisted = int(bl_mask.sum())
                if n_blacklisted > 0:
                    signal_matrix[bl_mask, :] = 0.0
                    print(
                        f"  {chrom}: zeroed {n_blacklisted:,} blacklisted bins",
                        file=sys.stderr,
                    )

            raw_abs = np.clip(signal_matrix, 0.0, None).astype(np.float32)

            # Per-track noise-floor masking (see apply_per_track_noise_floor).
            n_pertrack = 0
            if _min_signal_per_track:
                n_pertrack = apply_per_track_noise_floor(
                    signal_matrix, _cats, _min_signal_per_track
                )
            if n_pertrack:
                total_pertrack_zeroed += n_pertrack

            if normalize_tracks:
                signal_matrix = normalize_signal_matrix(
                    signal_matrix,
                    method=normalize_method,
                    cohort_ref=cohort_ref,
                )

            P = signals_to_prob(signal_matrix)
            raw_clipped = np.clip(signal_matrix, 0.0, None).astype(np.float32)
            # carry pre-normalization raw_abs so the min_signal floor gates true
            # signal (controlled bin set across normalized/un-normalized runs)
            prob_cache[chrom] = (bins, P, raw_clipped, raw_abs)
            if signal_plots:
                _raw_signal_cache[chrom]  = raw_abs
                _norm_signal_cache[chrom] = raw_clipped
            if bins_bed is not None:
                widths = np.array([e - s for (s, e) in bins], dtype=np.float64)
                Q_accumulator += (P * widths[:, None]).sum(axis=0)
                Q_denom += float(widths.sum())
            else:
                Q_accumulator += P.sum(axis=0)
                Q_denom += n
            total_bins += n
            print(f"  {chrom}: {n:,} bins")

        Q = Q_accumulator / Q_denom
        cats_dict = build_categories(num_states, _cats)
        print(f"\nBackground Q computed over {total_bins:,} bins.")
        print("  Top states by Q: " + ", ".join(
            f"{_cats[i][0]}={Q[i]:.4f}"
            for i in np.argsort(Q)[::-1][:5]
        ))

        for bw in bw_handles:
            bw.close()

        # ---- Pass 2: score and write qcat rows ----------------------------------
        print(f"\nWriting scored qcat rows -> {tmp_path}")

        rows_written = 0
        total_masked = 0
        # accumulate KL totals for summary
        chrom_score_totals = {}
        global_score_total = 0.0
        chrom_state_totals = {}       # chrom -> array of per-state sums
        global_state_totals = np.zeros(num_states, dtype=np.float64)
        with open(tmp_path, "w") as fh:
            bin_id = 1
            for chrom in prob_cache:
                bins, P, raw_clipped, raw_abs = prob_cache[chrom]
                scores, n_masked = kl_scores_per_bin(
                    P, Q,
                    raw_signal_matrix=raw_abs,
                    min_signal=min_signal,
                    normalize_score=normalize_score,
                    score_method=score_method,
                )
                total_masked += n_masked
                # accumulate chromosome-level and global KL sums
                chrom_sum = float(scores.sum())
                chrom_score_totals[chrom] = chrom_sum
                global_score_total += chrom_sum
                # accumulate per-state totals
                chrom_state_totals[chrom] = scores.sum(axis=0)
                global_state_totals += scores.sum(axis=0)
                for i, (s, e) in enumerate(bins):
                    row_scores = scores[i]
                    # qcat format: [score, state_index] pairs sorted descending by score
                    pairs = sorted(
                        ([float(f"{sc:.6g}"), state_idx]
                         for state_idx, sc in enumerate(row_scores, start=1)),
                        key=lambda x: -x[0],
                    )
                    raw_vec = raw_clipped[i].tolist()
                    qcat_col = "id:" + str(bin_id) + ",qcat:" + _json.dumps(pairs, separators=(",", ":"))
                    fh.write(chrom + "\t" + str(s) + "\t" + str(e) + "\t" + qcat_col + "\n")
                    rows_written += 1
                    bin_id += 1
    else:
        # Parallel per-chromosome processing using temporary cache files
        with tempfile.TemporaryDirectory() as tmp_dir:
            tasks = []
            for chrom in chroms:
                chrom_len = chrom_sizes[chrom]
                regions_for_chrom = regions.get(chrom, []) if regions is not None else None
                bins_override = bins_bed.get(chrom, []) if bins_bed is not None else None
                if bins_bed is not None and not bins_override:
                    continue   # no consensus bins on this chrom
                tasks.append((chrom, chrom_len, bw_paths, normalize_tracks,
                              normalize_method, _neg_strand,
                              regions_for_chrom, tmp_dir, blacklist,
                              _cats, _min_signal_per_track, cohort_ref,
                              bins_override))

            ctx = multiprocessing.get_context("spawn")
            with ctx.Pool(jobs) as pool:
                results = pool.starmap(_compute_chrom_cache, tasks)

            # Keep results in chrom order
            chrom_order = {c: i for i, c in enumerate(chroms)}
            results.sort(key=lambda x: chrom_order[x[0]])

            chrom_npz_paths = {}
            chrom_bin_counts = {}
            Q_denom = 0.0
            for chrom, n, P_sum, q_denom, npz_path, n_pertrack_zeroed in results:
                chrom_npz_paths[chrom] = npz_path
                chrom_bin_counts[chrom] = n
                Q_accumulator += P_sum
                Q_denom += q_denom
                total_bins += n
                total_pertrack_zeroed += int(n_pertrack_zeroed)
                print(f"  {chrom}: {n:,} bins")

            Q = Q_accumulator / Q_denom
            cats_dict = build_categories(num_states, _cats)
            print(f"\nBackground Q computed over {total_bins:,} bins.")
            print("  Top states by Q: " + ", ".join(
                f"{_cats[i][0]}={Q[i]:.4f}"
                for i in np.argsort(Q)[::-1][:5]
            ))

            # ---- Pass 2: score and write qcat rows ----------------------------------
            print(f"\nWriting scored qcat rows -> {tmp_path}")

            score_tasks = []
            start_id = 1
            for chrom in chroms:
                n_bins = chrom_bin_counts[chrom]
                out_tmp = os.path.join(tmp_dir, f"{chrom}.qcat.tmp")
                score_tasks.append((
                    chrom,
                    Q,
                    chrom_npz_paths[chrom],
                    min_signal,
                    normalize_score,
                    start_id,
                    out_tmp,
                    score_method,
                ))
                start_id += n_bins

            with ctx.Pool(jobs) as pool:
                score_results = pool.starmap(_score_chrom_from_cache, score_tasks)

            # Aggregate results and merge per-chrom qcat files
            rows_written = 0
            total_masked = 0
            chrom_score_totals = {}
            global_score_total = 0.0
            chrom_state_totals = {}
            global_state_totals = np.zeros(num_states, dtype=np.float64)

            with open(tmp_path, "w") as out_f:
                for chrom in chroms:
                    out_file = os.path.join(tmp_dir, f"{chrom}.qcat.tmp")
                    with open(out_file) as in_f:
                        out_f.write(in_f.read())

            for chrom, n_masked, chrom_sum, state_totals, n_bins, _ in score_results:
                rows_written += n_bins
                total_masked += n_masked
                chrom_score_totals[chrom] = chrom_sum
                global_score_total += chrom_sum
                chrom_state_totals[chrom] = state_totals
                global_state_totals += state_totals

            # Load prob_cache for stats/plots if requested
            if stats_metrics or signal_plots:
                for chrom in chroms:
                    data = np.load(chrom_npz_paths[chrom], allow_pickle=False)
                    bins = [tuple(map(int, b)) for b in data["bins"]]
                    raw_clipped = data["raw_clipped"]
                    prob_cache[chrom] = (bins, data["P"], raw_clipped)
                    if signal_plots:
                        _raw_signal_cache[chrom] = data["raw_abs"]
                        _norm_signal_cache[chrom] = raw_clipped

    print(f"  {rows_written:,} rows written.")
    print(f"  {total_masked:,} bins zeroed (total raw signal < {min_signal})")
    if _min_signal_per_track:
        print(f"  {total_pertrack_zeroed:,} per-track cells zeroed "
              f"(per-track floors: "
              f"{', '.join(f'{k}={v}' for k,v in _min_signal_per_track.items())})")

    # report KL score summary
    print("\nKL score summary:")
    print(f"  Total over all processed bins: {global_score_total:.2f}")
    # if only a single chromosome processed, note it is per-chrom and global are same
    for chrom, tot in chrom_score_totals.items():
        print(f"    {chrom}: {tot:.2f}")

    # per-state breakdown for the chromosome(s) of interest and global
    if chrom_state_totals:
        print("\nPer-state KL sums (background vs signal contributions):")
        def fmt_array(arr):
            return ", ".join(f"{_cats[i][0]}={arr[i]:.2f}"
                             for i in range(len(arr)))
        if summary_chrom and summary_chrom in chrom_state_totals:
            print(f"  {summary_chrom}: {fmt_array(chrom_state_totals[summary_chrom])}")
        elif len(chrom_state_totals) == 1:
            chrom, arr = next(iter(chrom_state_totals.items()))
            print(f"  {chrom}: {fmt_array(arr)}")
        else:
            for chrom, arr in chrom_state_totals.items():
                print(f"  {chrom}: {fmt_array(arr)}")
        print(f"  GLOBAL: {fmt_array(global_state_totals)}")

    # ---- Sort ---------------------------------------------------------------
    sorted_path = tmp_path + ".sorted"
    print("Sorting...")
    result = subprocess.run(
        ["sort", "-k1,1", "-k2,2n", tmp_path],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        with open(sorted_path, "w") as f:
            f.write(result.stdout)
    else:
        print("  (Falling back to Python sort)")
        with open(tmp_path) as f:
            lines = f.readlines()
        lines.sort(key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
        with open(sorted_path, "w") as f:
            f.writelines(lines)
    os.remove(tmp_path)

    # ---- Bgzip + tabix ------------------------------------------------------
    out_str = str(out_path)
    print(f"Bgzipping -> {out_str}")
    pysam.tabix_compress(sorted_path, out_str, force=True)
    os.remove(sorted_path)

    print(f"Indexing  -> {out_str}.tbi")
    pysam.tabix_index(out_str, preset="bed", force=True)

    # ---- Statistics (optional) ----------------------------------------------
    _kl_stats = None
    if stats_metrics:
        print(f"\nComputing quality statistics...")
        _kl_stats = compute_score_statistics(
            prob_cache,
            cats_dict,
            metrics=stats_metrics,
            normalize_score=normalize_score,
            score_method=score_method,
        )
        
        # Write TSV (includes q_background column automatically)
        base = str(out_path).replace(".qcat.bgz", "").replace(".bgz", "")
        stats_tsv_path = base + "_stats.tsv"
        write_statistics_tsv(_kl_stats, base, stats_metrics, stats_tsv_path)
        
        # Generate plots if requested
        if stats_plots:
            print("Generating statistics plots...")
            stats_plot_path = base + "_stats.pdf"
            plot_statistics(_kl_stats, stats_metrics, stats_plot_path)

    # ---- Signal distribution QC (optional) ----------------------------------
    if signal_plots:
        base = str(out_path).replace(".qcat.bgz", "").replace(".bgz", "")
        sample_label = Path(base).name
        dist_path = base + "_signal_distributions.pdf"
        print("\nGenerating signal distribution plots...")
        if normalize_tracks:
            print_normalization_qc(
                raw_cache=_raw_signal_cache,
                norm_cache=_norm_signal_cache,
                categories=_cats[:num_states],
            )
        plot_signal_distributions(
            raw_cache=_raw_signal_cache,
            norm_cache=_norm_signal_cache if normalize_tracks else None,
            categories=_cats[:num_states],
            out_path=dist_path,
            sample_name=sample_label,
        )

    print(f"\nDone!")
    print(f"  qcat file:   {out_str}")
    print(f"  tabix index: {out_str}.tbi")
    return out_str, cats_dict, prob_cache, Q, _kl_stats, _signal_stats, _percentile_used


# ---------------------------------------------------------------------------
# Companion files
# ---------------------------------------------------------------------------

def write_categories_json(path, categories, q_values=None, kl_stats=None,
                          normalize_score=False):
    """
    Write categories JSON.  Optional fields:
      q_values / q_values_by_assay : background probability fractions (sum ~1).
      kl_stats_by_assay            : per-assay KL score statistics (from --stats).
    """
    payload = {"categories": categories}
    if q_values is not None:
        payload["q_values"] = {
            str(i + 1): float(q_values[i])
            for i in range(min(len(q_values), len(categories)))
        }
        payload["q_values_by_assay"] = {
            categories[str(i + 1)][0]: float(q_values[i])
            for i in range(min(len(q_values), len(categories)))
        }
    if kl_stats is not None:
        # kl_stats is the stats_dict from compute_score_statistics
        # Strip q_background from here (already in q_values_by_assay) to avoid
        # duplication; keep only the KL score metrics.
        payload["kl_stats_by_assay"] = {
            assay: {k: v for k, v in metrics.items() if k != "q_background"}
            for assay, metrics in kl_stats.items()
        }
    payload["normalize_score"] = bool(normalize_score)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Categories JSON: {path}")


def write_debug_bigwigs(prob_cache, chrom_sizes, bw_paths, categories, out_dir,
                        normalize_tracks, normalize_method="nonzero-quantile"):
    """
    Write per-track BigWig files capturing two intermediate signal states:

      1. abs_signal  -- raw signal after abs() for negative-strand tracks,
                        before quantile normalization or pseudocount.
                        One file per track: <out_dir>/<name>_abs.bw

    2. norm_signal -- signal after normalization (only written
                        when normalize_tracks=True), before pseudocount.
                        One file per track: <out_dir>/<name>_norm.bw

    These are written at 200 bp bin resolution matching the scoring pipeline.
    Negative-strand tracks store the absolute value (magnitude), because that
    is what enters the KL computation.
    """
    import pyBigWig as pbw

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    num_states = len(bw_paths)
    chroms_in_cache = list(prob_cache.keys())

    # prob_cache[chrom] = (bins, P, raw_clipped)
    # raw_clipped is the pre-pseudocount signal matrix used for scoring.
    # When normalize_tracks=True, raw_clipped is already post-normalization.

    debug_paths = []

    for si in range(num_states):
        state_key = str(si + 1)
        name  = categories[state_key][0]
        safe  = name.replace(" ", "_").replace("/", "_").replace("+", "plus").replace("-", "minus")

        # --- abs signal bigwig ---
        abs_path = str(out_dir / (safe + "_abs.bw"))
        bw_out = pbw.open(abs_path, "w")
        # Build header: chrom sizes for chroms we have data for
        header = [(c, chrom_sizes[c]) for c in chroms_in_cache if c in chrom_sizes]
        bw_out.addHeader(header)
        for chrom in chroms_in_cache:
            bins, P, raw_clipped = prob_cache[chrom]
            col = raw_clipped[:, si].tolist()
            starts = [b[0] for b in bins]
            ends   = [b[1] for b in bins]
            chroms_list = [chrom] * len(bins)
            # pyBigWig addEntries expects lists
            bw_out.addEntries(chroms_list, starts, ends=ends, values=col)
        bw_out.close()
        print(f"  Debug abs signal: {abs_path}")
        debug_paths.append(("abs", si, abs_path))

        # --- norm signal bigwig (only if normalize_tracks) ---
        if normalize_tracks:
            norm_path = str(out_dir / (safe + "_norm.bw"))
            bw_out = pbw.open(norm_path, "w")
            bw_out.addHeader(header)
            for chrom in chroms_in_cache:
                bins, P, raw_clipped = prob_cache[chrom]
                # Avoid re-normalizing here; this matrix is already normalized.
                col = raw_clipped[:, si].tolist()
                starts = [b[0] for b in bins]
                ends   = [b[1] for b in bins]
                chroms_list = [chrom] * len(bins)
                bw_out.addEntries(chroms_list, starts, ends=ends, values=col)
            bw_out.close()
            print(f"  Debug norm signal: {norm_path}")
            debug_paths.append(("norm", si, norm_path))

    return debug_paths


def format_region_for_title(region_str):
    """
    Convert a region string from "chr6:40793981-41688054" to readable form "chr6:40.8-41.7Mb".
    """
    if not region_str or ':' not in region_str:
        return region_str
    try:
        chrom, coords = region_str.split(':')
        if '-' not in coords:
            return region_str
        start_s, end_s = coords.split('-')
        start, end = int(start_s), int(end_s)
        # Convert to Mb with one decimal place
        start_mb = f"{start / 1e6:.1f}"
        end_mb = f"{end / 1e6:.1f}"
        return f"{chrom}:{start_mb}-{end_mb}Mb"
    except (ValueError, AttributeError):
        return region_str


def extract_sample_name(qcat_path):
    """
    Extract sample name from qcat filename (e.g., "S3T3.qcat.bgz" -> "S3T3").
    """
    try:
        name = Path(qcat_path).name
        if '.qcat' in name:
            return name.split('.qcat')[0]
        return name
    except:
        return ""


def estimate_bigwig_abs_max(bw_path):
    """
    Estimate the maximum absolute signal present in a BigWig.

    Returns None if the file cannot be read. Uses per-chromosome max stats,
    which is fast and sufficient for display scaling decisions.
    """
    try:
        import pyBigWig
    except ImportError:
        return None

    try:
        with pyBigWig.open(str(bw_path)) as bw:
            chrom_sizes = bw.chroms()
            if not chrom_sizes:
                return None
            max_abs = 0.0
            for chrom, clen in chrom_sizes.items():
                v = bw.stats(chrom, 0, clen, type="max")[0]
                if v is None or math.isnan(v):
                    continue
                max_abs = max(max_abs, abs(float(v)))
            return max_abs
    except Exception:
        return None


def write_ini(qcat_path, cats_json_path, ini_path, bw_paths=None, categories=None,
              bed_files=None, vlines_bed=None, debug_bw_paths=None,
              sample_name=None, region_str=None, negative_strand_states=None,
              normalize_score=False):
    """
    Write a pyGenomeTracks .ini file containing:
      - epilogos stacked bar track
      - epilogos inverted track
      - one bigwig fill track per input file, colored to match its category
      - optional debug tracks (abs signal and/or quantile-normalized signal)
      - optional BED annotation tracks
      - optional vertical lines from a BED file
      - x-axis

    bw_paths     : list of Path objects for the input BigWig files (optional)
    categories   : dict returned by build_categories() (optional)
    bed_files    : list of (path, title) tuples for BED annotation tracks (optional)
    vlines_bed   : path to BED file used for vertical lines (optional)
    debug_bw_paths: list of (kind, state_idx, path) tuples from write_debug_bigwigs()
    sample_name  : optional sample name to include in track titles
    region_str   : optional region string (e.g., "chr6:40800000-41700000") to include in track titles
    """
    lines = []
    lines.append("# pyGenomeTracks configuration -- generated by bigwig_to_qcat.py")
    lines.append("# Edit file paths if you move the output files.")
    lines.append("#")
    lines.append(f"# normalize_score = {str(bool(normalize_score)).lower()}")
    lines.append("#")
    lines.append("# Bin size recommendations for different region scales:")
    lines.append("#   <500 kb:      200 bp (current default; shows fine detail)")
    lines.append("#   500 kb-5 Mb:  1-2 kb (typical whole-gene resolution)")
    lines.append("#   5-50 Mb:      5-10 kb (TAD-scale view)")
    lines.append("#   >50 Mb:       20-50 kb (very broad chromosomal overview)")
    lines.append("#")
    lines.append("# Usage:")
    lines.append(f"#   pyGenomeTracks --tracks {ini_path} --region chr1:10000000-11000000 -o out.png")
    lines.append("")

    # --- epilogos tracks ---
    lines.append("[epilogos]")
    lines.append(f"file = {qcat_path}")
    lines.append(f"categories_file = {cats_json_path}")
    lines.append("height = 5")
    lines.append("title = Epilogos (mm10)")
    lines.append("")
    lines.append("[spacer]")
    lines.append("height = 0.3")
    lines.append("")
    lines.append("[epilogos inverted]")
    lines.append(f"file = {qcat_path}")
    lines.append(f"categories_file = {cats_json_path}")
    lines.append("height = 5")
    lines.append("title = Epilogos inverted")
    lines.append("orientation = inverted")
    lines.append("")
    lines.append("[spacer]")
    lines.append("height = 0.5")
    lines.append("")

    # --- bigwig tracks ---
    if bw_paths and categories:
        lines.append("# --- Input BigWig tracks (colored to match epilogos categories) ---")
        lines.append("# For low-signal assays (e.g., RNA-seq, sparse focal-TF peaks),")
        lines.append(f"# positive-strand max_value is set to {MIN_SIGNAL} to avoid over-amplification of noise.")
        lines.append("")
        
        # Extract sample name if not provided
        sample_label = sample_name or extract_sample_name(qcat_path)
        region_label = format_region_for_title(region_str) if region_str else ""
        context_label = ""
        if sample_label or region_label:
            parts = []
            if sample_label:
                parts.append(sample_label)
            if region_label:
                parts.append(f"[{region_label}]")
            context_label = " | ".join(parts) + " | "
        
        for i, bw_path in enumerate(bw_paths):
            state_key = str(i + 1)
            name  = categories[state_key][0]
            color = categories[state_key][1]
            # Quiescent/Low is white -- use light gray so it is visible
            if color.lower() == "#ffffff":
                color = "#cccccc"
            safe_name = name.replace(" ", "_").replace("/", "_")
            _nss = negative_strand_states if negative_strand_states is not None else NEGATIVE_STRAND_STATES
            is_neg_strand = (i + 1) in _nss
            abs_max = estimate_bigwig_abs_max(bw_path)
            is_low_signal = (abs_max is not None) and (abs_max < MIN_SIGNAL)
            
            lines.append(f"[bigwig {safe_name}]")
            lines.append(f"file = {bw_path}")
            # Title with sample + region context
            title = f"{context_label}{name}"
            lines.append(f"title = {title}")
            lines.append(f"color = {color}")
            lines.append("height = 2")
            
            # Data-driven axis scaling based on track signal range
            if is_neg_strand:
                # Negative strand tracks are always shown from auto lower bound to 0.
                lines.append("min_value = auto")
                lines.append("max_value = 0")
            elif is_low_signal:
                # Weak positive track: keep a fixed baseline window
                lines.append("min_value = 0")
                lines.append(f"max_value = {MIN_SIGNAL}")
            else:
                # Non-weak positive track: allow dynamic range
                lines.append("min_value = 0")
                lines.append("max_value = auto")
            
            lines.append("type = fill")
            lines.append("nans_to_zeros = true")
            lines.append("show_data_range = true")
            lines.append("")
            lines.append("[spacer]")
            lines.append("height = 0.1")
            lines.append("")

    # --- Debug BigWig tracks (abs signal and quantile-normalized signal) ---
    if debug_bw_paths and categories:
        kind_labels = {"abs": "abs signal (pre-norm)", "norm": "quantile-normalized"}
        # Group by kind so we emit all abs tracks together, then all norm tracks
        for kind in ("abs", "norm"):
            kind_entries = [(si, p) for (k, si, p) in debug_bw_paths if k == kind]
            if not kind_entries:
                continue
            label = kind_labels[kind]
            lines.append(f"# --- Debug: {label} ---")
            lines.append("")
            for si, dbw_path in kind_entries:
                state_key = str(si + 1)
                name  = categories[state_key][0]
                color = categories[state_key][1]
                if color.lower() == "#ffffff":
                    color = "#cccccc"
                safe_name = name.replace(" ", "_").replace("/", "_").replace("+", "plus").replace("-", "minus")
                abs_max = estimate_bigwig_abs_max(dbw_path)
                is_low_signal = (abs_max is not None) and (abs_max < MIN_SIGNAL)
                lines.append(f"[bigwig {safe_name}_{kind}]")
                lines.append(f"file = {dbw_path}")
                lines.append(f"title = {name} ({label})")
                lines.append(f"color = {color}")
                lines.append("height = 1")        # smaller than raw tracks
                lines.append("alpha = 0.5")       # semi-transparent
                # Debug tracks store non-negative magnitudes; clamp weak tracks.
                if is_low_signal:
                    lines.append("min_value = 0")
                    lines.append(f"max_value = {MIN_SIGNAL}")
                else:
                    lines.append("min_value = 0")
                    lines.append("max_value = auto")
                lines.append("type = fill")
                lines.append("nans_to_zeros = true")
                lines.append("show_data_range = true")
                lines.append("")
                lines.append("[spacer]")
                lines.append("height = 0.05")
                lines.append("")

    # --- BED annotation tracks ---
    if bed_files:
        lines.append("# --- BED annotation tracks ---")
        lines.append("")
        for bed_path, bed_title in bed_files:
            safe = str(bed_path).replace(" ", "_")
            lines.append(f"[bed {bed_title}]")
            lines.append(f"file = {bed_path}")
            lines.append(f"title = {bed_title}")
            lines.append("height = 3")
            lines.append("fontsize = 10")
            lines.append("file_type = bed")
            lines.append("gene_rows = 5")
            lines.append("color = #1f78b4")
            lines.append("border_color = black")
            lines.append("display = stacked")
            lines.append("style = UCSC")
            lines.append("show_labels = true")
            lines.append("")
            lines.append("[spacer]")
            lines.append("height = 0.2")
            lines.append("")

    # --- Vertical lines from BED ---
    if vlines_bed:
        lines.append("# --- Vertical lines (positions taken from BED file) ---")
        lines.append("[vlines]")
        lines.append(f"file = {vlines_bed}")
        lines.append("type = vlines")
        lines.append("line_width = 1.5")
        lines.append("color = #e41a1c")
        lines.append("alpha = 0.7")
        lines.append("")

    lines.append("[x-axis]")
    lines.append("")

    content = "\n".join(lines)
    with open(ini_path, "w") as f:
        f.write(content)
    print(f"Tracks INI:      {ini_path}")
    print(f"\nTo plot, run:")
    print(f"  pyGenomeTracks --tracks {ini_path} --region chr1:10000000-11000000 -o out.png")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="BigWig files -> 200 bp KL-scored qcat.bgz for pyGenomeTracks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument(
        "--bw", nargs="+", metavar="FILE",
        help=(f"Between {MIN_STATES} and {MAX_STATES} BigWig files in category order. "
              "The first N categories from the hardcoded list will be used."),
    )
    grp.add_argument(
        "--bw-dir", metavar="DIR",
        help=("Directory of BigWig files (*.bw or *.bigwig). "
              "Files are sorted alphabetically -- prefix with 01_, 02_ etc. "
              "to control order."),
    )
    grp.add_argument(
        "--sheet", metavar="FILE",
        help=(
            "TSV sample sheet describing one or more qcat jobs. Columns:\n"
            "    bw           comma-separated list of BigWig paths (or use bw1,bw2,... columns)\n"
            "    out          output qcat filename (optional; defaults to <sample>.qcat.bgz)\n"
            "    sample/name  optional label used for progress messages\n"
            "Additional columns such as chroms, bed, vlines etc. are ignored."
        ),
    )

    parser.add_argument(
        "--out", required=True, metavar="FILE",
        help="Output path, e.g. experiment.qcat.bgz",
    )
    parser.add_argument(
        "--genome", default="mm10", choices=["mm10"],
        help="Genome assembly (default: mm10). "
             "Use --chrom-sizes to supply a custom assembly.",
    )
    parser.add_argument(
        "--chrom-sizes", metavar="FILE",
        help="Two-column chrom.sizes file. Overrides --genome built-in sizes.",
    )
    parser.add_argument(
        "--chroms", nargs="+", metavar="CHR",
        help="Process only these chromosomes, e.g. --chroms chr1 chr2 chrX",
    )
    parser.add_argument(
        "--regions", metavar="BED",
          help=("Regions file for testing mode. Accepts either standard BED "
              "(chrom start end) or TSV with a 'region' column containing "
              "chr:start-end (e.g. regions_template.tsv). Only 200 bp bins "
              "overlapping these intervals are scored. Background Q is "
              "computed from these bins only, so scores reflect local "
              "rather than genome-wide background."),
    )
    parser.add_argument(
        "--blacklist", metavar="BED", default=None,
        help=(
            "BED file of regions to exclude from scoring. Any 200 bp bin "
            "overlapping a blacklist interval has its signal set to zero "
            "before normalisation, P computation, and Q estimation. "
            "Accepts standard 3+ column BED."
        ),
    )
    parser.add_argument(
        "--detect-unmappable", action="store_true",
        help=(
            "Scan all sheet rows jointly to identify bins with near-zero "
            "signal across most track x sample combinations. Writes a BED "
            "file to --unmappable-out (or an automatic default). "
            "Detection-only by default; re-run with --blacklist to apply."
        ),
    )
    parser.add_argument(
        "--unmappable-out", metavar="BED", default=None,
        help=(
            "Output BED path for --detect-unmappable. Default: "
            "<out_prefix>_detected_blacklist.bed"
        ),
    )
    parser.add_argument(
        "--unmappable-zero-frac", type=float, default=0.90, metavar="FLOAT",
        help=(
            "Fraction of track x sample combinations that must be near-zero "
            "to flag a bin as unmappable (default: 0.90)."
        ),
    )
    parser.add_argument(
        "--unmappable-min-samples", type=int, default=None, metavar="N",
        help=(
            "Minimum number of samples that must be near-zero across all "
            "their tracks for a bin to be flagged. Default: all samples."
        ),
    )
    parser.add_argument(
        "--feature-set",
        choices=["promoters", "gene_bodies", "all_genes", "tss_2kb", "ctcf_sites"],
        help=(
            "Score only bins overlapping a predefined mm10 feature set. "
            "Equivalent to passing a bundled BED file to --regions. "
            "Ignored if --regions is also specified (--regions takes precedence). "
            "Available sets: promoters (TSS +/-1kb), gene_bodies (full gene spans), "
            "all_genes (union of above), tss_2kb (TSS +/-2kb), ctcf_sites "
            "(high-confidence mm10 CTCF binding sites from ENCODE)."
        ),
    )
    parser.add_argument(
        "--normalize-tracks", action="store_true",
        help=(
            "Normalize each BigWig track to the same global signal "
            "distribution before computing P and Q. By default this uses "
            "nonzero-only quantile normalization, which preserves low "
            "non-zero values better than the legacy method. This prevents "
            "broad tracks (e.g. H3K27ac, RNAseq) from dominating the "
            "background Q and ensures focal tracks (e.g. CTCF) are "
            "scored on a level playing field."
        ),
    )
    parser.add_argument(
        "--normalize-score", action="store_true",
        help=(
            "Divide all per-bin KL scores by log2(N) so the total score S(b) "
            "is bounded in [0, 1] regardless of the number of tracks. "
            "Scores are then interpretable as a fraction of the maximum "
            "possible divergence for this experiment."
        ),
    )
    parser.add_argument(
        "--score-method", choices=["kl", "jsd"], default="kl",
        help=(
            "Per-bin scoring divergence. 'kl' (default): clamped per-state KL "
            "contribution (unbounded above). 'jsd': bounded Jensen-Shannon "
            "contribution, per-state in [0,1] and per-bin total <= 1, using the "
            "same enrichment-only convention as KL -- a bounded analog of the "
            "KL score that is more robust to focal/outlier tracks."
        ),
    )
    parser.add_argument(
        "--bins-bed", default=None, metavar="BED",
        help=(
            "Opt-in adaptive (variable-width) binning. Score on the consensus "
            "segmentation in this BED (chrom start end ...) instead of the fixed "
            "grid; the background Q is computed width-weighted. Build the BED with "
            "build_adaptive_segmentation.py. When omitted, fixed binning is used "
            "(bit-identical to the default)."
        ),
    )
    parser.add_argument(
        "--normalize-method",
        choices=["nonzero-quantile", "quantile", "cohort-quantile"],
        default="nonzero-quantile",
        help=(
            "Normalization method used with --normalize-tracks. Default: "
            "nonzero-quantile (within-sample, across this sample's tracks). "
            "Use quantile for the legacy within-sample method, or "
            "cohort-quantile to map each track onto a COHORT-WIDE per-track "
            "reference (requires --cohort-reference), which equalizes "
            "between-sample dynamic-range differences for the same assay."
        ),
    )
    parser.add_argument(
        "--cohort-reference", default=None, metavar="NPZ",
        help=(
            "Per-track cohort reference (.npz from "
            "build_cohort_quantile_reference.py). Required for "
            "--normalize-method cohort-quantile."
        ),
    )
    parser.add_argument(
        "--bed", nargs="+", metavar="FILE",
        help=(
            "One or more BED annotation files to include as tracks in the "
            "generated .ini file (e.g. cbe_mm10.bed AgRgenes_mm10_s.bed). "
            "Each file is added as a stacked UCSC-style gene track below "
            "the BigWig tracks."
        ),
    )
    parser.add_argument(
        "--vlines", metavar="FILE",
        help=(
            "BED file whose interval start positions are drawn as vertical "
            "lines across all tracks (e.g. cbe_mm10.bed). "
            "Uses pyGenomeTracks [vlines] type."
        ),
    )
    parser.add_argument(
        "--min-signal", type=float, default=MIN_SIGNAL, metavar="FLOAT",
        help=(
            f"Bins where the total raw signal across all tracks is below "
            f"this threshold are zeroed out before KL scoring, preventing "
            f"spurious scores in low-signal or repeat-masked regions. "
            f"Default: {MIN_SIGNAL}. Set to 0 to disable."
        ),
    )
    parser.add_argument(
        "--pseudocount", type=float, default=PSEUDOCOUNT, metavar="FLOAT",
        help=(
            f"Zero-clamp epsilon added to the raw signal vector before "
            f"normalization to a probability distribution P (clip negatives, "
            f"add pseudocount, divide by row sum). Regularizes the per-bin "
            f"composition; recovery is expected to be insensitive to this over "
            f"a wide range. Default: {PSEUDOCOUNT}."
        ),
    )
    parser.add_argument(
        "--min-signal-per-track", type=str, default=None, metavar="SPEC",
        help=(
            "Optional per-track noise-floor masking applied BEFORE the "
            "compositional Q vector is computed. Tracks with sparse/noisy "
            "signal (e.g. RNA-seq) can contribute noise-floor values to Q "
            "when they share a bin with denser tracks (ATAC, ChIP-seq) "
            "that pass --min-signal. This flag zeros out per-bin entries "
            "whose individual track signal is below its per-track floor. "
            "Accepts either the literal 'default' (uses built-in floors: "
            "RNAseq +/-=0.1, ATAC/CTCF/Cohesin/H3K27ac=0.05) or a "
            "comma-separated list of NAME=VALUE pairs "
            "(e.g. 'RNAseq +=0.1,RNAseq -=0.1,ATAC=0.05'). "
            "Track names must match the active category list. "
            "Default: disabled. NOTE: prefer --min-signal-percentile for "
            "data-driven calibration."
        ),
    )
    parser.add_argument(
        "--min-signal-percentile", type=float, default=None, metavar="PCT",
        help=(
            "Data-driven per-track noise-floor masking. Computes the Nth "
            "percentile of NON-ZERO values for each track across the "
            "processed regions, then zeros out per-bin entries below their "
            "track's percentile floor before the compositional Q vector "
            "is computed. Recommended for mixed-normalization datasets "
            "(e.g. CPM RNA-seq + RPGC ChIP-seq) because the percentile "
            "rule is invariant to absolute scaling. Typical value: 5.0 "
            "(masks the bottom 5%% of non-zero signal per track). "
            "Computed floors are printed to stderr and applied uniformly "
            "across all chromosomes. If both --min-signal-per-track and "
            "--min-signal-percentile are given, --min-signal-percentile "
            "takes precedence. Default: disabled."
        ),
    )
    parser.add_argument(
        "--floors-tsv", type=str, default=None, metavar="PATH",
        help=(
            "Load pre-computed per-track noise floors from a TSV file. "
            "Use --sample-name to select which sample's rows to use. The "
            "TSV must have columns: sample, track, floor. Floors loaded "
            "this way SKIP the --min-signal-percentile sampling pass "
            "entirely (~4 min saved per run). Recommended for permutation "
            "null runs: compute floors once on the real BigWigs with "
            "--min-signal-percentile and --write-floors-tsv, then reuse "
            "those floors for all permutation rounds via --floors-tsv. "
            "Percentile-of-nonzero floors are invariant under circular "
            "shift so this reuse is statistically sound."
        ),
    )
    parser.add_argument(
        "--write-floors-tsv", type=str, default=None, metavar="PATH",
        help=(
            "After computing per-track floors via --min-signal-percentile, "
            "write them to a TSV at PATH for later reuse via --floors-tsv. "
            "Format: sample <TAB> track <TAB> floor <TAB> percentile "
            "<TAB> source."
        ),
    )
    parser.add_argument(
        "--sample-name", type=str, default=None, metavar="NAME",
        help=(
            "Sample name to associate with this run. Required when using "
            "--floors-tsv (selects which sample's floors to load) or "
            "--write-floors-tsv (used as the sample column value). If "
            "omitted, the output filename stem is used."
        ),
    )
    parser.add_argument(
        "--jobs", type=int, default=1, metavar="N",
        help=(
            "Number of parallel workers (chromosome-level) to use. "
            "Set to 1 for serial execution (default)."
        ),
    )
    parser.add_argument(
        "--debug-tracks", action="store_true",
        help=(
            "Write intermediate per-track BigWig files for debugging: "
            "one set showing abs() signal (after negative-strand correction, "
            "before normalization) and, if --normalize-tracks is active, one "
            "set showing quantile-normalized signal. Files are written to a "
            "'debug/' subdirectory next to the output file and added to the "
            ".ini as smaller, semi-transparent tracks below the raw BigWigs."
        ),
    )
    parser.add_argument(
        "--stats", nargs="+", metavar="METRIC",
        help=(
            "Compute per-assay quality statistics and write to *_stats.tsv. "
            "Metrics: mean, median, std, min, max, p10, p25, p75, p90. "
            "Example: --stats mean median p90"
        ),
    )
    parser.add_argument(
        "--stats-plots", action="store_true",
        help="Generate bar plots of statistics (requires --stats).",
    )
    parser.add_argument(
        "--signal-plots", action="store_true",
        help=(
            "Plot per-assay signal value distributions before and after "
            "quantile normalization (one PDF with violin plots). "
            "Both stages use the same x-axis scale so shifts are visible. "
            "Output: <out_base>_signal_distributions.pdf"
        ),
    )
    parser.add_argument(
        "--categories", metavar="YAML",
        help=(
            "YAML file defining category names, colors, and strand orientation. "
            "Overrides the built-in 15-state mm10 list. "
            "Example files are in categories/ (mm10_15state.yaml, "
            "mm10_5state_minimal.yaml, hg38_15state.yaml)."
        ),
    )
    parser.add_argument(
        "--no-extras", action="store_true",
        help="Skip writing the categories JSON and .ini files.",
    )

    args = parser.parse_args()

    # Opt-in adaptive binning: load the consensus segmentation once and pass it
    # to run(); when None, fixed binning is used (bit-identical to default).
    bins_bed = None
    if args.bins_bed:
        bins_bed = load_bins_bed(args.bins_bed)
        _nb = sum(len(v) for v in bins_bed.values())
        print(f"  Adaptive binning: {_nb:,} consensus bins from {args.bins_bed} "
              f"({len(bins_bed)} chroms); Q will be width-weighted.")

    # Load the cohort-wide per-track reference once if requested. It is a
    # (n_states, L) array keyed by track-column order, mapped onto each track.
    cohort_ref = None
    if args.normalize_method == "cohort-quantile":
        if not args.cohort_reference:
            sys.exit("ERROR: --normalize-method cohort-quantile requires "
                     "--cohort-reference <npz>.")
        _cref = np.load(args.cohort_reference, allow_pickle=True)
        cohort_ref = np.asarray(_cref["ref"], dtype=np.float64)
        print(f"  Cohort reference: {args.cohort_reference} "
              f"({cohort_ref.shape[0]} tracks x {cohort_ref.shape[1]} quantile points)")
    elif args.cohort_reference:
        print("  NOTE: --cohort-reference given but --normalize-method is not "
              "cohort-quantile; ignoring the reference.")

    # module-level PSEUDOCOUNT as a global at call time, so setting it here
    # (before run()) propagates to all scoring without rethreading it through
    # every call site.
    if args.pseudocount != PSEUDOCOUNT:
        globals()["PSEUDOCOUNT"] = args.pseudocount
        print(f"  Pseudocount (zero-clamp epsilon): {args.pseudocount}")

    if not (0.0 < args.unmappable_zero_frac <= 1.0):
        sys.exit("ERROR: --unmappable-zero-frac must be in (0, 1].")
    if args.unmappable_min_samples is not None and args.unmappable_min_samples < 1:
        sys.exit("ERROR: --unmappable-min-samples must be >= 1.")

    feature_set_tmp = None
    if args.feature_set and not args.regions:
        entries = BUILTIN_FEATURE_SETS[args.feature_set]
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
        with tmp as fh:
            for chrom, start, end in entries:
                fh.write(f"{chrom}\t{start}\t{end}\n")
        feature_set_tmp = tmp.name
        args.regions = feature_set_tmp
        print(f"Using built-in feature set: {args.feature_set} ({len(entries)} regions)")

    # Load categories from YAML if provided; otherwise use built-in defaults.
    _cli_cats = None
    _cli_neg_strand = None
    if args.categories:
        _cli_cats, _cli_neg_strand = load_categories_yaml(args.categories)
        print(f"Categories loaded from: {args.categories}  ({len(_cli_cats)} states)")

    # Helper to parse a sheet describing bigwig jobs
    def load_bw_sheet(path):
        samples = []
        with open(path) as fh:
            header = None
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if header is None:
                    header = [h.strip().lower() for h in parts]
                    required = {"bw"}
                    missing = required - set(header)
                    if missing:
                        sys.exit("ERROR: Sheet is missing column(s): " + ",".join(sorted(missing)))
                    continue
                row = dict(zip(header, [p.strip() for p in parts]))
                sample_label = row.get("sample") or row.get("name") or ""
                out_path = row.get("out")
                if not out_path:
                    if sample_label:
                        out_path = sample_label + ".qcat.bgz"
                    else:
                        sys.exit("ERROR: each sheet row must provide either 'out' or 'sample' to name output")
                # parse bigwig list from bw column or bw1,bw2,...
                bw_list = []
                if row.get("bw"):
                    bw_list = [p.strip() for p in row["bw"].split(",") if p.strip()]
                else:
                    # collect bw1,bw2,.. until missing
                    idx = 1
                    while True:
                        key = f"bw{idx}"
                        if key not in row:
                            break
                        if row[key]:
                            bw_list.append(row[key])
                        idx += 1
                if not bw_list:
                    sys.exit(f"ERROR: no bigwig paths found for sample {sample_label}")
                samples.append({
                    "sample": sample_label,
                    "out": out_path,
                    "bw_paths": [Path(p) for p in bw_list],
                })
        return samples

    # Resolve BigWig paths (or sheet of jobs)
    sheet_jobs = None
    if args.sheet:
        sheet_jobs = load_bw_sheet(args.sheet)
    elif args.bw_dir:
        bw_dir = Path(args.bw_dir)
        bw_paths = sorted(
            list(bw_dir.glob("*.bw")) + list(bw_dir.glob("*.bigwig"))
        )
        if not bw_paths:
            sys.exit(f"ERROR: No .bw / .bigwig files found in {bw_dir}")
        print(f"Found {len(bw_paths)} BigWig files in {bw_dir}:")
        for i, p in enumerate(bw_paths, 1):
            cat = ALL_CATEGORIES[i - 1][0] if i <= MAX_STATES else "?"
            print(f"  {i:>2}. {p.name}  ->  {cat}")
    else:
        bw_paths = [Path(p) for p in args.bw]

    # Resolve output path
    out_path = Path(args.out)
    if not str(out_path).endswith(".bgz"):
        out_path = Path(str(out_path) + ".bgz")

    # Resolve chromosome sizes
    chrom_sizes = load_chrom_sizes(args.chrom_sizes) if args.chrom_sizes else MM10_CHROM_SIZES

    # Filter to requested chromosomes
    chroms = args.chroms if args.chroms else None

    # Load regions BED if provided (testing mode)
    regions = None
    summary_chrom = None
    if args.regions:
        regions = load_regions_bed(args.regions)
        region_chroms = sorted(regions.keys())
        chroms = [c for c in chroms if c in regions] if chroms else region_chroms
        if len(region_chroms) == 1:
            summary_chrom = region_chroms[0]
        elif chroms and len(chroms) == 1:
            summary_chrom = chroms[0]
        total_intervals = sum(len(v) for v in regions.values())
        total_bp = sum(e - s for ivs in regions.values() for s, e in ivs)
        print(f"Testing mode: {args.regions}")
        print(f"  {total_intervals} intervals across {len(region_chroms)} chroms "
              f"({total_bp:,} bp covered)")

    # Detection-only pass support (requires sheet mode).
    detected_blacklist = None
    if args.detect_unmappable:
        if sheet_jobs is None:
            sys.exit("ERROR: --detect-unmappable requires --sheet")
        min_s = (
            args.unmappable_min_samples
            if args.unmappable_min_samples is not None
            else len(sheet_jobs)
        )
        detected_blacklist = detect_unmappable_bins(
            sheet_jobs,
            chrom_sizes,
            regions=regions,
            zero_frac=args.unmappable_zero_frac,
            min_samples=min_s,
            bin_size=BIN_SIZE,
            workers=args.jobs,
        )
        out_prefix = str(out_path).replace(".qcat.bgz", "").replace(".bgz", "")
        bl_out = args.unmappable_out or (out_prefix + "_detected_blacklist.bed")
        write_blacklist_bed(
            detected_blacklist,
            chrom_sizes,
            bl_out,
            bin_size=BIN_SIZE,
            regions=regions,
        )
        print(f"\nDetected blacklist written to: {bl_out}", file=sys.stderr)
        if regions:
            print(
                "WARNING: this was a region-restricted smoke test only; re-run without --regions for a genome-wide detected blacklist.",
                file=sys.stderr,
            )
        print("Re-run with --blacklist to apply this mask to scoring.\n", file=sys.stderr)

    user_blacklist = None
    if args.blacklist:
        user_blacklist = load_blacklist(args.blacklist)
        n_chroms = len(user_blacklist)
        n_intervals = sum(len(v) for v in user_blacklist.values())
        print(
            f"Loaded blacklist: {n_intervals:,} intervals on {n_chroms} chromosomes",
            file=sys.stderr,
        )

    combined_blacklist = user_blacklist
    if user_blacklist is not None and detected_blacklist is not None:
        combined_blacklist = {}
        all_chroms = set(user_blacklist.keys()) | set(detected_blacklist.keys())
        for chrom in all_chroms:
            intervals = []
            intervals.extend(user_blacklist.get(chrom, []))
            mask = detected_blacklist.get(chrom)
            if mask is not None:
                if regions and chrom in regions:
                    region_bins = bins_for_regions(regions[chrom])
                    for bi, flagged in enumerate(mask):
                        if flagged and bi < len(region_bins):
                            intervals.append(region_bins[bi])
                else:
                    for bi, flagged in enumerate(mask):
                        if flagged:
                            s = bi * BIN_SIZE
                            e = min(s + BIN_SIZE, chrom_sizes.get(chrom, s + BIN_SIZE))
                            intervals.append((s, e))
            intervals.sort()
            combined_blacklist[chrom] = intervals
        print(
            "Merging detected unmappable bins with user-provided blacklist for scoring.",
            file=sys.stderr,
        )

    # If sheet_jobs defined, process each row separately and exit
    if sheet_jobs is not None:
        # --- Global pre-flight: check every BigWig across all sheet rows ---
        print(f"\nPre-flight: checking all BigWig files across {len(sheet_jobs)} sheet row(s)...")
        all_failed = []
        for job in sheet_jobs:
            job_label = job.get("sample") or job.get("out")
            job_failed = check_bigwig_paths(job["bw_paths"])
            for path, reason in job_failed:
                all_failed.append((job_label, path, reason))
        if all_failed:
            print(f"\nERROR: {len(all_failed)} missing/unreadable BigWig file(s) found across the sheet:")
            for job_label, path, reason in all_failed:
                print(f"  [{job_label}]  {path}")
                if reason:
                    print(f"               ({reason})")
            sys.exit(
                f"\nAborted: fix all {len(all_failed)} missing/unreadable file(s) above "
                f"before re-running. No samples were processed."
            )
        else:
            print(f"  Pre-flight OK: all BigWig files accessible across all {len(sheet_jobs)} row(s).\n")

        for job in sheet_jobs:
            sample_label = job.get("sample") or job.get("out")
            bw_paths_job = job["bw_paths"]
            out_path_job = Path(job["out"])
            if not str(out_path_job).endswith(".bgz"):
                out_path_job = Path(str(out_path_job) + ".bgz")

            print(f"\n{'='*60}")
            _banner_cats = _cli_cats if _cli_cats is not None else ALL_CATEGORIES
            print(f"  BigWig -> qcat pipeline  (sheet row: {sample_label})")
            print(f"  States:  {len(bw_paths_job)} (of {MAX_STATES} possible)")
            print(f"  Using:   {', '.join(_banner_cats[i][0] for i in range(len(bw_paths_job)))}")
            print(f"  Genome:  {args.genome}")
            print(f"  Bin:     {BIN_SIZE} bp")
            print(f"  Output:  {out_path_job}")
            if args.regions:
                print(f"  Regions: {args.regions}  [testing mode]")
            if args.normalize_tracks:
                print(f"  Normalize: enabled ({args.normalize_method})")
            if args.debug_tracks:
                print(f"  Debug tracks: enabled")
            print(f"  Min signal threshold: {args.min_signal}")
            if args.min_signal_per_track:
                print(f"  Per-track noise floors: {args.min_signal_per_track}")
            if args.min_signal_percentile:
                print(f"  Per-track noise-floor percentile: p{args.min_signal_percentile:g}")
            if args.floors_tsv:
                print(f"  Per-track floors from TSV: {args.floors_tsv}"
                      + (f"  (sample={args.sample_name})" if args.sample_name else ""))
            if args.write_floors_tsv:
                print(f"  Will write floors TSV -> {args.write_floors_tsv}")
            if args.bed:
                print(f"  BED tracks: {', '.join(args.bed)}")
            if args.vlines:
                print(f"  Vlines:    {args.vlines}")
            print(f"{'='*60}\n")

            out_str, categories, prob_cache, q_values, kl_stats, signal_stats, signal_pct = run(
                bw_paths_job, out_path_job, chrom_sizes,
                chroms=chroms, regions=regions,
                normalize_tracks=args.normalize_tracks,
                normalize_method=args.normalize_method,
                cohort_ref=cohort_ref,
                score_method=args.score_method,
                normalize_score=args.normalize_score,
                min_signal=args.min_signal,
                min_signal_per_track=args.min_signal_per_track,
                min_signal_percentile=args.min_signal_percentile,
                floors_tsv=args.floors_tsv,
                write_floors_tsv=args.write_floors_tsv,
                sample_name=args.sample_name,
                stats_metrics=args.stats,
                stats_plots=args.stats_plots,
                signal_plots=args.signal_plots,
                categories=_cli_cats,
                negative_strand_states=_cli_neg_strand,
                jobs=args.jobs,
                summary_chrom=summary_chrom,
                skip_preflight=True,
                blacklist=combined_blacklist,
                bins_bed=bins_bed,
            )

            if not args.no_extras:
                base = str(out_path_job).replace(".qcat.bgz", "").replace(".bgz", "")
                cats_json = base + "_cats.json"
                ini_file  = base + "_tracks.ini"
                write_categories_json(
                    cats_json,
                    categories,
                    q_values=q_values,
                    normalize_score=args.normalize_score,
                )
                # Signal-distribution stats from the percentile-floor pass.
                if signal_stats:
                    sigstats_tsv = base + "_signal_stats.tsv"
                    sample_label = Path(base).name
                    write_signal_stats_tsv(
                        signal_stats, sample_label, sigstats_tsv,
                        percentile=signal_pct,
                    )
                    print(f"  Signal-distribution stats -> {sigstats_tsv}")
                    if args.stats_plots or args.signal_plots:
                        sigstats_pdf = base + "_signal_stats.pdf"
                        plot_signal_stats(
                            signal_stats, sigstats_pdf,
                            sample_name=sample_label,
                            percentile=signal_pct,
                        )
                        print(f"  Signal-distribution plot  -> {sigstats_pdf}")

                # Write debug BigWig files if requested
                debug_bw_paths = None
                if args.debug_tracks:
                    debug_dir = Path(base).parent / "debug"
                    print("Writing debug BigWig files -> " + str(debug_dir) + "/")
                    debug_bw_paths = write_debug_bigwigs(
                        prob_cache, MM10_CHROM_SIZES if not args.chrom_sizes
                            else load_chrom_sizes(args.chrom_sizes),
                        bw_paths_job, categories, debug_dir,
                        normalize_tracks=args.normalize_tracks,
                        normalize_method=args.normalize_method,
                    )

                # Build bed_files list: [(path, title), ...]
                bed_files = None
                if args.bed:
                    bed_files = []
                    for bp in args.bed:
                        title = Path(bp).stem.replace("_", " ").replace("-", " ")
                        bed_files.append((bp, title))

                vlines_bed = args.vlines if args.vlines else None

                write_ini(out_str, cats_json, ini_file,
                          bw_paths=bw_paths_job, categories=categories,
                          bed_files=bed_files, vlines_bed=vlines_bed,
                          debug_bw_paths=debug_bw_paths, sample_name=sample_label,
                          negative_strand_states=_cli_neg_strand,
                          normalize_score=args.normalize_score)
        if feature_set_tmp is not None:
            Path(feature_set_tmp).unlink(missing_ok=True)
        return

    # Startup banner
    _banner_cats = _cli_cats if _cli_cats is not None else ALL_CATEGORIES
    print(f"\n{'='*60}")
    print(f"  BigWig -> qcat pipeline")
    print(f"  States:  {len(bw_paths)} (of {MAX_STATES} possible)")
    print(f"  Using:   {', '.join(_banner_cats[i][0] for i in range(len(bw_paths)))}")
    print(f"  Genome:  {args.genome}")
    print(f"  Bin:     {BIN_SIZE} bp")
    print(f"  Output:  {out_path}")
    if args.regions:
        print(f"  Regions: {args.regions}  [testing mode]")
    if args.normalize_tracks:
        print(f"  Normalize: enabled ({args.normalize_method})")
    if args.debug_tracks:
        print(f"  Debug tracks: enabled")
    print(f"  Min signal threshold: {args.min_signal}")
    if args.min_signal_per_track:
        print(f"  Per-track noise floors: {args.min_signal_per_track}")
    if args.min_signal_percentile:
        print(f"  Per-track noise-floor percentile: p{args.min_signal_percentile:g}")
    if args.floors_tsv:
        print(f"  Per-track floors from TSV: {args.floors_tsv}"
              + (f"  (sample={args.sample_name})" if args.sample_name else ""))
    if args.write_floors_tsv:
        print(f"  Will write floors TSV -> {args.write_floors_tsv}")
    if args.bed:
        print(f"  BED tracks: {', '.join(args.bed)}")
    if args.vlines:
        print(f"  Vlines:    {args.vlines}")
    print(f"{'='*60}\n")

    out_str, categories, prob_cache, q_values, kl_stats, signal_stats, signal_pct = run(bw_paths, out_path, chrom_sizes,
                              chroms=chroms, regions=regions,
                              normalize_tracks=args.normalize_tracks,
                              normalize_method=args.normalize_method,
                              cohort_ref=cohort_ref,
                              score_method=args.score_method,
                              normalize_score=args.normalize_score,
                              min_signal=args.min_signal,
                              min_signal_per_track=args.min_signal_per_track,
                              min_signal_percentile=args.min_signal_percentile,
                              floors_tsv=args.floors_tsv,
                              write_floors_tsv=args.write_floors_tsv,
                              sample_name=args.sample_name,
                              stats_metrics=args.stats,
                              stats_plots=args.stats_plots,
                              signal_plots=args.signal_plots,
                              categories=_cli_cats,
                              negative_strand_states=_cli_neg_strand,
                              jobs=args.jobs,
                              summary_chrom=summary_chrom,
                              blacklist=combined_blacklist,
                              bins_bed=bins_bed)

    if not args.no_extras:
        base = str(out_path).replace(".qcat.bgz", "").replace(".bgz", "")
        cats_json = base + "_cats.json"
        ini_file  = base + "_tracks.ini"
        write_categories_json(
            cats_json,
            categories,
            q_values=q_values,
            kl_stats=kl_stats,
            normalize_score=args.normalize_score,
        )
        # Signal-distribution stats from the percentile-floor pass.
        if signal_stats:
            sigstats_tsv = base + "_signal_stats.tsv"
            sample_label = Path(base).name
            write_signal_stats_tsv(
                signal_stats, sample_label, sigstats_tsv,
                percentile=signal_pct,
            )
            print(f"  Signal-distribution stats -> {sigstats_tsv}")
            if args.stats_plots or args.signal_plots:
                sigstats_pdf = base + "_signal_stats.pdf"
                plot_signal_stats(
                    signal_stats, sigstats_pdf,
                    sample_name=sample_label,
                    percentile=signal_pct,
                )
                print(f"  Signal-distribution plot  -> {sigstats_pdf}")

        # Write debug BigWig files if requested
        debug_bw_paths = None
        if args.debug_tracks:
            debug_dir = Path(base).parent / "debug"
            print("Writing debug BigWig files -> " + str(debug_dir) + "/")
            debug_bw_paths = write_debug_bigwigs(
                prob_cache, MM10_CHROM_SIZES if not args.chrom_sizes
                    else load_chrom_sizes(args.chrom_sizes),
                bw_paths, categories, debug_dir,
                normalize_tracks=args.normalize_tracks,
                normalize_method=args.normalize_method,
            )

        # Build bed_files list: [(path, title), ...]
        bed_files = None
        if args.bed:
            bed_files = []
            for bp in args.bed:
                title = Path(bp).stem.replace("_", " ").replace("-", " ")
                bed_files.append((bp, title))

        vlines_bed = args.vlines if args.vlines else None

        # Extract sample name from output path for track title context
        sample_label = extract_sample_name(out_str)

        write_ini(out_str, cats_json, ini_file,
                  bw_paths=bw_paths, categories=categories,
                  bed_files=bed_files, vlines_bed=vlines_bed,
                  debug_bw_paths=debug_bw_paths, sample_name=sample_label,
                  negative_strand_states=_cli_neg_strand,
                  normalize_score=args.normalize_score)

    if feature_set_tmp is not None:
        Path(feature_set_tmp).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
