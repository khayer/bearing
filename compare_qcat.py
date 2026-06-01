#!/usr/bin/env python3
"""
compare_qcat.py
===============
Compare multiple qcat.bgz epilogos files across conditions and replicates.

Reads a sample sheet (TSV), parses the per-bin probability distributions
stored inside qcat.bgz files, and produces the following outputs:

  1. Pairwise Jensen-Shannon divergence heatmap  (*_jsd_heatmap.pdf)
  2. Per-state Spearman rho matrix               (*_spearman.pdf)
  3. PCA plot colored by condition               (*_pca.pdf)
  4. Differential epilogos qcat.bgz files        (*_diff_*.qcat.bgz)
     one per pair of conditions (mean of replicates vs mean of replicates),
     showing signed per-bin per-state KL divergence.
  5. Cross-sample stats comparison bar plot       (stats_comparison.pdf)
      if *_stats.tsv files are found beside sample qcat files.
  6. pyGenomeTracks comparison INI                (comparison_tracks.ini)
      with configurable qcat score limits.
  7. Per-region matplotlib comparison plots       (regions/*.pdf)
      one PDF per region in --regions-file; all tracks share fixed y-limits
      for direct cross-sample and cross-region comparability.
  8. Per-region JSD heatmaps                      (regions/*_jsd.pdf)
      one PDF per region in --regions-file; computed from the same aligned
      region score matrices used for PCA and Spearman.
  9. Per-sample Q-vector comparability heatmap     (*_q_pair_jsd_heatmap.pdf, *.png)
      with Q-pair JSD matrix and summary TSV (*_q_pair_jsd.tsv).

SAMPLE SHEET FORMAT (TSV, header required)
------------------------------------------
    sample          condition   replicate   qcat
    ctrl_rep1       ctrl        1           /path/to/ctrl_rep1.qcat.bgz
    ctrl_rep2       ctrl        2           /path/to/ctrl_rep2.qcat.bgz
    treat_rep1      treat       1           /path/to/treat_rep1.qcat.bgz
    treat_rep2      treat       2           /path/to/treat_rep2.qcat.bgz

  - sample     : unique label used in plots
  - condition  : groups replicates together (e.g. "ctrl", "KO", "treat")
  - replicate  : integer replicate number within the condition
  - qcat       : path to the qcat.bgz file produced by bigwig_to_qcat.py

USAGE
-----
    python compare_qcat.py --sheet samples.tsv --out comparison/

    # Only run specific analyses:
    python compare_qcat.py --sheet samples.tsv --out comparison/ --skip-diff

    # Restrict to specific chromosomes (faster for testing):
    python compare_qcat.py --sheet samples.tsv --out comparison/ --chroms chr1 chr2

    # Set qcat display range in generated INI (default 5.0):
    # epilogos tracks do not support min_value/max_value, so compare_qcat.py
    # writes clipped plotting copies of qcat files used by the INI:
    #   sample tracks: [0, qcat_max]
    #   diff tracks:   [-qcat_max, +qcat_max]
    python compare_qcat.py --sheet samples.tsv --out comparison/ --qcat-max 5

    # Add an optional gene BED/GTF track to the generated INI and region plots:
    python compare_qcat.py --sheet samples.tsv --out comparison/ --genes cbe_mm10.bed

    # Plot specific regions as matplotlib figures (one PDF per region):
    python compare_qcat.py --sheet samples.tsv --out comparison/ --regions-file regions.tsv

DEPENDENCIES
------------
    pip install numpy scipy matplotlib seaborn scikit-learn pysam
"""

import argparse
import concurrent.futures
import gzip
import json
import os
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

# scipy >= 1.14 removed SpearmanRConstantInputWarning; define a fallback
try:
    from scipy.stats import SpearmanRConstantInputWarning
except ImportError:
    class SpearmanRConstantInputWarning(Warning):
        pass

from bearing.utils import sanitize_token


def clustered_sample_order(distance_matrix):
    """
    Return row/column order from hierarchical clustering of a symmetric
    distance matrix. Falls back to identity order on any failure.
    """
    n = distance_matrix.shape[0]
    if n <= 2:
        return list(range(n))
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform

        dm = np.asarray(distance_matrix, dtype=np.float64)
        dm = np.nan_to_num(dm, nan=0.0, posinf=0.0, neginf=0.0)
        dm = 0.5 * (dm + dm.T)
        np.fill_diagonal(dm, 0.0)
        condensed = squareform(dm, checks=False)
        Z = linkage(condensed, method="average")
        order = leaves_list(Z).tolist()
        if len(order) == n:
            return order
    except Exception:
        pass
    return list(range(n))

# ---------------------------------------------------------------------------
# Category definitions (must match bigwig_to_qcat.py)
# ---------------------------------------------------------------------------
ALL_CATEGORIES = [
    ("ATAC",                       "#be92e0"),
    ("RNAseq +",                   "#6495ed"),
    ("RNAseq -",                   "#1a3a8f"),
    ("CTCF",                       "#ff2200"),
    ("Cohesin",                    "#8b0000"),
    ("H3K27ac",                    "#00e676"),
    ("Enhancer 4C",                "#c9a0ff"),
    ("DJ1 4C",                     "#7b2fff"),
    ("DJ2 4C",                     "#4b0082"),
    ("Bivalent/Poised TSS",        "#cd5c5c"),
    ("Flanking Bivalent TSS/Enh",  "#e9967a"),
    ("Bivalent Enhancer",          "#bdb76b"),
    ("Repressed PolyComb",         "#808080"),
    ("Weak Repressed PolyComb",    "#c0c0c0"),
    ("Quiescent/Low",              "#ffffff"),
]
MAX_STATES = 15


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


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def load_sample_sheet(tsv_path):
    """
    Parse the sample sheet TSV.

    Supports two flavours of sheet:
      * comparison sheet (existing format) with columns:
            sample, condition, replicate, qcat
      * raw BigWig sheet (new) with columns:
            sample (optional), condition, replicate (optional),
            bw  (comma-separated list) or bw1,bw2,...,
            out  (optional qcat output name)

    The returned list of dicts will always contain the keys used later by
    the workflow (sample, condition, replicate, qcat).  When the input
    sheet supplies BigWigs instead of qcat paths, the qcat field is left
    as a Path to the future output; conversion to qcat files is performed
    later in run().  Any intermediate fields ("bw_paths" or "out") are
    preserved as well.
    """
    samples = []
    with open(tsv_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = [f.strip().lower() for f in fields]
                continue
            row = dict(zip(header, [f.strip() for f in fields]))
            # always require condition
            if "condition" not in row or not row["condition"]:
                sys.exit("ERROR: each row must specify a condition")
            # sample label optional
            sample = row.get("sample") or row.get("name") or ""
            row["sample"] = sample
            # replicate optional
            if "replicate" in row and row["replicate"]:
                try:
                    row["replicate"] = int(row["replicate"])
                except ValueError:
                    row["replicate"] = 0
            else:
                row["replicate"] = 0
            # check for qcat vs bigwig columns
            if "qcat" in row and row["qcat"]:
                p = Path(row["qcat"])
                if not p.exists():
                    sys.exit(f"ERROR: qcat file not found: {p}")
                row["qcat"] = p
            else:
                # expect bigwig paths
                bw_list = []
                if "bw" in row and row["bw"]:
                    bw_list = [p.strip() for p in row["bw"].split(",") if p.strip()]
                else:
                    # collect bw1,bw2,...
                    idx = 1
                    while True:
                        key = f"bw{idx}"
                        if key not in row:
                            break
                        if row[key]:
                            bw_list.append(row[key])
                        idx += 1
                if not bw_list:
                    sys.exit(f"ERROR: no qcat or BigWig paths found for sample {sample}")
                row["bw_paths"] = [Path(p) for p in bw_list]
                # determine output name later, store if present
                if "out" in row and row["out"]:
                    row["out"] = row["out"]
            samples.append(row)

    if len(samples) < 2:
        sys.exit("ERROR: Need at least 2 samples in the sample sheet.")

    return samples


def parse_qcat_bgz(bgz_path, chroms=None, include_raw=False):
    """
    Parse a qcat.bgz file and return per-bin score vectors.

    Each row of the qcat file has the format:
        chr  start  end  id:N,qcat:[[score,state_idx],...]

    The (score, state_idx) pairs are KL scores sorted descending.
    We reconstruct the full num_states score vector per bin,
    then convert back to an approximate probability via softmax so
    JSD is meaningful even if state counts differ across files.

    Returns
    -------
    bins    : list of (chrom, start, end) tuples, in file order
    scores  : np.ndarray of shape (num_bins, num_states), float32
              KL scores (NOT probabilities -- we keep them as scores for
              Spearman and PCA; JSD uses a softmax-derived P)
    num_states : int
    raw_by_bin : dict {(chrom, start, end): np.ndarray|None}
                 Only returned when include_raw=True.
    """
    import pysam

    bins   = []
    rows   = []
    num_states_seen = set()
    raw_by_bin = {} if include_raw else None
    missing_raw_seen = False

    tbx = pysam.TabixFile(str(bgz_path))
    chroms_to_fetch = chroms if chroms else tbx.contigs

    def _parse_qcat_meta(meta):
        if meta.startswith("{"):
            payload = json.loads(meta)
            return payload, payload.get("qcat", []), payload.get("raw"), False

        qcat_start = meta.find("qcat:")
        if qcat_start < 0:
            return None, [], None, True
        raw_start = meta.find(",raw:", qcat_start)
        if raw_start >= 0:
            qcat_payload = meta[qcat_start + 5:raw_start]
            raw_payload = meta[raw_start + 5:]
        else:
            qcat_payload = meta[qcat_start + 5:]
            raw_payload = None
        return None, json.loads(qcat_payload), raw_payload, True

    for chrom in chroms_to_fetch:
        if chrom not in tbx.contigs:
            continue
        for rec in tbx.fetch(chrom):
            parts = rec.split("\t")
            if len(parts) < 4:
                continue
            c, s, e = parts[0], int(parts[1]), int(parts[2])
            meta = parts[3]
            _, pairs, raw_payload, legacy_format = _parse_qcat_meta(meta)
            # pairs: [[score, state_idx], ...]  (1-based state_idx)
            max_state = max(p[1] for p in pairs)
            num_states_seen.add(max_state)
            bins.append((c, s, e))
            rows.append((pairs, max_state))
            if include_raw:
                raw_arr = None
                if raw_payload is not None:
                    try:
                        raw_arr = np.asarray(
                            json.loads(raw_payload) if legacy_format else raw_payload,
                            dtype=np.float64,
                        )
                    except Exception:
                        raw_arr = None
                if raw_arr is None:
                    missing_raw_seen = True
                raw_by_bin[(c, s, e)] = raw_arr

    tbx.close()

    if not bins:
        sys.exit(f"ERROR: No bins parsed from {bgz_path}")

    # Determine consistent num_states across all rows
    num_states = max(num_states_seen)

    score_matrix = np.zeros((len(bins), num_states), dtype=np.float32)
    for i, (pairs, _) in enumerate(rows):
        for score_val, state_idx in pairs:
            si = int(state_idx) - 1   # to 0-based
            if 0 <= si < num_states:
                score_matrix[i, si] = float(score_val)

    if include_raw and missing_raw_seen:
        print(f"  WARNING: raw: field missing in one or more bins for {bgz_path}")

    if include_raw:
        return bins, score_matrix, num_states, raw_by_bin
    return bins, score_matrix, num_states


def load_q_from_cats_json(qcat_path):
    """
    Load per-track Q values from the _cats.json beside a qcat file.
    Returns numpy array of shape (n_tracks,) or None if not found.
    """
    base = str(qcat_path).replace(".qcat.bgz", "").replace(".bgz", "")
    cats_path = base + "_cats.json"
    if not Path(cats_path).exists():
        return None
    with open(cats_path) as f:
        cats = json.load(f)

    # New-style payload written by bigwig_to_qcat.py
    if isinstance(cats, dict) and "q_values_by_assay" in cats:
        q_map = cats.get("q_values_by_assay") or {}
        if not isinstance(q_map, dict) or not q_map:
            return None

        q_vals = []
        # Preserve numeric-state ordering when available via top-level categories.
        categories = cats.get("categories")
        if isinstance(categories, dict) and categories:
            ordered_names = []
            try:
                for k in sorted(categories.keys(), key=lambda x: int(x)):
                    entry = categories[k]
                    if isinstance(entry, dict):
                        ordered_names.append(entry.get("name"))
                    elif isinstance(entry, (list, tuple)) and entry:
                        ordered_names.append(entry[0])
                    else:
                        ordered_names.append(None)
            except Exception:
                ordered_names = []
            if ordered_names and all(name is not None and name in q_map for name in ordered_names):
                q_vals = [q_map[name] for name in ordered_names]
            else:
                q_vals = [q_map[key] for key in sorted(q_map.keys())]
        else:
            q_vals = [q_map[key] for key in sorted(q_map.keys())]

        if any(v is None for v in q_vals):
            return None
        return np.array(q_vals, dtype=np.float64)

    # Legacy numeric-key payloads: {"1": {"q": ...}, ...}
    if isinstance(cats, dict):
        try:
            q_vals = []
            for k in sorted(cats.keys(), key=lambda x: int(x)):
                entry = cats[k]
                if isinstance(entry, dict):
                    q_vals.append(entry.get("q", None))
                elif isinstance(entry, (list, tuple)) and len(entry) >= 3:
                    q_vals.append(entry[2])
                else:
                    q_vals.append(None)
            if any(v is None for v in q_vals):
                return None
            return np.array(q_vals, dtype=np.float64)
        except Exception:
            return None

    return None


def _compute_q_from_raw_bins(raw_bin_dict, pseudocount=1e-6):
    """Estimate per-sample Q from per-bin raw signal vectors."""
    rows = [v for v in raw_bin_dict.values() if v is not None]
    if not rows:
        return None
    raw = np.asarray(rows, dtype=np.float64)
    raw = np.clip(raw, 0.0, None)
    mat = raw + pseudocount
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums <= 0.0, 1.0, row_sums)
    P = mat / row_sums
    return P.mean(axis=0)


def compute_consensus_q(samples, raw_bin_dicts):
    """
    Compute element-wise median Q across all samples.

    Parameters
    ----------
    samples       : list of sample dicts with 'qcat' path
    raw_bin_dicts : list of {bin_key: raw_array} one per sample

    Returns
    -------
    Q_consensus : numpy array of shape (n_tracks,)
    q_source    : str describing how Q was derived ('cats_json' or 'raw_signal')
    per_sample_q: list of numpy arrays, one per sample (for JSD QC reporting)
    """
    per_sample_q = [load_q_from_cats_json(s["qcat"]) for s in samples]
    if all(q is not None for q in per_sample_q):
        n_tracks = len(per_sample_q[0])
        for i, q in enumerate(per_sample_q):
            if len(q) != n_tracks:
                sys.exit(
                    "ERROR: Inconsistent Q length in _cats.json files; "
                    f"expected {n_tracks}, got {len(q)} for {samples[i]['sample']}"
                )
        q_mat = np.vstack(per_sample_q)
        return np.median(q_mat, axis=0), "cats_json", per_sample_q

    print("  INFO: Falling back to raw-signal-derived Q (missing/incomplete _cats.json).")
    per_sample_q = []
    missing_raw_samples = []
    for sample, raw_bin_dict in zip(samples, raw_bin_dicts):
        q = _compute_q_from_raw_bins(raw_bin_dict)
        if q is None:
            missing_raw_samples.append(sample["sample"])
        per_sample_q.append(q)

    if missing_raw_samples:
        sys.exit(
            "ERROR: --consensus-q requires either complete _cats.json Q values "
            "for all samples, or raw:[...] fields in qcat bins. Missing raw for: "
            + ", ".join(missing_raw_samples)
        )

    n_tracks = len(per_sample_q[0])
    for i, q in enumerate(per_sample_q):
        if len(q) != n_tracks:
            sys.exit(
                "ERROR: Inconsistent raw-derived Q length across samples; "
                f"expected {n_tracks}, got {len(q)} for {samples[i]['sample']}"
            )
    q_mat = np.vstack(per_sample_q)
    return np.median(q_mat, axis=0), "raw_signal", per_sample_q


def rescore_bins_with_consensus_q(raw_bin_dict, Q_consensus,
                                  pseudocount=1e-6, min_signal=0.01):
    """
    Rescore all bins using the consensus Q vector.

    Parameters
    ----------
    raw_bin_dict : {bin_key: raw_array} (raw_array may be None for old-format bins)
    Q_consensus  : (n_tracks,) float array

    Returns
    -------
    score_dict : {bin_key: score_array}  (same shape as original parse output)
    n_rescored : int
    n_skipped  : int  (bins where raw was None or low-signal)
    """
    Q = np.asarray(Q_consensus, dtype=np.float64)
    Q_safe = np.clip(Q, 1e-300, None)
    score_dict = {}
    n_rescored = 0
    n_skipped = 0

    for bin_key, raw_arr in raw_bin_dict.items():
        if raw_arr is None:
            score_dict[bin_key] = np.zeros_like(Q_safe, dtype=np.float32)
            n_skipped += 1
            continue

        r = np.clip(np.asarray(raw_arr, dtype=np.float64), 0.0, None)
        if r.shape[0] != Q_safe.shape[0]:
            raise ValueError(
                f"Raw vector length {r.shape[0]} does not match consensus Q length {Q_safe.shape[0]}"
            )

        if float(r.sum()) < float(min_signal):
            score_dict[bin_key] = np.zeros_like(Q_safe, dtype=np.float32)
            n_skipped += 1
            continue

        r_pseudo = r + float(pseudocount)
        denom = float(r_pseudo.sum())
        if denom <= 0.0:
            score_dict[bin_key] = np.zeros_like(Q_safe, dtype=np.float32)
            n_skipped += 1
            continue

        P = r_pseudo / denom
        scores = P * np.log2((P + 1e-300) / Q_safe)
        scores = np.clip(scores, 0.0, None)
        score_dict[bin_key] = scores.astype(np.float32)
        n_rescored += 1

    return score_dict, n_rescored, n_skipped


def _normalize_dist(vec):
    """Normalize a non-negative vector to sum to 1 for QC diagnostics."""
    x = np.clip(np.asarray(vec, dtype=np.float64), 0.0, None)
    s = float(x.sum())
    if s <= 0.0:
        return np.full_like(x, 1.0 / max(1, x.size))
    return x / s


def scores_to_prob(score_matrix):
    """
    Convert KL score matrix to probability distributions via softmax.
    Used only for JSD computation.
    Clips negatives first (already 0 in KL scores), adds pseudocount.
    """
    eps = 1e-9
    P = np.clip(score_matrix, 0.0, None) + eps
    row_sums = P.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums < eps, 1.0, row_sums)
    return P / row_sums


def build_global_nonzero_mask(score_mats, mode="any"):
    """
    Build a bin-level mask across all samples based on aggregate non-zero signal.

    mode:
      - "any": keep bin if at least one sample has non-zero total score.
      - "all": keep bin only if all samples have non-zero total score.
    """
    if mode not in {"any", "all"}:
        raise ValueError(f"Unsupported nonzero mode: {mode}")

    # total KL signal per bin per sample (sum over states)
    totals = np.stack([np.sum(np.abs(m), axis=1) for m in score_mats], axis=1)
    active = totals > 0
    if mode == "any":
        return np.any(active, axis=1)
    return np.all(active, axis=1)


def write_nonzero_bin_summary(out_path, sample_names, score_mats, mask, mode):
    """
    Write summary diagnostics for global nonzero bin filtering.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    total_bins = int(mask.size)
    kept_bins = int(mask.sum())
    with open(out_path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"filter_mode\t{mode}\n")
        fh.write(f"bins_total\t{total_bins}\n")
        fh.write(f"bins_kept\t{kept_bins}\n")
        frac = (float(kept_bins) / float(total_bins)) if total_bins > 0 else 0.0
        fh.write(f"fraction_kept\t{frac:.6g}\n")

        fh.write("\n")
        fh.write("sample\tnonzero_bins\tfraction_nonzero\n")
        for name, mat in zip(sample_names, score_mats):
            sample_nz = int(np.sum(np.sum(np.abs(mat), axis=1) > 0))
            sample_frac = (float(sample_nz) / float(total_bins)) if total_bins > 0 else 0.0
            fh.write(f"{name}\t{sample_nz}\t{sample_frac:.6g}\n")


# ---------------------------------------------------------------------------
# Alignment: make sure all samples use the same bins
# ---------------------------------------------------------------------------

def align_bins(all_bins_list, all_score_list):
    """
    Given N lists of (chrom, start, end) tuples and N score matrices,
    find the intersection of bins present in ALL samples and return
    aligned score matrices.

    Returns
    -------
    common_bins  : list of (chrom, start, end) tuples
    aligned_mats : list of np.ndarray, each (n_common, num_states)
    """
    print("  Aligning bins across samples...")

    # Build index sets
    sets = [set(b) for b in all_bins_list]
    common = sets[0]
    for s in sets[1:]:
        common &= s
    common_bins = sorted(common, key=lambda x: (
        x[0],
        int(x[0].replace("chr", "").replace("X", "23").replace("Y", "24"))
        if x[0].replace("chr", "").replace("X", "23").replace("Y", "24").isdigit()
        else 0,
        x[1],
    ))
    print(f"  Common bins: {len(common_bins):,}  "
          f"(dropped {len(all_bins_list[0]) - len(common_bins):,} non-overlapping)")

    aligned = []
    for bins, mat in zip(all_bins_list, all_score_list):
        idx = {b: i for i, b in enumerate(bins)}
        row_ids = [idx[b] for b in common_bins]
        aligned.append(mat[row_ids, :])

    return common_bins, aligned


# ---------------------------------------------------------------------------
# Jensen-Shannon divergence
# ---------------------------------------------------------------------------

def js_divergence(P, Q):
    """
    Compute per-bin Jensen-Shannon divergence between two (n, k) probability
    matrices.  Returns scalar mean JSD over all bins.
    Uses base-2 logarithm so JSD in [0, 1].
    """
    eps = 1e-300
    M = 0.5 * (P + Q)
    kl_pm = np.sum(P * np.log2((P + eps) / (M + eps)), axis=1)
    kl_qm = np.sum(Q * np.log2((Q + eps) / (M + eps)), axis=1)
    jsd_per_bin = np.clip(0.5 * kl_pm + 0.5 * kl_qm, 0.0, 1.0)
    return float(np.mean(jsd_per_bin))


def build_jsd_matrix(prob_mats, sample_names):
    """
    Build symmetric (n_samples, n_samples) JSD matrix.
    """
    n = len(prob_mats)
    mat = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            jsd = js_divergence(prob_mats[i], prob_mats[j])
            mat[i, j] = jsd
            mat[j, i] = jsd
    return mat


# ---------------------------------------------------------------------------
# Spearman rho per state
# ---------------------------------------------------------------------------

def build_spearman_matrix(score_mats, sample_names, num_states,
                          nonzero_mode="all", min_bins=2,
                          return_counts=False):
    """
    For each state, compute pairwise Spearman rho across samples.
    Returns a (n_samples, n_samples, num_states) array.
    Also returns a (n_samples, n_samples) mean-rho matrix (mean over states).

    Parameters
    ----------
    nonzero_mode : str
        "all"  -> use all bins
        "any"  -> keep bins where at least one sample in pair is non-zero
        "both" -> keep bins where both samples in pair are non-zero
    min_bins : int
        Minimum bins required for a pair/state correlation.
    return_counts : bool
        If True, also return (bins_used, bins_total) arrays with
        shape (n_samples, n_samples, num_states).
    """
    from scipy.stats import spearmanr
    import warnings

    n = len(score_mats)
    rho_per_state = np.full((n, n, num_states), np.nan, dtype=np.float64)
    rho_mean      = np.zeros((n, n), dtype=np.float64)
    bins_used = np.zeros((n, n, num_states), dtype=np.int32)
    bins_total = np.zeros((n, n, num_states), dtype=np.int32)

    if nonzero_mode not in {"all", "any", "both"}:
        raise ValueError(f"Unsupported nonzero_mode: {nonzero_mode}")
    min_bins = max(2, int(min_bins))

    print("  Computing Spearman correlations...")
    for si in range(num_states):
        cols = np.stack([m[:, si] for m in score_mats], axis=1)  # (bins, n)
        n_total = int(cols.shape[0])

        for i in range(n):
            rho_per_state[i, i, si] = 1.0
            bins_used[i, i, si] = n_total
            bins_total[i, i, si] = n_total

        for i in range(n):
            for j in range(i + 1, n):
                x = cols[:, i]
                y = cols[:, j]
                if nonzero_mode == "all":
                    mask = np.ones(n_total, dtype=bool)
                elif nonzero_mode == "any":
                    mask = (x != 0) | (y != 0)
                else:
                    mask = (x != 0) & (y != 0)

                used = int(mask.sum())
                bins_used[i, j, si] = used
                bins_used[j, i, si] = used
                bins_total[i, j, si] = n_total
                bins_total[j, i, si] = n_total

                if used < min_bins:
                    continue

                # Suppress warning for constant vectors after filtering.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", SpearmanRConstantInputWarning)
                    rho_ij, _ = spearmanr(x[mask], y[mask])
                if np.isnan(rho_ij):
                    continue
                rho_per_state[i, j, si] = float(rho_ij)
                rho_per_state[j, i, si] = float(rho_ij)

    with np.errstate(invalid="ignore"):
        rho_mean = np.nanmean(rho_per_state, axis=2)
    rho_mean = np.nan_to_num(rho_mean, nan=0.0)
    np.fill_diagonal(rho_mean, 1.0)

    if return_counts:
        return rho_per_state, rho_mean, bins_used, bins_total
    return rho_per_state, rho_mean


def build_total_saliency_spearman(score_mats, sample_names,
                                  nonzero_mode="all", min_bins=2):
    """
    Compute pairwise Spearman rho on per-bin total saliency.

    Total saliency per sample is sum over states per bin.

    Parameters
    ----------
    nonzero_mode : str
        "all"  -> use all bins
        "any"  -> keep bins where at least one sample in pair is non-zero
        "both" -> keep bins where both samples in pair are non-zero
    min_bins : int
        Minimum bins required for pairwise correlation.
    """
    from scipy.stats import spearmanr
    import warnings

    n = len(score_mats)
    rho = np.zeros((n, n), dtype=np.float64)
    np.fill_diagonal(rho, 1.0)

    if nonzero_mode not in {"all", "any", "both"}:
        raise ValueError(f"Unsupported nonzero_mode: {nonzero_mode}")
    min_bins = max(2, int(min_bins))

    totals = np.stack([np.sum(m, axis=1) for m in score_mats], axis=1)
    n_total = int(totals.shape[0])

    for i in range(n):
        for j in range(i + 1, n):
            x = totals[:, i]
            y = totals[:, j]

            if nonzero_mode == "all":
                mask = np.ones(n_total, dtype=bool)
            elif nonzero_mode == "any":
                mask = (x != 0) | (y != 0)
            else:
                mask = (x != 0) & (y != 0)

            if int(mask.sum()) < min_bins:
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", SpearmanRConstantInputWarning)
                rho_ij, _ = spearmanr(x[mask], y[mask])

            if np.isnan(rho_ij):
                continue

            rho[i, j] = float(rho_ij)
            rho[j, i] = float(rho_ij)

    return rho


def write_spearman_nonzero_diagnostics(out_path, sample_names, categories,
                                       rho_all, rho_nonzero,
                                       bins_used, bins_total):
    """
    Write per-pair, per-state Spearman comparison for all-bins vs nonzero-only.
    """
    n = len(sample_names)
    n_states = rho_all.shape[2]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(
            "sample_a\tsample_b\tstate_index\tstate_name\tbins_total"
            "\tbins_used_nonzero\tfraction_used_nonzero\trho_all\trho_nonzero\n"
        )
        for i in range(n):
            for j in range(i + 1, n):
                for si in range(n_states):
                    state_name = categories[si][0] if si < len(categories) else f"State {si+1}"
                    total = int(bins_total[i, j, si])
                    used = int(bins_used[i, j, si])
                    frac = (float(used) / float(total)) if total > 0 else 0.0
                    fh.write(
                        f"{sample_names[i]}\t{sample_names[j]}\t{si+1}\t{state_name}"
                        f"\t{total}\t{used}\t{frac:.6g}\t{float(rho_all[i, j, si]):.6g}"
                        f"\t{float(rho_nonzero[i, j, si]):.6g}\n"
                    )


def write_pairwise_matrix_tsv(out_path, sample_names, matrix):
    """Write a square sample x sample matrix TSV with sample-name labels."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("sample\t" + "\t".join(sample_names) + "\n")
        for i, name in enumerate(sample_names):
            vals = "\t".join(f"{float(matrix[i, j]):.6g}" for j in range(len(sample_names)))
            fh.write(f"{name}\t{vals}\n")


def summarize_spearman_nonzero_effect(rho_all_mean, rho_nonzero_mean,
                                      bins_used, bins_total):
    """
    Summarize how much nonzero-only filtering changed Spearman inputs/results.

    Returns dict with fraction-of-bins-used stats and delta stats for
    mean-rho matrices (upper triangle, excluding diagonal).
    """
    n = int(rho_all_mean.shape[0])
    frac_vals = []
    for i in range(n):
        for j in range(i + 1, n):
            used = bins_used[i, j, :].astype(np.float64)
            total = bins_total[i, j, :].astype(np.float64)
            good = total > 0
            if np.any(good):
                frac_vals.extend((used[good] / total[good]).tolist())

    if frac_vals:
        frac_arr = np.asarray(frac_vals, dtype=np.float64)
        frac_min = float(np.min(frac_arr))
        frac_med = float(np.median(frac_arr))
        frac_max = float(np.max(frac_arr))
    else:
        frac_min = np.nan
        frac_med = np.nan
        frac_max = np.nan

    iu = np.triu_indices(n, k=1)
    delta = rho_nonzero_mean[iu] - rho_all_mean[iu]
    abs_delta = np.abs(delta)
    if abs_delta.size:
        mean_abs_delta = float(np.mean(abs_delta))
        max_abs_delta = float(np.max(abs_delta))
    else:
        mean_abs_delta = 0.0
        max_abs_delta = 0.0

    return {
        "frac_min": frac_min,
        "frac_med": frac_med,
        "frac_max": frac_max,
        "mean_abs_delta": mean_abs_delta,
        "max_abs_delta": max_abs_delta,
    }


# ---------------------------------------------------------------------------
# Differential epilogos (signed KL)
# ---------------------------------------------------------------------------

def diff_qcat(bins, mean_score_A, mean_score_B, out_path):
    """
    Compute signed differential KL score per bin per state:
        diff_i = KL_A_i - KL_B_i

    Positive = state i is more active in condition A.
    Negative = state i is more active in condition B.

    Writes a qcat.bgz with the same format as bigwig_to_qcat.py output.
    The scores are signed floats; pyGenomeTracks epilogos will render
    positive bars above the axis and negative below.
    """
    try:
        import pysam
    except ImportError:
        print("  WARNING: pysam not available, skipping diff qcat output.")
        return

    diff = mean_score_A - mean_score_B   # (bins, num_states)

    tmp_path = str(out_path) + ".tmp.tsv"
    with open(tmp_path, "w") as fh:
        for idx, (chrom, start, end) in enumerate(bins):
            row_diff = diff[idx]
            pairs = sorted(
                ([float(f"{v:.6g}"), si + 1]
                 for si, v in enumerate(row_diff)),
                key=lambda x: -abs(x[0]),
            )
            qcat_col = "id:" + str(idx + 1) + ",qcat:" + json.dumps(pairs, separators=(",", ":"))
            fh.write(chrom + "\t" + str(start) + "\t" + str(end) +
                     "\t" + qcat_col + "\n")

    # Sort
    import subprocess
    sorted_path = tmp_path + ".sorted"
    with open(sorted_path, "w") as sf:
        subprocess.run(["sort", "-k1,1", "-k2,2n", tmp_path],
                       stdout=sf, check=True)
    os.remove(tmp_path)

    bgz_path = str(out_path)
    if not bgz_path.endswith(".bgz"):
        bgz_path += ".bgz"

    pysam.tabix_compress(sorted_path, bgz_path, force=True)
    pysam.tabix_index(bgz_path, preset="bed", force=True)
    os.remove(sorted_path)

    print(f"  Diff qcat:  {bgz_path}")


# ---------------------------------------------------------------------------
# Absolute-signal floor for differential bins
# ---------------------------------------------------------------------------

def _bin_total_signal(raw_arr):
    """Total raw signal (coverage) at a bin = sum of the per-state raw vector."""
    if raw_arr is None:
        return 0.0
    a = np.asarray(raw_arr, dtype=np.float64)
    if a.size == 0:
        return 0.0
    return float(np.clip(a, 0.0, None).sum())


def condition_signal_vectors(common_bins, all_raw_by_bin, cond_groups):
    """
    Per-condition mean per-bin total raw signal, aligned to common_bins.

    Signal at a bin is the sum of its raw per-state vector (total coverage).
    The per-condition value is the mean across that condition's replicates,
    matching how cond_means averages scores. Bins lacking a raw: field count
    as 0 signal (cannot be assessed -> treated as below any positive floor).

    Returns {condition: np.ndarray of shape (n_common_bins,)}.
    """
    n = len(common_bins)
    sample_sig = []
    for raw_by_bin in all_raw_by_bin:
        sums = {b: _bin_total_signal(a) for b, a in raw_by_bin.items()}
        sample_sig.append(
            np.fromiter((sums.get(b, 0.0) for b in common_bins),
                        dtype=np.float64, count=n)
        )
    sig_by_cond = {}
    for cond, idxs in cond_groups.items():
        sig_by_cond[cond] = np.mean([sample_sig[i] for i in idxs], axis=0)
    return sig_by_cond


def condition_signal_vectors_bigwig(common_bins, samples, cond_groups,
                                    track_aggregate="sum", rep_aggregate="mean",
                                    sheet_dir=None):
    """
    Per-condition per-bin total signal read from each sample's BigWigs,
    aligned to common_bins. Used when the qcat files carry no raw: field
    (the canonical pyGenomeTracks-compatible format), so the absolute-signal
    floor still has a real magnitude to threshold on.

    Within a sample, the per-track mean coverage at a bin is combined by
    track_aggregate ('sum' = total activity across assays; 'max' = strongest).
    Within a condition, samples are combined by rep_aggregate ('mean' matches
    how cond_means averages scores; 'max' keeps a bin if any replicate has
    signal). Bigwig paths come from each sample row's 'bw' field; relative
    paths resolve against sheet_dir.

    Returns {condition: np.ndarray (n_common_bins,)}.
    """
    import pyBigWig

    # Group common_bins by chromosome with their row index, for batched reads.
    by_chrom = {}
    for i, (c, s, e) in enumerate(common_bins):
        by_chrom.setdefault(c, []).append((i, int(s), int(e)))
    n = len(common_bins)

    def sample_bw_paths(srow):
        raw = srow.get("bw")
        if raw:
            paths = [p.strip() for p in str(raw).split(",") if p.strip()]
        else:
            paths = [str(p) for p in srow.get("bw_paths", [])]
        out = []
        for p in paths:
            pp = p
            if sheet_dir and not os.path.isabs(p):
                pp = os.path.join(sheet_dir, p)
            out.append(pp)
        return out

    def read_sample(srow):
        paths = sample_bw_paths(srow)
        if not paths:
            return None
        handles = [pyBigWig.open(p) for p in paths]
        try:
            tot = np.zeros(n, dtype=np.float64)
            for h in handles:
                chd = h.chroms()
                col = np.zeros(n, dtype=np.float64)
                for c, items in by_chrom.items():
                    if c not in chd:
                        continue
                    clen = int(chd[c])
                    for i, s, e in items:
                        if s >= clen:
                            continue
                        ee = min(e, clen)
                        try:
                            v = h.stats(c, s, ee, type="mean")
                        except RuntimeError:
                            v = None
                        val = v[0] if (v and v[0] is not None) else 0.0
                        if val > 0:
                            col[i] = val
                tot = np.maximum(tot, col) if track_aggregate == "max" else tot + col
            return tot
        finally:
            for h in handles:
                try:
                    h.close()
                except Exception:
                    pass

    sample_sig = [None] * len(samples)
    sig_by_cond = {}
    for cond, idxs in cond_groups.items():
        arrs = []
        for i in idxs:
            if sample_sig[i] is None:
                sample_sig[i] = read_sample(samples[i])
            if sample_sig[i] is not None:
                arrs.append(sample_sig[i])
        if not arrs:
            sig_by_cond[cond] = np.zeros(n, dtype=np.float64)
            continue
        stacked = np.vstack(arrs)
        if rep_aggregate == "max":
            sig_by_cond[cond] = stacked.max(axis=0)
        elif rep_aggregate == "min":
            sig_by_cond[cond] = stacked.min(axis=0)
        else:
            sig_by_cond[cond] = stacked.mean(axis=0)
    return sig_by_cond


def apply_diff_signal_floor(common_bins, mean_A, mean_B, sig_A, sig_B,
                            floor, combiner):
    """
    Keep only bins whose combined per-condition signal clears the floor.

    combiner 'max': keep where max(sig_A, sig_B) >= floor. This retains
                    genuine on/off events (signal present in only one
                    condition) while removing double-low noise-vs-noise bins.
    combiner 'sum': keep where sig_A + sig_B >= floor.

    Returns (kept_bins, mean_A[keep], mean_B[keep], keep_mask).
    """
    if floor is None or floor <= 0.0:
        keep = np.ones(len(common_bins), dtype=bool)
    elif combiner == "sum":
        keep = (sig_A + sig_B) >= float(floor)
    else:
        keep = np.maximum(sig_A, sig_B) >= float(floor)
    if keep.all():
        return common_bins, mean_A, mean_B, keep
    kept_bins = [b for b, k in zip(common_bins, keep) if k]
    return kept_bins, mean_A[keep], mean_B[keep], keep


def write_tested_bins_bed(bed_path, bins):
    """Write the surviving (tested) bin set as a BED so the permutation null
    can be restricted to the identical fixed bin set (blacklist-style)."""
    with open(bed_path, "w") as fh:
        for chrom, start, end in bins:
            fh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\n")


def _dump_signal_hist(out_dir, common_bins, sig_by_cond, cond_groups, combiner):
    """
    Dump per-bin combined signal and percentile/threshold summaries to help
    choose --diff-min-signal from the actual genome-wide distribution.

    Writes two files in out_dir:
      signal_hist.tsv     : per (cond_A, cond_B) condition pair, one row per
                            percentile of max(sig_A,sig_B) (or sum), and the
                            implied bins-kept count for that threshold.
      signal_per_bin.tsv  : per (chrom, start, end), the per-condition signal
                            and the combined value for every condition pair.
                            (Genome-wide; expect millions of rows.)
    Both are TSVs you can load with pandas / awk.
    """
    out_dir = Path(out_dir)
    conds = sorted(cond_groups.keys())
    n = len(common_bins)
    PCTS = [1, 5, 10, 25, 50, 75, 90, 95, 99, 99.5, 99.9]

    hist_path = out_dir / "signal_hist.tsv"
    with open(hist_path, "w") as fh:
        fh.write("\t".join(["cond_A", "cond_B", "combiner", "n_bins",
                            "percentile", "threshold",
                            "bins_at_or_above", "frac_at_or_above"]) + "\n")
        for ci, cA in enumerate(conds):
            for cB in conds[ci + 1:]:
                sA = sig_by_cond[cA]
                sB = sig_by_cond[cB]
                if combiner == "sum":
                    comb = sA + sB
                else:
                    comb = np.maximum(sA, sB)
                # percentile rows
                for p in PCTS:
                    thr = float(np.percentile(comb, p))
                    kept = int(np.sum(comb >= thr))
                    fh.write("\t".join([cA, cB, combiner, str(n),
                                        "{:.3g}".format(p),
                                        "{:.6g}".format(thr),
                                        str(kept),
                                        "{:.4g}".format(kept / max(n, 1))]) + "\n")
                # exemplar fixed thresholds people actually try
                for thr in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0,
                            2.0, 5.0, 10.0]:
                    if thr > float(np.max(comb)):
                        continue
                    kept = int(np.sum(comb >= thr))
                    fh.write("\t".join([cA, cB, combiner, str(n),
                                        "fixed",
                                        "{:.6g}".format(thr),
                                        str(kept),
                                        "{:.4g}".format(kept / max(n, 1))]) + "\n")
    print("  Signal histogram: " + str(hist_path))

    # Per-bin dump (genome-wide). Big; useful for downstream eyeballing.
    perbin_path = out_dir / "signal_per_bin.tsv"
    with open(perbin_path, "w") as fh:
        header = ["chrom", "start", "end"]
        for c in conds:
            header.append("sig_" + c)
        pair_cols = []
        for ci, cA in enumerate(conds):
            for cB in conds[ci + 1:]:
                pair_cols.append((cA, cB))
                header.append("comb_" + cA + "_vs_" + cB)
        fh.write("\t".join(header) + "\n")
        cond_arrays = [sig_by_cond[c] for c in conds]
        comb_arrays = []
        for cA, cB in pair_cols:
            sA = sig_by_cond[cA]; sB = sig_by_cond[cB]
            comb_arrays.append(sA + sB if combiner == "sum"
                               else np.maximum(sA, sB))
        for i, (ch, s, e) in enumerate(common_bins):
            row = [ch, str(s), str(e)]
            for arr in cond_arrays:
                row.append("{:.6g}".format(float(arr[i])))
            for arr in comb_arrays:
                row.append("{:.6g}".format(float(arr[i])))
            fh.write("\t".join(row) + "\n")
    print("  Signal per-bin: " + str(perbin_path))


# ---------------------------------------------------------------------------
# INI writer
# ---------------------------------------------------------------------------

def write_compare_ini(ini_path, out_dir, samples, conditions,
                      qcat_paths, diff_pairs, cats_json_path, num_states,
                      genes_path=None, beds=None):
    """
    Write a pyGenomeTracks INI file for the comparison output.

    Layout per condition pair (interleaved):
      [spacer]
      # --- Diff: A vs B ---
      [epilogos diff_A_vs_B]          <- positive scores (A > B)
      [epilogos diff_A_vs_B inverted] <- negative scores (B > A)
      [spacer]
      # Condition A samples
      [epilogos sample_A_rep1]
      [epilogos sample_A_rep1 inv]
      [spacer]
      # Condition B samples
      [epilogos sample_B_rep1]
      [epilogos sample_B_rep1 inv]
      [spacer]
    [x-axis]

    diff_pairs     : list of (cond_A, cond_B, bgz_path) tuples
    cats_json_path : path to the shared categories JSON file
    """
    from collections import defaultdict as _dd

    lines = []
    lines.append("# pyGenomeTracks configuration -- generated by compare_qcat.py")
    lines.append("# Differential and per-sample epilogos tracks.")
    lines.append("# Usage:")
    lines.append("#   pyGenomeTracks --tracks " + str(ini_path)
                 + " --region chr1:10000000-11000000 -o out.png")
    lines.append("")

    # Group samples by condition for easy lookup
    cond_to_samples = _dd(list)
    for s, c, q in zip(
        [s["sample"] for s in samples],
        [s["condition"] for s in samples],
        qcat_paths,
    ):
        cond_to_samples[c].append((s, q))

    def _safe(name):
        return sanitize_token(name)

    def _epilogos_block(track_id, bgz_path, title, height=4, inverted=False):
        block = []
        block.append("[" + track_id + "]")
        block.append("file = " + str(bgz_path))
        block.append("categories_file = " + str(cats_json_path))
        block.append("height = " + str(height))
        block.append("title = " + title)
        if inverted:
            block.append("orientation = inverted")
        block.append("")
        return block

    def _spacer(h=0.3):
        return ["[spacer]", "height = " + str(h), ""]

    # --- Interleaved: per condition pair ---
    for cond_A, cond_B, diff_bgz in diff_pairs:
        safe_A = _safe(cond_A)
        safe_B = _safe(cond_B)
        pair_label = cond_A + " vs " + cond_B

        lines.append("# " + "=" * 56)
        lines.append("# " + pair_label)
        lines.append("# " + "=" * 56)
        lines.append("")

        # Differential track
        lines.append("# Positive bars = enriched in " + cond_A)
        lines += _epilogos_block(
            "epilogos diff_" + safe_A + "_vs_" + safe_B,
            diff_bgz,
            "Diff " + pair_label + " (+" + cond_A + ")",
            height=4,
        )
        lines.append("# Positive bars = enriched in " + cond_B)
        lines += _epilogos_block(
            "epilogos diff_" + safe_A + "_vs_" + safe_B + "_inv",
            diff_bgz,
            "Diff " + pair_label + " (+" + cond_B + ")",
            height=4,
            inverted=True,
        )
        lines += _spacer(0.5)

        # Per-sample tracks, condition A then condition B
        for cond in (cond_A, cond_B):
            lines.append("# --- " + cond + " samples ---")
            lines.append("")
            for sname, qpath in cond_to_samples[cond]:
                safe_s = _safe(sname)
                lines += _epilogos_block(
                    "epilogos " + safe_s,
                    qpath,
                    sname + " (" + cond + ")",
                    height=3,
                )
                lines += _epilogos_block(
                    "epilogos " + safe_s + "_inv",
                    qpath,
                    sname + " (" + cond + ") inv",
                    height=3,
                    inverted=True,
                )
                lines += _spacer(0.2)
            lines += _spacer(0.5)

    if genes_path:
        gp = Path(genes_path)
        gene_label = gp.stem.replace("_", " ")
        lines.append("# Gene annotation track")
        lines.append("[bed " + gene_label + "]")
        lines.append("file = " + str(genes_path))
        lines.append("title = " + gene_label)
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

    # Additional BED tracks (generic BED; user can pass multiple via --bed)
    if beds:
        for bedp in beds:
            bp = Path(bedp)
            bed_label = bp.stem.replace("_", " ")
            lines.append("# BED track")
            lines.append("[bed " + bed_label + "]")
            lines.append("file = " + str(bedp))
            lines.append("title = " + bed_label)
            lines.append("height = 3")
            lines.append("fontsize = 10")
            lines.append("file_type = bed")
            lines.append("color = #000000")
            lines.append("border_color = black")
            lines.append("display = stacked")
            lines.append("style = UCSC")
            lines.append("show_labels = false")
            lines.append("")

    lines.append("[x-axis]")
    lines.append("")

    ini_path = Path(ini_path)
    with open(ini_path, "w") as fh:
        fh.write("\n".join(lines))
    print("  Comparison INI: " + str(ini_path))


def clip_qcat_for_plot(src_bgz, dst_bgz, qcat_max=5.0, signed=False):
    """
    Create a plotting-only qcat.bgz with scores clipped to a fixed range.

    For sample tracks (signed=False), scores are clipped to [0, qcat_max].
    For differential tracks (signed=True), scores are clipped to
    [-qcat_max, +qcat_max].
    """
    try:
        import pysam
    except ImportError:
        print("  WARNING: pysam not available, skipping clipped qcat generation.")
        return src_bgz

    src_bgz = Path(src_bgz)
    dst_bgz = Path(dst_bgz)
    dst_bgz.parent.mkdir(parents=True, exist_ok=True)

    tmp_tsv = str(dst_bgz) + ".tmp.tsv"
    with gzip.open(src_bgz, "rt") as fin, open(tmp_tsv, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            meta = parts[3]
            if meta.startswith("{"):
                payload = json.loads(meta)
                pairs = payload.get("qcat", [])
            else:
                qcat_start = meta.find("qcat:")
                if qcat_start < 0:
                    fout.write(line + "\n")
                    continue
                raw_start = meta.find(",raw:", qcat_start)
                if raw_start >= 0:
                    qcat_payload = meta[qcat_start + 5:raw_start]
                    raw_suffix = meta[raw_start:]
                else:
                    qcat_payload = meta[qcat_start + 5:]
                    raw_suffix = ""
                pairs = json.loads(qcat_payload)

            clipped = []
            for score, state_idx in pairs:
                v = float(score)
                if signed:
                    v = max(-qcat_max, min(qcat_max, v))
                else:
                    v = max(0.0, min(qcat_max, v))
                clipped.append([float(f"{v:.6g}"), int(state_idx)])

            if meta.startswith("{"):
                payload["qcat"] = clipped
                parts[3] = json.dumps(payload, separators=(",", ":"))
            else:
                prefix = meta[:qcat_start + 5]
                parts[3] = prefix + json.dumps(clipped, separators=(",", ":")) + raw_suffix
            fout.write("\t".join(parts) + "\n")

    pysam.tabix_compress(tmp_tsv, str(dst_bgz), force=True)
    pysam.tabix_index(str(dst_bgz), preset="bed", force=True)
    os.remove(tmp_tsv)
    return dst_bgz


# ---------------------------------------------------------------------------
# Per-region matplotlib comparison figure
# ---------------------------------------------------------------------------

def plot_qcat_region(region_str, samples, qcat_paths, diff_pairs,
                     out_path, genes_path=None, highlights_path=None,
                     qcat_max=5.0, categories=None, region_label=None):
    """
    Plot all sample qcat tracks and diff tracks for a single genomic region
    using matplotlib.

    All sample tracks share a fixed y-limit of qcat_max * 1.15 and all diff
    tracks share a symmetric y-limit of ±qcat_max * 1.15, enabling direct
    visual comparison across tracks and regions.

    Parameters
    ----------
    region_str      : "chr:start-end"
    samples         : list of sample dicts from load_sample_sheet
    qcat_paths      : list of Paths (original, unclipped), parallel to samples
    diff_pairs      : list of (cond_A, cond_B, bgz_path)
    out_path        : output PDF path
    genes_path      : optional BED6 gene annotation file
    highlights_path : optional BED4 highlight regions (col 4 = hex color)
    qcat_max        : shared y-axis ceiling for sample tracks; ± for diffs
    categories      : list of (name, color) tuples; defaults to ALL_CATEGORIES
    region_label    : display label for figure title (defaults to region_str)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import matplotlib.colors as mcolors
    from bearing_hic_plot import (
        parse_region, load_qcat_scores, load_genes, load_genes_gtf, load_highlights,
        draw_epilogos_horizontal, draw_diff_horizontal,
        draw_gene_track, draw_genomic_axis,
    )

    _cats = categories if categories is not None else ALL_CATEGORIES
    chrom, region_start, region_end = parse_region(region_str)
    print(f"    {chrom}:{region_start:,}-{region_end:,}")

    # --- Load qcat scores ---
    sample_tracks = []
    ns_max = 1
    for s, qp in zip(samples, qcat_paths):
        pos, mat, ns = load_qcat_scores(qp, chrom, region_start, region_end)
        ns_max = max(ns_max, ns)
        sample_tracks.append({"label": s["sample"], "condition": s["condition"],
                               "pos": pos, "mat": mat})
    for t in sample_tracks:
        if t["mat"].shape[1] < ns_max:
            t["mat"] = np.pad(t["mat"], ((0, 0), (0, ns_max - t["mat"].shape[1])))

    # --- Load diff tracks ---
    diff_tracks = []
    for cond_A, cond_B, bgz_path in diff_pairs:
        pos, mat, ns = load_qcat_scores(bgz_path, chrom, region_start, region_end)
        if mat.shape[1] < ns_max:
            mat = np.pad(mat, ((0, 0), (0, ns_max - mat.shape[1])))
        diff_tracks.append({"cond_A": cond_A, "cond_B": cond_B,
                             "pos": pos, "mat": mat})

    # --- Load optional annotations ---
    genes = None
    if genes_path:
        gene_path_lower = str(genes_path).lower()
        if gene_path_lower.endswith(('.gtf', '.gff', '.gff3')):
            genes = load_genes_gtf(genes_path, chrom, region_start, region_end)
        else:
            genes = load_genes(genes_path, chrom, region_start, region_end)
    highlights = (load_highlights(highlights_path, chrom, region_start, region_end)
                  if highlights_path else None)
    has_genes = bool(genes)

    # --- Figure layout ---
    n_samples = len(sample_tracks)
    n_diffs   = len(diff_tracks)
    sample_h  = 1.0
    diff_h    = 0.8
    gene_h    = 0.30
    axis_h    = 0.28
    heights = (
        [sample_h] * n_samples
        + [diff_h]  * n_diffs
        + ([gene_h] if has_genes else [])
        + [axis_h]
    )
    n_rows = len(heights)
    fig_h  = sum(heights) + 0.8
    fig = plt.figure(figsize=(10.0, fig_h), dpi=150)
    gs  = gridspec.GridSpec(
        n_rows, 1, figure=fig,
        height_ratios=heights,
        hspace=0.05,
        left=0.14, right=0.97, top=0.96, bottom=0.06,
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_rows)]

    # Title
    title = (f"{region_label}  ({region_str})"
             if (region_label and region_label != region_str) else region_str)
    fig.suptitle(title, fontsize=9, y=0.99)

    # Condition color map for background tinting
    cmap_cond = condition_color_map([t["condition"] for t in sample_tracks])
    y_limit = qcat_max * 1.15   # sample track y ceiling
    d_limit = qcat_max           # diff track ± limit

    # --- Sample tracks ---
    for i, t in enumerate(sample_tracks):
        draw_epilogos_horizontal(
            axes[i], t["pos"], t["mat"], ns_max,
            region_start, region_end,
            _cats[:ns_max],
            highlights=highlights,
            label=f"{t['label']}\n({t['condition']})",
            y_max=y_limit,
        )
        axes[i].set_facecolor(
            mcolors.to_rgba(cmap_cond[t["condition"]], alpha=0.06)
        )

    # --- Diff tracks ---
    for i, t in enumerate(diff_tracks):
        draw_diff_horizontal(
            axes[n_samples + i], t["pos"], t["mat"], ns_max,
            region_start, region_end,
            _cats[:ns_max],
            highlights=highlights,
            label=f"Diff\n{t['cond_A']}\u2212{t['cond_B']}",
            diff_max=d_limit,
        )

    # --- Optional gene track ---
    row = n_samples + n_diffs
    if has_genes:
        draw_gene_track(axes[row], genes, region_start, region_end,
                        highlights=highlights, label="Genes")
        row += 1

    # --- Genomic axis ---
    draw_genomic_axis(axes[row], region_start, region_end, chrom)

    # --- Save ---
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"    Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

PALETTE = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#999999",
    "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
]


def condition_color_map(conditions):
    unique = sorted(set(conditions))
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(unique)}


def plot_jsd_heatmap(jsd_mat, sample_names, conditions, out_path, title=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Cluster samples by JSD distance and reorder matrix accordingly.
    order = clustered_sample_order(jsd_mat)
    jsd_mat = jsd_mat[np.ix_(order, order)]
    sample_names = [sample_names[i] for i in order]
    conditions = [conditions[i] for i in order]

    cmap = condition_color_map(conditions)
    row_colors = [cmap[c] for c in conditions]

    fig, ax = plt.subplots(figsize=(max(6, len(sample_names) * 0.8 + 2),
                                    max(5, len(sample_names) * 0.8 + 1.5)))
    im = ax.imshow(jsd_mat, cmap="viridis", vmin=0, vmax=jsd_mat.max())
    plt.colorbar(im, ax=ax, label="Mean Jensen-Shannon divergence")

    ax.set_xticks(range(len(sample_names)))
    ax.set_yticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(sample_names, fontsize=9)

    # Annotate cells
    for i in range(len(sample_names)):
        for j in range(len(sample_names)):
            val = jsd_mat[i, j]
            color = "white" if val > jsd_mat.max() * 0.6 else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=7, color=color)

    # Color strips on axes for condition
    for tick, col in zip(ax.get_xticklabels(), row_colors):
        tick.set_color(col)
    for tick, col in zip(ax.get_yticklabels(), row_colors):
        tick.set_color(col)

    # Legend for conditions
    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=k) for k, c in cmap.items()]
    ax.legend(handles=handles, title="Condition",
              bbox_to_anchor=(1.25, 1), loc="upper left", fontsize=8)

    ax.set_title(title if title else "Pairwise Jensen-Shannon Divergence\n(lower = more similar)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  JSD heatmap:   {out_path}")


def plot_spearman(rho_mean, rho_per_state, sample_names, conditions,
                  num_states, out_path, categories=None, title=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    _cats = categories if categories is not None else ALL_CATEGORIES

    # Cluster samples using a distance derived from mean Spearman rho.
    # Similarity rho in [-1, 1] -> distance in [0, 2].
    dist = np.clip(1.0 - rho_mean, 0.0, 2.0)

    # Attempt hierarchical clustering to drive a dendrogram. Fall back to
    # the existing clustered_sample_order if SciPy is unavailable.
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
        from scipy.spatial.distance import squareform

        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        order = leaves_list(Z).tolist()
    except Exception:
        order = clustered_sample_order(dist)

    rho_mean = rho_mean[np.ix_(order, order)]
    rho_per_state = rho_per_state[np.ix_(order, order, np.arange(num_states))]
    sample_names = [sample_names[i] for i in order]
    conditions = [conditions[i] for i in order]

    cmap_cond = condition_color_map(conditions)
    n = len(sample_names)
    state_names = [_cats[i][0] for i in range(num_states)]

    # Build a figure with a nested left panel so the dendrogram spans only the
    # mean-heatmap square, not the colorbar/legend area.
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(max(10, n * 0.9 + 6), max(5, n * 0.8 + 1.5)))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.6], wspace=0.5)
    left = gs[0, 0].subgridspec(2, 2, width_ratios=[1.0, 0.06],
                                height_ratios=[0.22, 1.0], wspace=0.05,
                                hspace=0.02)

    ax_dendro = fig.add_subplot(left[0, 0])
    ax = fig.add_subplot(left[1, 0])
    cax = fig.add_subplot(left[:, 1])
    ax2 = fig.add_subplot(gs[0, 1])

    # Plot dendrogram if linkage was computed; hide axis on failure.
    try:
        dendrogram(Z, ax=ax_dendro, orientation="top", no_labels=True,
                   color_threshold=None)
        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])
    except Exception:
        ax_dendro.set_visible(False)

    if title:
        fig.suptitle(title, fontsize=11, y=0.995)

    # Left: mean rho heatmap
    im = ax.imshow(rho_mean, cmap="RdYlGn", vmin=-1, vmax=1)
    plt.colorbar(im, cax=cax, label="Spearman rho (mean over states)")
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(sample_names, fontsize=8)
    for i in range(n):
        for j in range(n):
            val = rho_mean[i, j]
            color = "white" if abs(val) > 0.7 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6, color=color)
    for tick, col in zip(ax.get_xticklabels(),
                         [cmap_cond[c] for c in conditions]):
        tick.set_color(col)
    for tick, col in zip(ax.get_yticklabels(),
                         [cmap_cond[c] for c in conditions]):
        tick.set_color(col)
    ax.set_title("Mean Spearman rho\n(over all states)", fontsize=10)

    # Right: per-state rho for each pair, shown as a heatmap
    pairs = list(combinations(range(n), 2))
    pair_order = {pair: idx for idx, pair in enumerate(pairs)}
    sample_pos = {idx: pos for pos, idx in enumerate(order)}
    pairs = sorted(
        pairs,
        key=lambda ij: (
            min(sample_pos[ij[0]], sample_pos[ij[1]]),
            max(sample_pos[ij[0]], sample_pos[ij[1]]),
            pair_order[ij],
        ),
    )
    pair_labels = [sample_names[i] + " vs " + sample_names[j]
                   for i, j in pairs]
    per_state_vals = np.array([
        [rho_per_state[i, j, si] for si in range(num_states)]
        for i, j in pairs
    ])

    im2 = ax2.imshow(per_state_vals, cmap="RdYlGn", vmin=-1, vmax=1,
                     aspect="auto")
    plt.colorbar(im2, ax=ax2, label="Spearman rho")
    ax2.set_xticks(range(num_states))
    ax2.set_xticklabels(state_names, rotation=45, ha="right", fontsize=7)
    ax2.set_yticks(range(len(pairs)))
    ax2.set_yticklabels(pair_labels, fontsize=7)
    ax2.set_title("Per-state Spearman rho\n(each row = one sample pair)",
                  fontsize=10)

    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Spearman plot: {out_path}")


def plot_total_saliency_spearman(rho_mat, sample_names, conditions, out_path, title=None):
    """Plot pairwise Spearman rho heatmap for total-saliency vectors."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dist = np.clip(1.0 - rho_mat, 0.0, 2.0)
    order = clustered_sample_order(dist)
    rho_mat = rho_mat[np.ix_(order, order)]
    sample_names = [sample_names[i] for i in order]
    conditions = [conditions[i] for i in order]

    cmap_cond = condition_color_map(conditions)
    n = len(sample_names)

    fig, ax = plt.subplots(figsize=(max(6, n * 0.85 + 2), max(5, n * 0.85 + 1.5)))
    im = ax.imshow(rho_mat, cmap="RdYlGn", vmin=-1, vmax=1, aspect="equal")
    plt.colorbar(im, ax=ax, label="Spearman rho (total saliency)")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(sample_names, fontsize=9)

    for i in range(n):
        for j in range(n):
            val = float(rho_mat[i, j])
            color = "white" if abs(val) > 0.7 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color=color)

    for tick, col in zip(ax.get_xticklabels(), [cmap_cond[c] for c in conditions]):
        tick.set_color(col)
    for tick, col in zip(ax.get_yticklabels(), [cmap_cond[c] for c in conditions]):
        tick.set_color(col)

    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=k) for k, c in cmap_cond.items()]
    ax.legend(handles=handles, title="Condition",
              bbox_to_anchor=(1.25, 1), loc="upper left", fontsize=8)

    ax.set_title(
        title if title else "Pairwise Spearman rho on total saliency\n(sum over states per bin)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Total-saliency Spearman plot: {out_path}")


def plot_pca(score_mats, sample_names, conditions, replicates, out_path,
             title=None, categories=None, assay_arrows=True, top_assays=5,
             contributors_tsv=None):
    """
    Stack all samples as rows of a flat score matrix, run PCA,
    plot PC1 vs PC2 colored by condition, shaped by replicate.
    """
    from sklearn.decomposition import PCA
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Each sample: flatten the (bins, states) matrix to a 1D vector
    # Subsample bins if too many (keep up to 200k bins for speed)
    n_bins = score_mats[0].shape[0]
    max_bins = 200_000
    if n_bins > max_bins:
        rng = np.random.default_rng(42)
        idx = rng.choice(n_bins, size=max_bins, replace=False)
        idx.sort()
        print(f"  PCA: subsampling {max_bins:,} of {n_bins:,} bins")
        X = np.stack([m[idx].ravel() for m in score_mats])
    else:
        X = np.stack([m.ravel() for m in score_mats])

    pca = PCA(n_components=min(len(score_mats), 10))
    coords = pca.fit_transform(X)   # (n_samples, n_components)
    var = pca.explained_variance_ratio_

    cmap_cond = condition_color_map(conditions)
    unique_reps = sorted(set(replicates))
    markers = ["o", "s", "^", "D", "v", "P", "*"]
    rep_marker = {r: markers[i % len(markers)] for i, r in enumerate(unique_reps)}

    fig, ax = plt.subplots(figsize=(7, 6))
    for i, name in enumerate(sample_names):
        ax.scatter(
            coords[i, 0], coords[i, 1],
            c=cmap_cond[conditions[i]],
            marker=rep_marker[replicates[i]],
            s=120, edgecolors="black", linewidths=0.6,
            zorder=3, label=None,
        )
        ax.annotate(
            name,
            (coords[i, 0], coords[i, 1]),
            textcoords="offset points",
            xytext=(6, 4),
            fontsize=8,
        )

    # Legends: condition (color) and replicate (shape)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    cond_handles = [
        Patch(color=c, label=k) for k, c in cmap_cond.items()
    ]
    rep_handles = [
        Line2D([0], [0], marker=rep_marker[r], color="gray",
               markerfacecolor="gray", markersize=8, label=f"rep {r}",
               linestyle="None")
        for r in unique_reps
    ]
    legend1 = ax.legend(
        handles=cond_handles,
        title="Condition",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        fontsize=8,
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=rep_handles,
        title="Replicate",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.52),
        borderaxespad=0.0,
        fontsize=8,
    )

    ax.set_xlabel(f"PC1 ({var[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var[1]*100:.1f}% variance)")
    ax.set_title(title if title else "PCA of genome-wide KL score profiles", fontsize=11)
    ax.axhline(0, color="lightgray", lw=0.5, zorder=0)
    ax.axvline(0, color="lightgray", lw=0.5, zorder=0)

    # Optional assay-loading summary (biplot arrows + TSV export).
    # Each feature is a (bin, assay-state) value; aggregate component loadings
    # across bins to get one vector per assay-state.
    assay_v1 = None
    assay_v2 = None
    assay_mag = None
    n_states = score_mats[0].shape[1]
    if categories is not None and pca.components_.shape[0] >= 2:
        n_states = score_mats[0].shape[1]
        if n_states > 0 and pca.components_.shape[1] % n_states == 0:
            comp1 = pca.components_[0].reshape(-1, n_states)
            comp2 = pca.components_[1].reshape(-1, n_states)
            assay_v1 = comp1.mean(axis=0)
            assay_v2 = comp2.mean(axis=0)
            assay_mag = np.hypot(assay_v1, assay_v2)

    if contributors_tsv and assay_v1 is not None and assay_v2 is not None:
        pc1_abs_rank = np.argsort(np.abs(assay_v1))[::-1]
        pc2_abs_rank = np.argsort(np.abs(assay_v2))[::-1]
        mag_rank = np.argsort(assay_mag)[::-1]
        rank_pc1_map = {int(si): int(r + 1) for r, si in enumerate(pc1_abs_rank)}
        rank_pc2_map = {int(si): int(r + 1) for r, si in enumerate(pc2_abs_rank)}
        rank_mag_map = {int(si): int(r + 1) for r, si in enumerate(mag_rank)}

        out_tsv = Path(contributors_tsv)
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(out_tsv, "w") as fh:
            fh.write(
                "assay\tstate_index\tpc1_loading\tpc2_loading\tvector_magnitude"
                "\trank_pc1_abs\trank_pc2_abs\trank_vector\n"
            )
            for si in range(n_states):
                label = categories[si][0] if si < len(categories) else f"State {si+1}"
                fh.write(
                    f"{label}\t{si+1}\t{float(assay_v1[si]):.8g}\t{float(assay_v2[si]):.8g}"
                    f"\t{float(assay_mag[si]):.8g}\t{rank_pc1_map[si]}\t{rank_pc2_map[si]}\t{rank_mag_map[si]}\n"
                )
        print(f"  PCA assay loadings: {out_tsv}")

    if assay_arrows and assay_v1 is not None and assay_v2 is not None:
            k = max(1, min(int(top_assays), n_states))
            top_idx = np.argsort(assay_mag)[::-1][:k]

            coord_span = max(
                np.max(np.abs(coords[:, 0])) if coords.size else 1.0,
                np.max(np.abs(coords[:, 1])) if coords.size else 1.0,
                1.0,
            )
            max_mag = float(np.max(assay_mag[top_idx])) if top_idx.size else 1.0
            scale = (0.45 * coord_span / max(max_mag, 1e-12))

            for si in top_idx:
                dx = float(assay_v1[si]) * scale
                dy = float(assay_v2[si]) * scale
                label = categories[si][0] if si < len(categories) else f"State {si+1}"
                color = categories[si][1] if si < len(categories) else "#333333"
                if str(color).lower() == "#ffffff":
                    color = "#888888"

                ax.annotate(
                    "",
                    xy=(dx, dy),
                    xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.3),
                    zorder=2,
                )
                ax.text(
                    dx * 1.08,
                    dy * 1.08,
                    label,
                    fontsize=7,
                    color=color,
                    ha="left" if dx >= 0 else "right",
                    va="bottom" if dy >= 0 else "top",
                    zorder=2,
                )

    fig.tight_layout(rect=[0, 0, 0.80, 1])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  PCA plot:      {out_path}")


def build_region_aligned_scores(region_str, qcat_paths):
    """
    Load qcat scores restricted to one genomic region and align bins across
    all samples.

    Returns
    -------
    common_bins     : list of (chrom, start, end)
    aligned_scores  : list of np.ndarray (n_common_bins, num_states)
    num_states      : int (minimum consistent state count across samples)
    """
    from bearing_hic_plot import parse_region, load_qcat_scores

    chrom, region_start, region_end = parse_region(region_str)
    all_bins = []
    all_scores = []
    num_states = None

    for qp in qcat_paths:
        pos, score_mat, ns = load_qcat_scores(qp, chrom, region_start, region_end)
        bins = [(chrom, int(s), int(e)) for s, e in pos]
        if len(bins) == 0:
            raise ValueError(f"No bins found in region {region_str} for {qp}")
        if num_states is None:
            num_states = ns
        else:
            num_states = min(num_states, ns)
        all_bins.append(bins)
        all_scores.append(score_mat[:, :num_states])

    common_bins, aligned_scores = align_bins(all_bins, all_scores)
    return common_bins, aligned_scores, num_states


def load_stats_tsv(path):
    """
    Read a *_stats.tsv file produced by bigwig_to_qcat --stats.

    Returns dict { category_name : { metric : float } }.
    """
    import csv
    stats = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cat = row.pop("category")
            stats[cat] = {k: float(v) for k, v in row.items()}
    return stats


def plot_stats_comparison(stats_by_sample, sample_names, conditions, out_path,
                          title_suffix=""):
    """
    Bar chart comparing per-category quality statistics across all samples.

    Layout: one subplot per metric, categories on the X-axis, one bar cluster
    per category with one bar per sample, bars colored by condition.

    Parameters
    ----------
    stats_by_sample : list of dicts  { category : { metric : float } }
                      one entry per sample (None if stats file not found)
    sample_names    : list[str]
    conditions      : list[str]  -- same order as sample_names
    out_path        : str/Path
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # Filter to samples that have stats
    valid = [(i, s) for i, s in enumerate(stats_by_sample) if s is not None]
    if not valid:
        print("  Stats comparison: no stats TSV files found, skipping.")
        return

    valid_idx = [i for i, _ in valid]
    valid_stats = [s for _, s in valid]
    valid_names = [sample_names[i] for i in valid_idx]
    valid_conds = [conditions[i] for i in valid_idx]

    cmap_cond = condition_color_map(valid_conds)

    # Union of categories (preserve order from first sample with data)
    cat_order = list(valid_stats[0].keys())
    for s in valid_stats[1:]:
        for c in s:
            if c not in cat_order:
                cat_order.append(c)

    # Union of metrics
    metric_order = list(next(iter(valid_stats[0].values())).keys())

    n_cats    = len(cat_order)
    n_samples = len(valid_names)
    n_metrics = len(metric_order)

    bar_w = 0.8 / max(n_samples, 1)
    x = np.arange(n_cats)

    fig, axes = plt.subplots(
        n_metrics, 1,
        figsize=(max(10, n_cats * n_samples * bar_w * 1.4 + 3), 4 * n_metrics),
        squeeze=False,
    )

    for mi, metric in enumerate(metric_order):
        ax = axes[mi, 0]
        for si, (name, cond, stats) in enumerate(
            zip(valid_names, valid_conds, valid_stats)
        ):
            vals = [stats.get(cat, {}).get(metric, 0.0) for cat in cat_order]
            offset = (si - (n_samples - 1) / 2) * bar_w
            ax.bar(
                x + offset, vals,
                width=bar_w * 0.9,
                color=cmap_cond[cond],
                edgecolor="none",
                alpha=0.85,
                label=name if mi == 0 else None,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(cat_order, rotation=40, ha="right", fontsize=8)
        ax.set_ylabel(metric, fontsize=9)
        ax.set_title(f"Quality stats: {metric}{title_suffix}", fontsize=10)
        ax.set_xlim(-0.6, n_cats - 0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Condition legend on top subplot
    cond_handles = [Patch(color=c, label=k) for k, c in cmap_cond.items()]
    # Sample legend (individual bars) only when conditions differ from samples
    unique_conds = sorted(set(valid_conds))
    if len(valid_names) > len(unique_conds):
        # add per-sample hatching legend
        hatches = ["", "//", "xx", "..", "++", "oo"]
        for si, (name, cond) in enumerate(zip(valid_names, valid_conds)):
            offset = (si - (n_samples - 1) / 2) * bar_w
            for mi, metric in enumerate(metric_order):
                vals = [valid_stats[si].get(cat, {}).get(metric, 0.0) for cat in cat_order]
                axes[mi, 0].bar(
                    x + offset, vals,
                    width=bar_w * 0.9,
                    color="none",
                    edgecolor="black",
                    linewidth=0.4,
                    hatch=hatches[si % len(hatches)],
                    label=None,
                )
        from matplotlib.patches import Patch as _Patch
        sample_handles = [
            _Patch(facecolor=cmap_cond[valid_conds[si]],
                   hatch=hatches[si % len(hatches)],
                   edgecolor="black", linewidth=0.4,
                   label=valid_names[si])
            for si in range(n_samples)
        ]
        axes[0, 0].legend(
            handles=sample_handles,
            title="Sample",
            bbox_to_anchor=(1.01, 1), loc="upper left",
            fontsize=7, framealpha=0.8,
        )
    else:
        axes[0, 0].legend(
            handles=cond_handles,
            title="Condition",
            bbox_to_anchor=(1.01, 1), loc="upper left",
            fontsize=8, framealpha=0.8,
        )

    if len(valid_names) < len(sample_names):
        missing = [sample_names[i] for i in range(len(sample_names))
                   if i not in valid_idx]
        fig.text(
            0.5, 0.0,
            "Note: stats TSV not found for: " + ", ".join(missing),
            ha="center", fontsize=7, color="gray",
        )

    fig.tight_layout(rect=[0, 0.02, 0.88, 1.0])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Stats comparison: {out_path}")


# ---------------------------------------------------------------------------
# Q-pair JSD: Per-sample Q-vector comparability
# ---------------------------------------------------------------------------

def compute_q_pair_jsd_matrix(per_sample_q, sample_names):
    """
    Compute pairwise JSD between per-sample Q vectors.

    Parameters
    ----------
    per_sample_q : list of numpy arrays, shape (n_tracks,) each
    sample_names : list of str, length n_samples

    Returns
    -------
    jsd_matrix   : (n_samples, n_samples) numpy array, symmetric, diagonal = 0
    pair_records : list of dicts with keys: sample_A, sample_B, jsd, status
                   where status in ["ok", "warn", "error"]
    """
    n = len(per_sample_q)
    jsd_matrix = np.zeros((n, n), dtype=np.float64)
    pair_records = []

    eps = 1e-300
    for i in range(n):
        qi = _normalize_dist(per_sample_q[i])
        for j in range(i + 1, n):
            qj = _normalize_dist(per_sample_q[j])
            # JSD using 1-row matrix trick (same as existing code)
            q_jsd = js_divergence(qi[np.newaxis, :], qj[np.newaxis, :])
            jsd_matrix[i, j] = q_jsd
            jsd_matrix[j, i] = q_jsd

            # Status thresholds
            if q_jsd < 0.05:
                status = "ok"
            elif q_jsd < 0.15:
                status = "warn"
            else:
                status = "error"

            pair_records.append({
                "sample_A": sample_names[i],
                "sample_B": sample_names[j],
                "jsd": q_jsd,
                "status": status,
            })

    return jsd_matrix, pair_records


def write_q_pair_jsd_tsv(pair_records, out_path):
    """
    Write Q-pair JSD summary TSV.

    Parameters
    ----------
    pair_records : list of dicts with keys: sample_A, sample_B, jsd, status
    out_path     : Path or str
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Sort by descending JSD
    sorted_records = sorted(pair_records, key=lambda r: -r["jsd"])

    with open(out_path, "w") as fh:
        fh.write("sample_A\tsample_B\tjsd\tstatus\n")
        for rec in sorted_records:
            fh.write(f"{rec['sample_A']}\t{rec['sample_B']}\t{rec['jsd']:.6g}\t{rec['status']}\n")

    print(f"  Q-pair JSD TSV:  {out_path}")


def plot_q_pair_jsd_heatmap(jsd_matrix, sample_names, conditions, out_path, title=None):
    """
    Plot Q-pair JSD heatmap with traffic-light colormap (green/amber/red).

    Parameters
    ----------
    jsd_matrix   : (n_samples, n_samples) numpy array
    sample_names : list of str
    conditions   : list of str
    out_path     : Path or str (PDF output)
    title        : str, optional
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    n = len(sample_names)
    cmap_cond = condition_color_map(conditions)

    # Traffic-light colormap: green (ok) -> amber (warn) -> red (error)
    # The display range is 0.0-0.20, so the control points are normalized to
    # 0-1 while preserving the 0.05 and 0.15 threshold positions.
    cmap_tl = LinearSegmentedColormap.from_list(
        "traffic_light",
        [
            (0.00, "#1a7e3a"),  # dark green (ok zone)
            (0.25, "#f0c44f"),  # amber (warn threshold at 0.05 / 0.20)
            (0.75, "#c2452f"),  # red (error threshold at 0.15 / 0.20)
            (1.00, "#5a1414"),  # dark red (deep error)
        ],
    )

    fig, ax = plt.subplots(figsize=(max(8, n * 0.9 + 2), max(7, n * 0.9 + 1.5)))
    im = ax.imshow(jsd_matrix, cmap=cmap_tl, vmin=0.0, vmax=0.20, aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(sample_names, fontsize=9)

    # Annotate cells with JSD values (use em dash on diagonal)
    for i in range(n):
        for j in range(n):
            if i == j:
                text = u"\u2014"  # em dash
            else:
                val = jsd_matrix[i, j]
                text = f"{val:.3f}"
            # Text color: white in deep green or deep red, black in amber
            cell_val = jsd_matrix[i, j] if i != j else 0.0
            if cell_val < 0.025 or cell_val > 0.125:
                text_color = "white"
            else:
                text_color = "black"
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=8, color=text_color, weight="normal")

    # Color x and y tick labels by condition
    for tick, col in zip(ax.get_xticklabels(),
                         [cmap_cond[c] for c in conditions]):
        tick.set_color(col)
    for tick, col in zip(ax.get_yticklabels(),
                         [cmap_cond[c] for c in conditions]):
        tick.set_color(col)

    # Colorbar with threshold lines
    cbar = plt.colorbar(im, ax=ax, label="JSD", pad=0.02)
    # Draw dotted line at 0.05 (warn threshold)
    cbar.ax.axhline(y=0.05, color="black", linestyle=":", linewidth=1.5, zorder=10)
    cbar.ax.text(1.5, 0.05, "0.05 warn", fontsize=8, va="center")
    # Draw solid line at 0.15 (error threshold)
    cbar.ax.axhline(y=0.15, color="black", linestyle="-", linewidth=1.5, zorder=10)
    cbar.ax.text(1.5, 0.15, "0.15 error", fontsize=8, va="center")

    # Condition legend
    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=k) for k, c in cmap_cond.items()]
    ax.legend(handles=handles, title="Condition",
              bbox_to_anchor=(1.4, 1), loc="upper left", fontsize=8)

    # Title and subtitle
    if title is None:
        title = "Per-pair JSD between sample Q vectors\n(Q-comparability check)"
    subtitle = "Green: ok (JSD<0.05), Amber: warn (0.05<=JSD<0.15), Red: error (JSD>=0.15)"
    fig.suptitle(title, fontsize=11, y=0.98)
    ax.text(0.5, -0.18, subtitle, transform=ax.transAxes,
            fontsize=8, ha="center", style="italic")

    fig.tight_layout(rect=[0, 0.05, 0.88, 0.96])
    pdf_path = Path(out_path)
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")

    # Also save as PNG
    png_path = pdf_path.with_suffix(".png")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")

    plt.close(fig)
    print(f"  Q-pair JSD heatmap (PDF): {pdf_path}")
    print(f"  Q-pair JSD heatmap (PNG): {png_path}")


def get_default_track_names(n_tracks):
    """
    Get default track names for common numbers of tracks.
    For K=6: ATAC, RNAseq+, RNAseq-, CTCF, Cohesin, H3K27ac
    """
    if n_tracks == 6:
        return ["ATAC", "RNAseq_pos", "RNAseq_neg", "CTCF", "RAD21", "H3K27ac"]
    elif n_tracks == 15:
        return [name for name, _ in ALL_CATEGORIES]
    else:
        return [f"Track_{i}" for i in range(n_tracks)]


def compute_per_track_jsd_decomposition(per_sample_q, sample_names, track_names=None):
    """
    Compute per-track JSD decomposition for all pairwise sample comparisons.

    For each pair (i, j), decomposes the total JSD into contributions from each track:
    JSD_k(i, j) = 0.5 * Q_i[k] * log2(Q_i[k] / M[k]) + 0.5 * Q_j[k] * log2(Q_j[k] / M[k])
    where M[k] = (Q_i[k] + Q_j[k]) / 2

    Parameters
    ----------
    per_sample_q : list of numpy arrays, shape (K,) each
    sample_names : list of str, length n_samples
    track_names  : list of str, length K (optional; uses defaults if None)

    Returns
    -------
    per_track_records : list of dicts with keys:
        sample_A, sample_B, jsd_total, jsd_per_track (dict keyed by track name),
        dominant_track, dominant_fraction
    """
    n = len(per_sample_q)
    if n == 0:
        return []

    n_tracks = len(per_sample_q[0])
    if track_names is None:
        track_names = get_default_track_names(n_tracks)

    eps = 1e-300
    per_track_records = []

    for i in range(n):
        qi = _normalize_dist(per_sample_q[i])
        for j in range(i + 1, n):
            qj = _normalize_dist(per_sample_q[j])
            M = 0.5 * (qi + qj)
            M_safe = np.clip(M, eps, None)

            # Compute per-track JSD
            jsd_per_track_array = np.zeros(n_tracks)
            for k in range(n_tracks):
                qi_k = np.clip(qi[k], eps, None)
                qj_k = np.clip(qj[k], eps, None)
                jsd_per_track_array[k] = (
                    0.5 * qi_k * np.log2(qi_k / M_safe[k]) +
                    0.5 * qj_k * np.log2(qj_k / M_safe[k])
                )

            jsd_total = float(np.sum(jsd_per_track_array))
            jsd_per_track_dict = {
                track_names[k]: float(jsd_per_track_array[k])
                for k in range(n_tracks)
            }

            # Find dominant track
            dominant_idx = int(np.argmax(jsd_per_track_array))
            dominant_track = track_names[dominant_idx]
            dominant_value = float(jsd_per_track_array[dominant_idx])
            if jsd_total > 0:
                dominant_fraction = dominant_value / jsd_total
            else:
                dominant_fraction = 0.0

            per_track_records.append({
                "sample_A": sample_names[i],
                "sample_B": sample_names[j],
                "jsd_total": jsd_total,
                "jsd_per_track": jsd_per_track_dict,
                "dominant_track": dominant_track,
                "dominant_fraction": dominant_fraction,
            })

    return per_track_records


def write_per_track_jsd_tsv(per_track_records, track_names, out_path):
    """
    Write per-track JSD decomposition to TSV.

    Parameters
    ----------
    per_track_records : list of dicts (output of compute_per_track_jsd_decomposition)
    track_names       : list of str (track names)
    out_path          : Path or str
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Sort by descending JSD total
    sorted_records = sorted(per_track_records, key=lambda r: -r["jsd_total"])

    with open(out_path, "w") as fh:
        # Header
        header_cols = [
            "sample_A", "sample_B", "JSD_total",
            *[f"JSD_{name}" for name in track_names],
            "dominant_track", "dominant_fraction"
        ]
        fh.write("\t".join(header_cols) + "\n")

        # Rows
        for rec in sorted_records:
            cols = [
                rec["sample_A"],
                rec["sample_B"],
                f"{rec['jsd_total']:.6g}",
                *[f"{rec['jsd_per_track'][name]:.6g}" for name in track_names],
                rec["dominant_track"],
                f"{rec['dominant_fraction']:.4f}",
            ]
            fh.write("\t".join(cols) + "\n")

    print(f"  Per-track JSD TSV: {out_path}")


def plot_per_track_jsd_heatmaps(per_sample_q, sample_names, conditions,
                                track_names, out_path):
    """
    Generate small-multiples plot: one heatmap per track showing
    per-track JSD contributions in the same N×N layout as the overall JSD.

    Parameters
    ----------
    per_sample_q : list of numpy arrays, shape (K,) each
    sample_names : list of str
    conditions   : list of str
    track_names  : list of str (track names)
    out_path     : Path or str (PDF output)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    n = len(per_sample_q)
    n_tracks = len(track_names)

    # Compute per-track JSD matrices (n x n x n_tracks)
    jsd_matrices = np.zeros((n_tracks, n, n))
    eps = 1e-300

    for i in range(n):
        qi = _normalize_dist(per_sample_q[i])
        for j in range(n):
            qj = _normalize_dist(per_sample_q[j])
            if i == j:
                jsd_matrices[:, i, j] = 0.0
            else:
                M = 0.5 * (qi + qj)
                M_safe = np.clip(M, eps, None)
                for k in range(n_tracks):
                    qi_k = np.clip(qi[k], eps, None)
                    qj_k = np.clip(qj[k], eps, None)
                    jsd_matrices[k, i, j] = (
                        0.5 * qi_k * np.log2(qi_k / M_safe[k]) +
                        0.5 * qj_k * np.log2(qj_k / M_safe[k])
                    )

    # Determine layout for subplots
    ncols = 3
    nrows = int(np.ceil(n_tracks / ncols))

    # Use normalized color scale: per-track JSD thresholds are 1/K of global thresholds
    # Global: <0.05 ok, 0.05-0.15 warn, >=0.15 error
    # Per-track (K=6): <0.0083 ok, 0.0083-0.025 warn, >=0.025 error
    per_track_warn_thresh = 0.05 / max(1, n_tracks)
    per_track_error_thresh = 0.15 / max(1, n_tracks)
    vmax_color = 0.20 / max(1, n_tracks)

    cmap_tl = LinearSegmentedColormap.from_list(
        "traffic_light",
        [
            (0.00, "#1a7e3a"),  # dark green
            (0.25, "#f0c44f"),  # amber
            (0.75, "#c2452f"),  # red
            (1.00, "#5a1414"),  # dark red
        ],
    )

    cmap_cond = condition_color_map(conditions)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 3.5 + 0.5, nrows * 3 + 0.5),
        squeeze=False
    )

    for k in range(n_tracks):
        row = k // ncols
        col = k % ncols
        ax = axes[row, col]

        jsd_mat = jsd_matrices[k, :, :]
        im = ax.imshow(jsd_mat, cmap=cmap_tl, vmin=0.0, vmax=vmax_color, aspect="equal")

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(sample_names, fontsize=8)
        ax.set_title(track_names[k], fontsize=10, weight="bold")

        # Annotate cells
        for i in range(n):
            for j in range(n):
                if i == j:
                    text = u"\u2014"
                else:
                    val = jsd_mat[i, j]
                    text = f"{val:.3f}"
                cell_val = jsd_mat[i, j] if i != j else 0.0
                if cell_val < 0.025 or cell_val > 0.125:
                    text_color = "white"
                else:
                    text_color = "black"
                ax.text(j, i, text, ha="center", va="center",
                        fontsize=7, color=text_color, weight="normal")

        # Color tick labels by condition
        for tick, col in zip(ax.get_xticklabels(),
                             [cmap_cond[c] for c in conditions]):
            tick.set_color(col)
        for tick, col in zip(ax.get_yticklabels(),
                             [cmap_cond[c] for c in conditions]):
            tick.set_color(col)

    # Hide unused subplots
    for k in range(n_tracks, nrows * ncols):
        row = k // ncols
        col = k % ncols
        axes[row, col].set_visible(False)

    # Add a colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label="JSD per track")

    fig.suptitle(
        f"Per-track JSD decomposition (K={n_tracks} tracks)\nEach heatmap shows one track's contribution",
        fontsize=12, y=0.98
    )

    subtitle = (
        f"Color scale per-track thresholds: ok<{per_track_warn_thresh:.4f}, "
        f"warn {per_track_warn_thresh:.4f}-{per_track_error_thresh:.4f}, "
        f"error>={per_track_error_thresh:.4f}\n"
        f"(Scaled from global thresholds <0.05, 0.05-0.15, >=0.15 by factor 1/{n_tracks})"
    )
    fig.text(0.5, 0.01, subtitle, ha="center", fontsize=8, style="italic", wrap=True)

    fig.tight_layout(rect=[0, 0.04, 0.91, 0.96])
    pdf_path = Path(out_path)
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Per-track JSD heatmaps (PDF): {pdf_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(sheet_path, out_dir, chroms=None, skip_diff=False, skip_pca=False,
    diff_only=False, plots_only=False,
        bw_normalize=False, bw_min_signal=None,
        hic_a=None, hic_b=None, region=None, resolution=10000,
        hic_out=None, label_a=None, label_b=None,
        loops=None, genes=None, highlights=None,
    categories=None, qcat_max=5.0, regions_file=None, workers=1,
        clip_for_ini=True,
        consensus_q=False,
        diff_min_signal=0.0,
        diff_floor_combiner="max",
        diff_signal_source="auto",
        diff_signal_track_aggregate="sum",
        diff_signal_rep_aggregate="mean",
        dump_signal_hist=False,
        spearman_nonzero_only=False,
        spearman_nonzero_mode="any",
        spearman_min_bins=2,
        jsd_nonzero_only=False,
        pca_nonzero_only=False,
        global_nonzero_mode="any",
        global_nonzero_min_bins=2,
        skip_q_pair_jsd=False,
        per_track_jsd=False, beds=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = max(1, int(workers))

    _cats = categories if categories is not None else ALL_CATEGORIES

    # --- PLOTS-ONLY MODE ---
    # Skip computation, just regenerate PDFs from cached pickle files
    if plots_only:
        print(f"\n{'='*60}")
        print(f"  compare_qcat.py [PLOTS-ONLY MODE]")
        print(f"  Regenerating PDFs from cached metrics")
        print(f"  Output: {out_dir}")
        print(f"{'='*60}\n")
        
        if regions_file is None:
            sys.exit("ERROR: --plots-only requires --regions-file")
        
        import pickle
        import sys
        from bearing_hic_plot import load_regions_file as _load_rf
        
        regions_dir = out_dir / "regions"
        region_metrics_dir = out_dir / "regions" / "metrics"
        region_metrics_dir.mkdir(parents=True, exist_ok=True)
        
        # Load regions using the same function as normal mode
        print(f"Loading regions from {regions_file}...")
        regions = _load_rf(regions_file)
        print(f"  {len(regions)} regions found")
        
        # Re-render each region's PDFs from pickle cache
        for reg in regions:
            reg_name = reg["name"]
            reg_label = reg.get("label", reg_name)
            pickle_path = region_metrics_dir / f"{reg_name}.pkl"
            
            if not pickle_path.exists():
                print(f"  [{reg_name}] skipped: no cached pickle at {pickle_path}")
                continue
            
            print(f"  [{reg_name}] rendering from {pickle_path.name}...")
            try:
                with open(pickle_path, "rb") as fh:
                    payload = pickle.load(fh)
            except Exception as e:
                print(f"    ERROR: could not load pickle: {e}")
                continue
            
            # Re-render region comparison plot
            out_fig = regions_dir / f"{reg_name}_compare.pdf"
            try:
                plot_qcat_region(
                    region_str=payload["reg_str"],
                    samples=payload.get("samples"),
                    qcat_paths=payload.get("qcat_paths"),
                    diff_pairs=payload.get("diff_pairs"),
                    out_path=out_fig,
                    genes_path=genes,
                    highlights_path=highlights,
                    qcat_max=qcat_max,
                    categories=_cats,
                    region_label=reg_label,
                )
                print(f"    ✓ {out_fig.name}")
            except Exception as e:
                print(f"    ✗ ERROR rendering {out_fig.name}: {e}")
                continue
            
            # Re-render JSD if available
            if "region_jsd" in payload:
                out_jsd = region_metrics_dir / f"{reg_name}_jsd.pdf"
                try:
                    plot_jsd_heatmap(
                        payload["region_jsd"],
                        payload.get("sample_names", []),
                        payload.get("conditions", []),
                        out_jsd,
                        title=f"Per-region JSD: {reg_label}",
                    )
                    print(f"    ✓ {out_jsd.name}")
                except Exception as e:
                    print(f"    ✗ ERROR rendering {out_jsd.name}: {e}")
            
            # Re-render Spearman if available
            if "rho_mean_r" in payload:
                out_spear = region_metrics_dir / f"{reg_name}_spearman.pdf"
                try:
                    plot_spearman(
                        payload["rho_mean_r"],
                        payload.get("rho_per_state_r"),
                        payload.get("sample_names", []),
                        payload.get("conditions", []),
                        payload.get("region_states", []),
                        out_spear,
                        categories=_cats,
                        title=f"Per-region Spearman: {reg_label}",
                    )
                    print(f"    ✓ {out_spear.name}")
                except Exception as e:
                    print(f"    ✗ ERROR rendering {out_spear.name}: {e}")
            
            # Re-render PCA if available and not skipped
            if not skip_pca and "aligned_region_scores" in payload:
                out_pca = region_metrics_dir / f"{reg_name}_pca.pdf"
                try:
                    plot_pca(
                        payload["aligned_region_scores"],
                        payload.get("sample_names", []),
                        payload.get("conditions", []),
                        payload.get("replicates", []),
                        out_pca,
                        title=f"Per-region PCA: {reg_label}",
                        categories=_cats,
                    )
                    print(f"    ✓ {out_pca.name}")
                except Exception as e:
                    print(f"    ✗ ERROR rendering {out_pca.name}: {e}")
        
        print(f"\n{'='*60}")
        print("Plots regenerated successfully")
        print(f"{'='*60}\n")
        return

    # --- NORMAL MODE (with computation) ---
    samples = load_sample_sheet(sheet_path)

    # 1a. If sheet provided BigWigs instead of qcats, convert them now.
    import bigwig_to_qcat
    if any("bw_paths" in s for s in samples):
        print("\nConverting BigWig -> qcat for sheet rows...")
        for s in samples:
            if "bw_paths" in s:
                label = s.get("sample") or "<unnamed>"
                outname = s.get("out") or (label + ".qcat.bgz" if label else "sample.qcat.bgz")
                outpath = out_dir / outname
                tbi_path = Path(str(outpath) + ".tbi")

                if outpath.exists() and tbi_path.exists():
                    print(f"  {label}: using existing {outpath.name}")
                else:
                    why = []
                    if not outpath.exists():
                        why.append("qcat missing")
                    if not tbi_path.exists():
                        why.append("index missing")
                    reason = ", ".join(why)
                    print(f"  {label}: {len(s['bw_paths'])} bigwigs -> {outpath.name}  ({reason})")
                    bigwig_to_qcat.run(
                        s["bw_paths"], outpath,
                        bigwig_to_qcat.MM10_CHROM_SIZES,
                        chroms=chroms,
                        normalize_tracks=bw_normalize,
                        min_signal=(bw_min_signal if bw_min_signal is not None else bigwig_to_qcat.MIN_SIGNAL),
                        categories=categories,
                    )
                s["qcat"] = outpath

    sample_names = [s["sample"]    for s in samples]
    conditions   = [s["condition"] for s in samples]
    replicates   = [s["replicate"] for s in samples]
    qcat_paths   = [s["qcat"]      for s in samples]
    n_samples    = len(samples)

    print(f"  {n_samples} samples found:")
    for s in samples:
        print(f"    {s['sample']:20s}  cond={s['condition']}  rep={s['replicate']}  "
              f"file={s['qcat'].name}")

    # 2. Parse all qcat files
    print("\nParsing qcat.bgz files...")
    all_bins   = []
    all_scores = []
    all_raw_by_bin = []
    num_states = None

    # Always parse raw data to support Q extraction (for Q-pair JSD and consensus Q)
    if workers > 1 and len(samples) > 1:
        parse_workers = min(workers, len(samples))
        print(f"  Parallel parse enabled ({parse_workers} worker(s))")

        def _parse_one(idx, sample):
            bins, score_mat, ns, raw_by_bin = parse_qcat_bgz(
                sample["qcat"], chroms=chroms, include_raw=True
            )
            return idx, sample["sample"], bins, score_mat, ns, raw_by_bin

        ordered = [None] * len(samples)
        with concurrent.futures.ThreadPoolExecutor(max_workers=parse_workers) as ex:
            futs = [ex.submit(_parse_one, i, s) for i, s in enumerate(samples)]
            for fut in concurrent.futures.as_completed(futs):
                idx, sample_name, bins, score_mat, ns, raw_by_bin = fut.result()
                ordered[idx] = (sample_name, bins, score_mat, ns, raw_by_bin)

        for sample_name, bins, score_mat, ns, raw_by_bin in ordered:
            print(f"  {sample_name}...")
            print(f"    {len(bins):,} bins, {ns} states")
            if num_states is None:
                num_states = ns
            elif ns != num_states:
                print(f"    WARNING: {sample_name} has {ns} states, expected {num_states}. "
                      "Will use min.")
                num_states = min(num_states, ns)
            all_bins.append(bins)
            all_scores.append(score_mat[:, :num_states])
            all_raw_by_bin.append(raw_by_bin)
    else:
        for s in samples:
            print(f"  {s['sample']}...")
            bins, score_mat, ns, raw_by_bin = parse_qcat_bgz(
                s["qcat"], chroms=chroms, include_raw=True
            )
            print(f"    {len(bins):,} bins, {ns} states")
            if num_states is None:
                num_states = ns
            elif ns != num_states:
                print(f"    WARNING: {s['sample']} has {ns} states, expected {num_states}. "
                      "Will use min.")
                num_states = min(num_states, ns)
            all_bins.append(bins)
            all_scores.append(score_mat[:, :num_states])
            all_raw_by_bin.append(raw_by_bin)

    # 2b. Extract per-sample Q vectors (early, before consensus or differential analysis)
    per_sample_q = None
    Q_consensus = None
    q_source = None

    if consensus_q:
        print("\nComputing consensus Q and rescoring bins...")
        Q_consensus, q_source, per_sample_q = compute_consensus_q(samples, all_raw_by_bin)
        if len(Q_consensus) < num_states:
            print(
                "  WARNING: consensus Q has fewer tracks than parsed states; "
                f"truncating states from {num_states} to {len(Q_consensus)}"
            )
            num_states = len(Q_consensus)
        elif len(Q_consensus) > num_states:
            Q_consensus = Q_consensus[:num_states]
            per_sample_q = [q[:num_states] for q in per_sample_q]

        print(f"  Consensus Q source: {q_source}")
        print("\n  Per-sample Q and consensus (first rows shown by track):")
        header = ["track", "consensus"] + sample_names
        print("    " + "\t".join(header))
        for ti in range(num_states):
            track_name = _cats[ti][0] if ti < len(_cats) else f"State {ti+1}"
            row = [
                track_name,
                f"{float(Q_consensus[ti]):.6g}",
            ] + [f"{float(q[ti]):.6g}" for q in per_sample_q]
            print("    " + "\t".join(row))

        rescored_all_scores = []
        for s, bins, raw_by_bin in zip(samples, all_bins, all_raw_by_bin):
            score_dict, n_rescored, n_skipped = rescore_bins_with_consensus_q(
                raw_by_bin,
                Q_consensus,
                pseudocount=1e-6,
                min_signal=0.01,
            )
            mat = np.zeros((len(bins), num_states), dtype=np.float32)
            for bi, bin_key in enumerate(bins):
                row = score_dict.get(bin_key)
                if row is not None:
                    mat[bi, :] = row[:num_states]
            rescored_all_scores.append(mat)
            print(
                f"  {s['sample']}: rescored={n_rescored:,}, skipped={n_skipped:,}"
            )
        all_scores = rescored_all_scores
    elif not skip_q_pair_jsd:
        # Extract per-sample Q vectors for Q-pair JSD analysis (non-consensus mode)
        per_sample_q = [load_q_from_cats_json(s["qcat"]) for s in samples]
        if all(q is not None for q in per_sample_q):
            n_tracks = len(per_sample_q[0])
            for i, q in enumerate(per_sample_q):
                if len(q) != n_tracks:
                    print(
                        "  WARNING: Inconsistent Q length in _cats.json files; "
                        f"expected {n_tracks}, got {len(q)} for {samples[i]['sample']}. "
                        "Q-pair JSD will be skipped."
                    )
                    per_sample_q = None
                    break
            if per_sample_q:
                per_sample_q = [q[:num_states] for q in per_sample_q]
        else:
            per_sample_q_raw = []
            missing_raw_samples = []
            for sample, raw_bin_dict in zip(samples, all_raw_by_bin):
                q = _compute_q_from_raw_bins(raw_bin_dict)
                if q is None:
                    missing_raw_samples.append(sample["sample"])
                per_sample_q_raw.append(q)
            if missing_raw_samples:
                print(
                    "  WARNING: Missing raw: fields for Q-pair JSD; "
                    "skipping Q-pair analysis. Missing in: " + ", ".join(missing_raw_samples)
                )
                per_sample_q = None
            else:
                n_tracks = len(per_sample_q_raw[0])
                for i, q in enumerate(per_sample_q_raw):
                    if len(q) != n_tracks:
                        print(
                            "  WARNING: Inconsistent raw-derived Q length; "
                            "skipping Q-pair analysis."
                        )
                        per_sample_q = None
                        break
                if per_sample_q is None and per_sample_q_raw[0] is not None:
                    per_sample_q = [q[:num_states] for q in per_sample_q_raw]

    # 3. Align bins
    print("\nAligning bins...")
    common_bins, aligned_scores = align_bins(all_bins, all_scores)
    n_bins = len(common_bins)

    # 3b. Q-pair JSD analysis (per-sample Q-vector comparability check)
    if per_sample_q is not None and not skip_q_pair_jsd:
        print("\nAnalyzing per-sample Q-vector comparability...")
        q_jsd_matrix, q_pair_records = compute_q_pair_jsd_matrix(per_sample_q, sample_names)

        # Write TSV summary
        q_jsd_tsv = out_dir / "q_pair_jsd.tsv"
        write_q_pair_jsd_tsv(q_pair_records, q_jsd_tsv)

        # Plot heatmap
        q_jsd_pdf = out_dir / "q_pair_jsd_heatmap.pdf"
        plot_q_pair_jsd_heatmap(q_jsd_matrix, sample_names, conditions, q_jsd_pdf)

        # Print summary
        n_ok = sum(1 for r in q_pair_records if r["status"] == "ok")
        n_warn = sum(1 for r in q_pair_records if r["status"] == "warn")
        n_error = sum(1 for r in q_pair_records if r["status"] == "error")
        print(f"  Q-pair JSD summary: {n_ok} ok (<0.05), {n_warn} warn (0.05-0.15), {n_error} error (>=0.15)")

        # 3c. Per-track JSD decomposition (optional)
        if per_track_jsd:
            print("\nAnalyzing per-track JSD decomposition...")
            n_tracks = len(per_sample_q[0])
            track_names = get_default_track_names(n_tracks)

            per_track_records = compute_per_track_jsd_decomposition(
                per_sample_q, sample_names, track_names
            )

            # Write per-track TSV
            per_track_tsv = out_dir / "q_pair_jsd_per_track.tsv"
            write_per_track_jsd_tsv(per_track_records, track_names, per_track_tsv)

            # Plot per-track heatmaps (small multiples)
            per_track_pdf = out_dir / "q_pair_jsd_per_track.pdf"
            plot_per_track_jsd_heatmaps(
                per_sample_q, sample_names, conditions, track_names, per_track_pdf
            )

            # Print summary of dominant tracks
            warn_or_error = [r for r in q_pair_records if r["status"] in ("warn", "error")]
            n_flagged = len(warn_or_error)

            dominant_track_counts = {}
            dominant_majority = 0
            for rec in per_track_records:
                # Find corresponding q_pair record
                for q_rec in q_pair_records:
                    if (q_rec["sample_A"] == rec["sample_A"] and
                        q_rec["sample_B"] == rec["sample_B"]):
                        if q_rec["status"] in ("warn", "error"):
                            dominant_track = rec["dominant_track"]
                            frac = rec["dominant_fraction"]
                            if frac > 0.5:
                                dominant_track_counts[dominant_track] = dominant_track_counts.get(dominant_track, 0) + 1
                                dominant_majority += 1
                        break

            print(f"  Per-track JSD summary:")
            print(f"    Flagged pairs (warn/error): {n_flagged}")
            print(f"    Pairs with single dominant track (>50%): {dominant_majority}")
            if dominant_track_counts:
                print(f"    Dominant tracks:")
                for track, count in sorted(dominant_track_counts.items(), key=lambda x: -x[1]):
                    print(f"      {track}: {count} pair(s)")

            # Sanity check: verify sum of per-track JSD matches total JSD
            print(f"\n  Sanity check: per-track JSD sum vs. total JSD:")
            max_error = 0.0
            example_pair = None
            for rec, per_track_rec in zip(q_pair_records, per_track_records):
                per_track_sum = sum(per_track_rec["jsd_per_track"].values())
                error = abs(per_track_sum - rec["jsd"])
                if error > max_error:
                    max_error = error
                    example_pair = (rec, per_track_rec, per_track_sum)

            print(f"    Max difference: {max_error:.3e}")
            if max_error < 1e-10:
                print(f"    ✓ Decomposition is correct (tolerance: 1e-10)")
            else:
                print(f"    ✗ WARNING: Decomposition error exceeds 1e-10")

            if example_pair:
                rec, per_track_rec, per_track_sum = example_pair
                print(f"\n  Example (highest-JSD pair): {rec['sample_A']} vs {rec['sample_B']}")
                print(f"    Total JSD: {rec['jsd']:.6g}")
                print(f"    Per-track breakdown:")
                for track in track_names:
                    val = per_track_rec["jsd_per_track"][track]
                    frac = 100.0 * val / rec["jsd"] if rec["jsd"] > 0 else 0
                    print(f"      {track:20s}: {val:.6g} ({frac:5.1f}%)")
                print(f"    Sum of per-track: {per_track_sum:.6g}")
                print(f"    Difference: {abs(per_track_sum - rec['jsd']):.3e}")

    # Global nonzero (union) bin set. Per Methods M.9 this is the shared bin set
    # for the PRIMARY reported Spearman correlations: the union of bins with
    # nonzero total score in at least one sample (global_nonzero_mode="any"),
    # then all bins of that shared set. A single shared set across all pairs
    # keeps heatmap off-diagonals comparable; per-pair filters do not. Built
    # whenever not diff_only so the Spearman, JSD, and PCA paths use one set.
    nonzero_mask = None
    aligned_scores_nonzero = None
    if not diff_only:
        nonzero_mask = build_global_nonzero_mask(aligned_scores, mode=global_nonzero_mode)
        kept = int(nonzero_mask.sum())
        print(
            "  Global nonzero filter (shared bin set for reported Spearman): "
            f"mode={global_nonzero_mode}, kept {kept:,}/{len(nonzero_mask):,} bins"
        )
        if kept >= max(2, int(global_nonzero_min_bins)):
            aligned_scores_nonzero = [m[nonzero_mask, :] for m in aligned_scores]
            nz_summary = out_dir / "nonzero_bin_filter_summary.tsv"
            write_nonzero_bin_summary(
                nz_summary,
                sample_names,
                aligned_scores,
                nonzero_mask,
                global_nonzero_mode,
            )
            print(f"  Nonzero-bin summary: {nz_summary}")
        else:
            print(
                "  WARNING: too few bins after global nonzero filtering "
                f"({kept} < {max(2, int(global_nonzero_min_bins))}); "
                "falling back to all aligned bins for the reported Spearman."
            )

    jsd_mat = None
    if not diff_only:
        # 4. Derive probability matrices for JSD
        prob_mats = [scores_to_prob(m) for m in aligned_scores]

        # 5. JSD heatmap
        print("\nComputing pairwise JSD...")
        jsd_mat = build_jsd_matrix(prob_mats, sample_names)
        plot_jsd_heatmap(
            jsd_mat, sample_names, conditions,
            out_dir / "jsd_heatmap.pdf",
        )
        if jsd_nonzero_only and aligned_scores_nonzero is not None:
            print("\nComputing nonzero-only JSD...")
            prob_mats_nz = [scores_to_prob(m) for m in aligned_scores_nonzero]
            jsd_mat_nz = build_jsd_matrix(prob_mats_nz, sample_names)
            plot_jsd_heatmap(
                jsd_mat_nz,
                sample_names,
                conditions,
                out_dir / "jsd_heatmap_nonzero_only.pdf",
                title=(
                    "Pairwise Jensen-Shannon Divergence "
                    f"(global nonzero bins, mode={global_nonzero_mode})"
                ),
            )

        # 6. Spearman rho
        # PRIMARY (reported): all bins of the global-nonzero shared set (M.9).
        spearman_scores = (aligned_scores_nonzero
                           if aligned_scores_nonzero is not None
                           else aligned_scores)
        spearman_set_desc = (
            f"global nonzero bins, mode={global_nonzero_mode}"
            if aligned_scores_nonzero is not None
            else "all aligned bins (fallback)"
        )
        print(f"\nComputing Spearman correlations (primary: {spearman_set_desc})...")
        rho_per_state, rho_mean = build_spearman_matrix(
            spearman_scores, sample_names, num_states,
            nonzero_mode="all", min_bins=spearman_min_bins,
        )
        plot_spearman(
            rho_mean, rho_per_state, sample_names, conditions,
            num_states, out_dir / "spearman.pdf",
            categories=_cats,
            title=f"Spearman ({spearman_set_desc})",
        )

        rho_total = build_total_saliency_spearman(
            spearman_scores,
            sample_names,
            nonzero_mode="all",
            min_bins=spearman_min_bins,
        )
        plot_total_saliency_spearman(
            rho_total,
            sample_names,
            conditions,
            out_dir / "total_saliency_spearman.pdf",
        )
        write_pairwise_matrix_tsv(
            out_dir / "total_saliency_spearman.tsv",
            sample_names,
            rho_total,
        )

        # SUPPLEMENTARY: Spearman over ALL aligned bins (incl. jointly-zero),
        # reported for completeness (M.9). Skipped when the primary already used
        # all aligned bins (fallback) to avoid duplicate output.
        if aligned_scores_nonzero is not None:
            print("\nComputing supplementary Spearman (all aligned bins)...")
            rho_ps_all, rho_mean_all = build_spearman_matrix(
                aligned_scores, sample_names, num_states,
                nonzero_mode="all", min_bins=spearman_min_bins,
            )
            plot_spearman(
                rho_mean_all, rho_ps_all, sample_names, conditions,
                num_states, out_dir / "spearman_all_aligned_bins.pdf",
                categories=_cats,
                title="Spearman (all aligned bins, incl. jointly-zero; supplementary)",
            )
            rho_total_all = build_total_saliency_spearman(
                aligned_scores, sample_names,
                nonzero_mode="all", min_bins=spearman_min_bins,
            )
            write_pairwise_matrix_tsv(
                out_dir / "total_saliency_spearman_all_aligned_bins.tsv",
                sample_names,
                rho_total_all,
            )

        if spearman_nonzero_only:
            print("\nComputing nonzero-only Spearman correlations...")
            rho_nz_per_state, rho_nz_mean, bins_used_nz, bins_total_nz = build_spearman_matrix(
                aligned_scores,
                sample_names,
                num_states,
                nonzero_mode=spearman_nonzero_mode,
                min_bins=spearman_min_bins,
                return_counts=True,
            )
            plot_spearman(
                rho_nz_mean,
                rho_nz_per_state,
                sample_names,
                conditions,
                num_states,
                out_dir / "spearman_nonzero_only.pdf",
                categories=_cats,
                title=(
                    f"Spearman (nonzero-only; mode={spearman_nonzero_mode}; "
                    f"min bins={max(2, int(spearman_min_bins))})"
                ),
            )
            diag_tsv = out_dir / "spearman_nonzero_diagnostics.tsv"
            write_spearman_nonzero_diagnostics(
                out_path=diag_tsv,
                sample_names=sample_names,
                categories=_cats,
                rho_all=rho_per_state,
                rho_nonzero=rho_nz_per_state,
                bins_used=bins_used_nz,
                bins_total=bins_total_nz,
            )
            print(f"  Spearman nonzero diagnostics: {diag_tsv}")

            nz_summary = summarize_spearman_nonzero_effect(
                rho_all_mean=rho_mean,
                rho_nonzero_mean=rho_nz_mean,
                bins_used=bins_used_nz,
                bins_total=bins_total_nz,
            )
            print(
                "  Spearman nonzero coverage: "
                f"min={nz_summary['frac_min']:.4g}, "
                f"median={nz_summary['frac_med']:.4g}, "
                f"max={nz_summary['frac_max']:.4g}"
            )
            print(
                "  Spearman nonzero impact (mean-rho matrix): "
                f"mean |delta|={nz_summary['mean_abs_delta']:.4g}, "
                f"max |delta|={nz_summary['max_abs_delta']:.4g}"
            )
            if (
                np.isfinite(nz_summary["frac_med"]) and nz_summary["frac_med"] >= 0.999
                and nz_summary["max_abs_delta"] <= 1e-6
            ):
                print(
                    "  WARNING: nonzero-only Spearman appears effectively identical to all-bin Spearman. "
                    "Consider --spearman-nonzero-mode both for stricter filtering."
                )

            rho_total_nz = build_total_saliency_spearman(
                aligned_scores,
                sample_names,
                nonzero_mode=spearman_nonzero_mode,
                min_bins=spearman_min_bins,
            )
            plot_total_saliency_spearman(
                rho_total_nz,
                sample_names,
                conditions,
                out_dir / "total_saliency_spearman_nonzero_only.pdf",
                title=(
                    f"Pairwise Spearman rho on total saliency "
                    f"(nonzero-only; mode={spearman_nonzero_mode}; "
                    f"min bins={max(2, int(spearman_min_bins))})"
                ),
            )
            write_pairwise_matrix_tsv(
                out_dir / "total_saliency_spearman_nonzero_only.tsv",
                sample_names,
                rho_total_nz,
            )

        # 7. PCA
        if not skip_pca:
            print("\nRunning PCA...")
            plot_pca(
                aligned_scores, sample_names, conditions, replicates,
                out_dir / "pca.pdf",
                categories=_cats,
                contributors_tsv=out_dir / "pca_assay_loadings.tsv",
            )
            if pca_nonzero_only and aligned_scores_nonzero is not None:
                print("\nRunning nonzero-only PCA...")
                plot_pca(
                    aligned_scores_nonzero,
                    sample_names,
                    conditions,
                    replicates,
                    out_dir / "pca_nonzero_only.pdf",
                    title=(
                        "PCA of genome-wide KL score profiles "
                        f"(global nonzero bins, mode={global_nonzero_mode})"
                    ),
                    categories=_cats,
                    contributors_tsv=out_dir / "pca_assay_loadings_nonzero.tsv",
                )
    else:
        print("\nDiff-only mode: skipping JSD, Spearman, and PCA outputs.")

    # 8. Differential epilogos
    diff_pairs = []   # (cond_A, cond_B, bgz_path) -- used later for INI
    if not skip_diff:
        print("\nComputing differential epilogos...")
        # Group samples by condition
        cond_groups = defaultdict(list)
        for i, s in enumerate(samples):
            cond_groups[s["condition"]].append(i)

        # Mean score matrix per condition
        cond_means = {}
        for cond, idxs in cond_groups.items():
            cond_means[cond] = np.mean(
                [aligned_scores[i] for i in idxs], axis=0
            )
            print("  Condition '" + cond + "': mean of " + str(len(idxs)) + " replicate(s)")

        # Absolute-signal floor (per comparison). Computed once from observed
        # signal; the resulting tested-bin set is fixed and exported so the
        # permutation null can be restricted to the identical bins.
        floor_active = (diff_min_signal is not None and float(diff_min_signal) > 0.0)
        want_hist = bool(dump_signal_hist)
        sig_by_cond = None
        if floor_active or want_hist:
            # Decide the signal source. qcat files in the canonical
            # pyGenomeTracks format carry no raw: field, so fall back to the
            # per-sample BigWigs named in the sheet.
            raw_has_signal = any(
                _bin_total_signal(a) > 0
                for rbb in all_raw_by_bin for a in rbb.values()
            ) if all_raw_by_bin else False
            have_bw = any((s.get("bw") or s.get("bw_paths")) for s in samples)
            source = diff_signal_source
            if source == "auto":
                source = "raw" if raw_has_signal else ("bigwig" if have_bw else "raw")

            if source == "bigwig":
                if not have_bw:
                    raise SystemExit(
                        "ERROR: --diff-min-signal/--dump-signal-hist with bigwig "
                        "source needs a 'bw' column in the sample sheet.")
                msg_prefix = ("Absolute-signal floor" if floor_active
                              else "Signal histogram")
                tail = ("per-condition signal from BigWigs, tracks="
                        + diff_signal_track_aggregate + ", reps="
                        + diff_signal_rep_aggregate)
                if floor_active:
                    print("  " + msg_prefix + ": keep "
                          + ("max" if diff_floor_combiner != "sum" else "sum")
                          + "(sig_A, sig_B) >= " + str(diff_min_signal)
                          + "  (" + tail + ")")
                else:
                    print("  " + msg_prefix + ": (" + tail + ")")
                sig_by_cond = condition_signal_vectors_bigwig(
                    common_bins, samples, cond_groups,
                    track_aggregate=diff_signal_track_aggregate,
                    rep_aggregate=diff_signal_rep_aggregate,
                    sheet_dir=os.path.dirname(os.path.abspath(str(sheet_path))),
                )
            else:
                if floor_active:
                    print("  Absolute-signal floor: keep "
                          + ("max" if diff_floor_combiner != "sum" else "sum")
                          + "(sig_A, sig_B) >= " + str(diff_min_signal)
                          + "  (per-condition mean total raw signal)")
                else:
                    print("  Signal histogram: (per-condition mean total raw signal)")
                sig_by_cond = condition_signal_vectors(
                    common_bins, all_raw_by_bin, cond_groups
                )
            _max_sig = max(
                (float(np.max(v)) if getattr(v, "size", 0) else 0.0)
                for v in sig_by_cond.values()
            )
            if _max_sig <= 0.0:
                print("  WARNING: all per-bin signal is 0 -- no raw: fields and "
                      "no usable 'bw' column. Add a 'bw' column to the sheet, "
                      "or run without --diff-min-signal / --dump-signal-hist.")

        if want_hist and sig_by_cond is not None:
            _dump_signal_hist(out_dir, common_bins, sig_by_cond, cond_groups,
                              diff_floor_combiner)

        def _emit_diff(cond_A, cond_B, bgz_out):
            bins_use = common_bins
            A_use = cond_means[cond_A]
            B_use = cond_means[cond_B]
            if floor_active:
                bins_use, A_use, B_use, keep = apply_diff_signal_floor(
                    common_bins, cond_means[cond_A], cond_means[cond_B],
                    sig_by_cond[cond_A], sig_by_cond[cond_B],
                    diff_min_signal, diff_floor_combiner,
                )
                print("  " + cond_A + " vs " + cond_B + ": kept "
                      + format(int(keep.sum()), ",") + " / "
                      + format(len(common_bins), ",") + " bins above floor")
            diff_qcat(bins_use, A_use, B_use, bgz_out)
            if floor_active:
                bed_path = str(bgz_out)
                if bed_path.endswith(".qcat.bgz"):
                    bed_path = bed_path[:-len(".qcat.bgz")]
                bed_path += ".tested_bins.bed"
                write_tested_bins_bed(bed_path, bins_use)
                print("  Tested bins: " + bed_path)

        # Write one diff qcat per condition pair
        cond_pairs = list(combinations(sorted(cond_groups.keys()), 2))

        if workers > 1 and len(cond_pairs) > 1:
            diff_workers = min(workers, len(cond_pairs))
            print(f"  Parallel diff generation enabled ({diff_workers} worker(s))")

            def _build_diff(pair):
                cond_A, cond_B = pair
                safe_A = cond_A.replace(" ", "_")
                safe_B = cond_B.replace(" ", "_")
                out_name = "diff_" + safe_A + "_vs_" + safe_B + ".qcat.bgz"
                bgz_out = out_dir / out_name
                _emit_diff(cond_A, cond_B, bgz_out)
                return cond_A, cond_B, bgz_out

            by_pair = {}
            with concurrent.futures.ThreadPoolExecutor(max_workers=diff_workers) as ex:
                fut_to_pair = {
                    ex.submit(_build_diff, pair): pair for pair in cond_pairs
                }
                for fut in concurrent.futures.as_completed(fut_to_pair):
                    cond_A, cond_B, bgz_out = fut.result()
                    by_pair[(cond_A, cond_B)] = bgz_out

            for cond_A, cond_B in cond_pairs:
                diff_pairs.append((cond_A, cond_B, by_pair[(cond_A, cond_B)]))
        else:
            for cond_A, cond_B in cond_pairs:
                safe_A = cond_A.replace(" ", "_")
                safe_B = cond_B.replace(" ", "_")
                out_name = "diff_" + safe_A + "_vs_" + safe_B + ".qcat.bgz"
                bgz_out = out_dir / out_name
                _emit_diff(cond_A, cond_B, bgz_out)
                diff_pairs.append((cond_A, cond_B, bgz_out))

    if diff_only:
        print("\nAll outputs written to: " + str(out_dir) + "/")
        return

    # 9. Print summary table
    print("\n" + "=" * 60)
    print("JSD summary (lower = more similar):")
    print("-" * 60)
    # Print within-condition (replicate) vs between-condition stats
    within_vals  = []
    between_vals = []
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            v = jsd_mat[i, j]
            if conditions[i] == conditions[j]:
                within_vals.append(v)
            else:
                between_vals.append(v)

    if within_vals:
        print(f"  Within-condition  (replicates): "
              f"mean={np.mean(within_vals):.4f}  "
              f"max={np.max(within_vals):.4f}")
    if between_vals:
        print(f"  Between-condition             : "
              f"mean={np.mean(between_vals):.4f}  "
              f"min={np.min(between_vals):.4f}")
    if within_vals and between_vals:
        sep = np.mean(between_vals) / (np.mean(within_vals) + 1e-9)
        print(f"  Separation ratio (between/within): {sep:.1f}x")
        if sep < 1.5:
            print("  NOTE: Low separation -- conditions may be very similar "
                  "or data quality may need review.")
    print("=" * 60)

    print("\nAll outputs written to: " + str(out_dir) + "/")

    # 9b. Cross-sample stats comparison (if *_stats.tsv files exist)
    stats_by_sample = []
    for s in samples:
        qp = s["qcat"]
        stats_tsv = Path(str(qp).replace(".qcat.bgz", "").replace(".bgz", "") + "_stats.tsv")
        if stats_tsv.exists():
            try:
                stats_by_sample.append(load_stats_tsv(stats_tsv))
            except Exception as e:
                print(f"  WARNING: could not read {stats_tsv}: {e}")
                stats_by_sample.append(None)
        else:
            stats_by_sample.append(None)

    if any(s is not None for s in stats_by_sample):
        print("\nGenerating cross-sample stats comparison plot...")
        plot_stats_comparison(
            stats_by_sample, sample_names, conditions,
            out_dir / "stats_comparison.pdf",
            title_suffix=(" (consensus Q)" if consensus_q else ""),
        )

    # 10. Write categories JSON (shared across all tracks)
    # Use the first sample qcat to infer num_states; cats JSON maps
    # state index -> [name, color] matching ALL_CATEGORIES.
    cats_json_path = out_dir / "categories.json"
    cats_inner = {}
    for si in range(num_states):
        name  = _cats[si][0] if si < len(_cats) else ("State " + str(si + 1))
        color = _cats[si][1] if si < len(_cats) else "#cccccc"
        cats_inner[str(si + 1)] = [name, color]
    import json as _json
    with open(cats_json_path, "w") as _fh:
        _json.dump({"categories": cats_inner}, _fh, indent=2)
    print("  Categories JSON: " + str(cats_json_path))

    # 11. Build qcat paths for INI tracks (optionally clipped for plotting)
    def _strip_qcat_suffix(p):
        s = str(Path(p).name)
        if s.endswith(".qcat.bgz"):
            return s[:-9]
        if s.endswith(".bgz"):
            return s[:-4]
        return Path(s).stem

    if clip_for_ini:
        clip_dir = out_dir / "clipped_for_ini"
        clip_dir.mkdir(parents=True, exist_ok=True)

        qcat_paths_for_ini = []
        sample_clip_tasks = []
        for qp in qcat_paths:
            base = _strip_qcat_suffix(qp)
            clipped_qp = clip_dir / (base + f".clip{qcat_max:g}.qcat.bgz")
            sample_clip_tasks.append((qp, clipped_qp))

        if workers > 1 and len(sample_clip_tasks) > 1:
            clip_workers = min(workers, len(sample_clip_tasks))
            print(f"  Parallel clipping (sample tracks) enabled ({clip_workers} worker(s))")

            def _clip_sample(task):
                src, dst = task
                return clip_qcat_for_plot(src, dst, qcat_max=qcat_max, signed=False)

            with concurrent.futures.ThreadPoolExecutor(max_workers=clip_workers) as ex:
                qcat_paths_for_ini = list(ex.map(_clip_sample, sample_clip_tasks))
        else:
            for qp, clipped_qp in sample_clip_tasks:
                qcat_paths_for_ini.append(
                    clip_qcat_for_plot(qp, clipped_qp, qcat_max=qcat_max, signed=False)
                )

        diff_pairs_for_ini = []
        diff_clip_tasks = []
        for cond_A, cond_B, diff_bgz in diff_pairs:
            base = _strip_qcat_suffix(diff_bgz)
            clipped_diff = clip_dir / (base + f".clip{qcat_max:g}.qcat.bgz")
            diff_clip_tasks.append((cond_A, cond_B, diff_bgz, clipped_diff))

        if workers > 1 and len(diff_clip_tasks) > 1:
            dclip_workers = min(workers, len(diff_clip_tasks))
            print(f"  Parallel clipping (diff tracks) enabled ({dclip_workers} worker(s))")

            def _clip_diff(task):
                cond_A, cond_B, src, dst = task
                clipped = clip_qcat_for_plot(src, dst, qcat_max=qcat_max, signed=True)
                return cond_A, cond_B, clipped

            with concurrent.futures.ThreadPoolExecutor(max_workers=dclip_workers) as ex:
                diff_pairs_for_ini = list(ex.map(_clip_diff, diff_clip_tasks))
        else:
            for cond_A, cond_B, diff_bgz, clipped_diff in diff_clip_tasks:
                diff_pairs_for_ini.append(
                    (
                        cond_A,
                        cond_B,
                        clip_qcat_for_plot(diff_bgz, clipped_diff, qcat_max=qcat_max, signed=True),
                    )
                )
    else:
        print("  INI clipping disabled; using original qcat files in comparison_tracks.ini")
        qcat_paths_for_ini = list(qcat_paths)
        diff_pairs_for_ini = [(cond_A, cond_B, diff_bgz) for cond_A, cond_B, diff_bgz in diff_pairs]

    # 11b. Matplotlib region comparison plots (--regions-file)
    if regions_file:
        print("\nGenerating matplotlib region comparison plots...")
        from bearing_hic_plot import load_regions_file as _load_rf
        regions_list = _load_rf(regions_file)
        regions_dir = out_dir / "regions"
        regions_dir.mkdir(parents=True, exist_ok=True)
        print(f"  {len(regions_list)} region(s) from {regions_file}")

        # Store per-region metric pickles and metric PDFs under 'regions/metrics'
        region_metrics_dir = out_dir / "regions" / "metrics"
        region_metrics_dir.mkdir(parents=True, exist_ok=True)
        region_summaries = []

        def _pairwise_group_stats(mat):
            """Return within/between means and separation ratio for a pairwise matrix."""
            within_vals = []
            between_vals = []
            for i in range(n_samples):
                for j in range(i + 1, n_samples):
                    v = float(mat[i, j])
                    if conditions[i] == conditions[j]:
                        within_vals.append(v)
                    else:
                        between_vals.append(v)
            within_mean = float(np.mean(within_vals)) if within_vals else np.nan
            between_mean = float(np.mean(between_vals)) if between_vals else np.nan
            if np.isfinite(within_mean) and within_mean != 0 and np.isfinite(between_mean):
                sep_ratio = float(between_mean / (within_mean + 1e-9))
            else:
                sep_ratio = np.nan
            return within_mean, between_mean, sep_ratio

        def _compute_region_metrics(reg_name, reg_str, reg_label):
            common_bins, aligned_region_scores, region_states = build_region_aligned_scores(
                reg_str, qcat_paths
            )
            if len(common_bins) < 2:
                raise ValueError(
                    f"Too few common bins ({len(common_bins)}) for PCA/Spearman in {reg_name}"
                )

            rho_per_state_r, rho_mean_r = build_spearman_matrix(
                aligned_region_scores, sample_names, region_states
            )
            region_prob_mats = [scores_to_prob(m) for m in aligned_region_scores]
            region_jsd = build_jsd_matrix(region_prob_mats, sample_names)
            spear_within_all, spear_between_all, spear_sep_all = _pairwise_group_stats(rho_mean_r)
            jsd_within_all, jsd_between_all, jsd_sep_all = _pairwise_group_stats(region_jsd)

            summary = {
                "region": reg_name,
                "region_label": reg_label,
                "region_str": reg_str,
                "bins_total": int(len(common_bins)),
                "bins_nonzero": "",
                "fraction_bins_nonzero": "",
                "spearman_within_all": spear_within_all,
                "spearman_between_all": spear_between_all,
                "spearman_sep_between_over_within_all": spear_sep_all,
                "spearman_within_nonzero": "",
                "spearman_between_nonzero": "",
                "spearman_sep_between_over_within_nonzero": "",
                "jsd_within_all": jsd_within_all,
                "jsd_between_all": jsd_between_all,
                "jsd_sep_between_over_within_all": jsd_sep_all,
                "jsd_within_nonzero": "",
                "jsd_between_nonzero": "",
                "jsd_sep_between_over_within_nonzero": "",
            }

            payload = {
                "reg_name": reg_name,
                "reg_str": reg_str,
                "reg_label": reg_label,
                "common_bins": common_bins,
                "region_states": region_states,
                "aligned_region_scores": aligned_region_scores,
                "rho_per_state_r": rho_per_state_r,
                "rho_mean_r": rho_mean_r,
                "region_jsd": region_jsd,
                "summary": summary,
                "region_nonzero_mask": None,
                "aligned_region_scores_nz": None,
                "region_jsd_nz": None,
                "rho_nz_per_state_r": None,
                "rho_nz_mean_r": None,
                "bins_used_r_nz": None,
                "bins_total_r_nz": None,
                "nonzero_kept": None,
                "nonzero_min_keep": None,
                "nonzero_mode": global_nonzero_mode,
                "nonzero_warning": None,
            }

            # Optional region-level nonzero-only sensitivity analyses.
            if spearman_nonzero_only or jsd_nonzero_only or pca_nonzero_only:
                region_nonzero_mask = build_global_nonzero_mask(
                    aligned_region_scores,
                    mode=global_nonzero_mode,
                )
                kept = int(region_nonzero_mask.sum())
                min_keep = max(2, int(global_nonzero_min_bins))
                payload["region_nonzero_mask"] = region_nonzero_mask
                payload["nonzero_kept"] = kept
                payload["nonzero_min_keep"] = min_keep

                if kept >= min_keep:
                    aligned_region_scores_nz = [m[region_nonzero_mask, :] for m in aligned_region_scores]
                    payload["aligned_region_scores_nz"] = aligned_region_scores_nz
                    summary["bins_nonzero"] = kept
                    summary["fraction_bins_nonzero"] = kept / max(1, len(common_bins))

                    if jsd_nonzero_only:
                        region_prob_mats_nz = [scores_to_prob(m) for m in aligned_region_scores_nz]
                        region_jsd_nz = build_jsd_matrix(region_prob_mats_nz, sample_names)
                        payload["region_jsd_nz"] = region_jsd_nz
                        jsd_within_nz, jsd_between_nz, jsd_sep_nz = _pairwise_group_stats(region_jsd_nz)
                        summary["jsd_within_nonzero"] = jsd_within_nz
                        summary["jsd_between_nonzero"] = jsd_between_nz
                        summary["jsd_sep_between_over_within_nonzero"] = jsd_sep_nz

                    if spearman_nonzero_only:
                        rho_nz_per_state_r, rho_nz_mean_r, bins_used_r_nz, bins_total_r_nz = build_spearman_matrix(
                            aligned_region_scores,
                            sample_names,
                            region_states,
                            nonzero_mode=spearman_nonzero_mode,
                            min_bins=spearman_min_bins,
                            return_counts=True,
                        )
                        payload["rho_nz_per_state_r"] = rho_nz_per_state_r
                        payload["rho_nz_mean_r"] = rho_nz_mean_r
                        payload["bins_used_r_nz"] = bins_used_r_nz
                        payload["bins_total_r_nz"] = bins_total_r_nz
                        spear_within_nz, spear_between_nz, spear_sep_nz = _pairwise_group_stats(rho_nz_mean_r)
                        summary["spearman_within_nonzero"] = spear_within_nz
                        summary["spearman_between_nonzero"] = spear_between_nz
                        summary["spearman_sep_between_over_within_nonzero"] = spear_sep_nz
                else:
                    payload["nonzero_warning"] = (
                        f"    [{reg_name}] WARNING: too few bins after nonzero filtering "
                        f"({kept} < {min_keep}); skipping region nonzero-only outputs."
                    )
                    summary["bins_nonzero"] = kept
                    summary["fraction_bins_nonzero"] = kept / max(1, len(common_bins))

            return payload

        def _render_region_metrics(payload):
            reg_name = payload["reg_name"]
            reg_str = payload["reg_str"]
            reg_label = payload["reg_label"]

            out_fig = regions_dir / f"{reg_name}_compare.pdf"
            out_pca = region_metrics_dir / f"{reg_name}_pca.pdf"
            out_spear = region_metrics_dir / f"{reg_name}_spearman.pdf"
            out_jsd = region_metrics_dir / f"{reg_name}_jsd.pdf"
            out_pca_assay = region_metrics_dir / f"{reg_name}_pca_assay_loadings.tsv"

            plot_qcat_region(
                region_str=reg_str,
                samples=samples,
                qcat_paths=qcat_paths,
                diff_pairs=diff_pairs,
                out_path=out_fig,
                genes_path=genes,
                highlights_path=highlights,
                qcat_max=qcat_max,
                categories=_cats,
                region_label=reg_label,
            )

            plot_jsd_heatmap(
                payload["region_jsd"],
                sample_names,
                conditions,
                out_jsd,
                title=f"Per-region JSD: {reg_label}",
            )
            plot_spearman(
                payload["rho_mean_r"],
                payload["rho_per_state_r"],
                sample_names,
                conditions,
                payload["region_states"],
                out_spear,
                categories=_cats,
                title=f"Per-region Spearman: {reg_label}",
            )
            if not skip_pca:
                plot_pca(
                    payload["aligned_region_scores"],
                    sample_names,
                    conditions,
                    replicates,
                    out_pca,
                    title=f"Per-region PCA: {reg_label}",
                    categories=_cats,
                    contributors_tsv=out_pca_assay,
                )

            if payload["region_nonzero_mask"] is not None:
                kept = int(payload["nonzero_kept"])
                print(
                    f"    [{reg_name}] region nonzero filter "
                    f"mode={payload['nonzero_mode']}, kept {kept}/{len(payload['region_nonzero_mask'])} bins"
                )
                if payload["aligned_region_scores_nz"] is not None:
                    nz_summary = region_metrics_dir / f"{reg_name}_nonzero_bin_filter_summary.tsv"
                    write_nonzero_bin_summary(
                        nz_summary,
                        sample_names,
                        payload["aligned_region_scores"],
                        payload["region_nonzero_mask"],
                        payload["nonzero_mode"],
                    )

                    if jsd_nonzero_only and payload["region_jsd_nz"] is not None:
                        plot_jsd_heatmap(
                            payload["region_jsd_nz"],
                            sample_names,
                            conditions,
                            region_metrics_dir / f"{reg_name}_jsd_nonzero_only.pdf",
                            title=(
                                f"Per-region JSD (nonzero-only): {reg_label} "
                                f"(mode={payload['nonzero_mode']})"
                            ),
                        )

                    if spearman_nonzero_only and payload["rho_nz_mean_r"] is not None:
                        plot_spearman(
                            payload["rho_nz_mean_r"],
                            payload["rho_nz_per_state_r"],
                            sample_names,
                            conditions,
                            payload["region_states"],
                            region_metrics_dir / f"{reg_name}_spearman_nonzero_only.pdf",
                            categories=_cats,
                            title=(
                                f"Per-region Spearman (nonzero-only): {reg_label} "
                                f"(mode={spearman_nonzero_mode}; min bins={max(2, int(spearman_min_bins))})"
                            ),
                        )
                        write_spearman_nonzero_diagnostics(
                            out_path=region_metrics_dir / f"{reg_name}_spearman_nonzero_diagnostics.tsv",
                            sample_names=sample_names,
                            categories=_cats,
                            rho_all=payload["rho_per_state_r"],
                            rho_nonzero=payload["rho_nz_per_state_r"],
                            bins_used=payload["bins_used_r_nz"],
                            bins_total=payload["bins_total_r_nz"],
                        )

                    if pca_nonzero_only and not skip_pca:
                        plot_pca(
                            payload["aligned_region_scores_nz"],
                            sample_names,
                            conditions,
                            replicates,
                            region_metrics_dir / f"{reg_name}_pca_nonzero_only.pdf",
                            title=(
                                f"Per-region PCA (nonzero-only): {reg_label} "
                                f"(mode={payload['nonzero_mode']})"
                            ),
                            categories=_cats,
                            contributors_tsv=region_metrics_dir / f"{reg_name}_pca_assay_loadings_nonzero.tsv",
                        )
                elif payload["nonzero_warning"]:
                    print(payload["nonzero_warning"])

            return payload["summary"]

        if workers > 1 and len(regions_list) > 1:
            compute_workers = min(workers, len(regions_list))
            print(
                "  Parallel region metric computation enabled "
                f"({compute_workers} worker(s)); plotting remains serial for stability."
            )
            with concurrent.futures.ThreadPoolExecutor(max_workers=compute_workers) as ex:
                futures = []
                for reg in regions_list:
                    reg_name = reg["name"]
                    reg_str = reg["region"]
                    reg_label = reg["label"] or reg_name
                    futures.append((reg_name, reg_str, ex.submit(
                        _compute_region_metrics, reg_name, reg_str, reg_label
                    )))

                for reg_name, reg_str, fut in futures:
                    print(f"  [{reg_name}] {reg_str}...")
                    try:
                        payload = fut.result()
                        # Save pickle for --plots-only mode (include metadata for rendering)
                        import pickle
                        payload["samples"] = samples
                        payload["qcat_paths"] = qcat_paths
                        payload["sample_names"] = sample_names
                        payload["conditions"] = conditions
                        payload["replicates"] = replicates
                        payload["diff_pairs"] = diff_pairs
                        pickle_path = region_metrics_dir / f"{reg_name}.pkl"
                        with open(pickle_path, "wb") as fh:
                            pickle.dump(payload, fh)
                        summary = _render_region_metrics(payload)
                        region_summaries.append(summary)
                    except Exception as e:
                        print(f"    WARNING: failed for {reg_name}: {e}")
        else:
            for reg in regions_list:
                reg_name = reg["name"]
                reg_str = reg["region"]
                reg_label = reg["label"] or reg_name
                print(f"  [{reg_name}] {reg_str}...")
                try:
                    payload = _compute_region_metrics(reg_name, reg_str, reg_label)
                    # Save pickle for --plots-only mode (include metadata for rendering)
                    import pickle
                    payload["samples"] = samples
                    payload["qcat_paths"] = qcat_paths
                    payload["sample_names"] = sample_names
                    payload["conditions"] = conditions
                    payload["replicates"] = replicates
                    payload["diff_pairs"] = diff_pairs
                    pickle_path = region_metrics_dir / f"{reg_name}.pkl"
                    with open(pickle_path, "wb") as fh:
                        pickle.dump(payload, fh)
                    summary = _render_region_metrics(payload)
                    region_summaries.append(summary)
                except Exception as e:
                    print(f"    WARNING: failed for {reg_name}: {e}")

        if region_summaries:
            import csv
            summary_tsv = region_metrics_dir / "region_nonzero_comparison.tsv"
            fields = [
                "region",
                "region_label",
                "region_str",
                "bins_total",
                "bins_nonzero",
                "fraction_bins_nonzero",
                "spearman_within_all",
                "spearman_between_all",
                "spearman_sep_between_over_within_all",
                "spearman_within_nonzero",
                "spearman_between_nonzero",
                "spearman_sep_between_over_within_nonzero",
                "jsd_within_all",
                "jsd_between_all",
                "jsd_sep_between_over_within_all",
                "jsd_within_nonzero",
                "jsd_between_nonzero",
                "jsd_sep_between_over_within_nonzero",
            ]
            with open(summary_tsv, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for row in region_summaries:
                    writer.writerow(row)
            print(f"  Region nonzero comparison TSV: {summary_tsv}")

            def _to_float_or_nan(v):
                try:
                    if v == "" or v is None:
                        return np.nan
                    return float(v)
                except Exception:
                    return np.nan

            shift_rows = []
            for row in region_summaries:
                s_all = _to_float_or_nan(row.get("spearman_sep_between_over_within_all"))
                s_nz = _to_float_or_nan(row.get("spearman_sep_between_over_within_nonzero"))
                j_all = _to_float_or_nan(row.get("jsd_sep_between_over_within_all"))
                j_nz = _to_float_or_nan(row.get("jsd_sep_between_over_within_nonzero"))
                f_nz = _to_float_or_nan(row.get("fraction_bins_nonzero"))

                d_s = s_nz - s_all if np.isfinite(s_all) and np.isfinite(s_nz) else np.nan
                d_j = j_nz - j_all if np.isfinite(j_all) and np.isfinite(j_nz) else np.nan

                abs_d_s = abs(d_s) if np.isfinite(d_s) else np.nan
                abs_d_j = abs(d_j) if np.isfinite(d_j) else np.nan
                if np.isfinite(abs_d_s) and np.isfinite(abs_d_j):
                    max_abs = max(abs_d_s, abs_d_j)
                elif np.isfinite(abs_d_s):
                    max_abs = abs_d_s
                elif np.isfinite(abs_d_j):
                    max_abs = abs_d_j
                else:
                    max_abs = np.nan

                shift_rows.append({
                    "region": row.get("region", ""),
                    "region_label": row.get("region_label", ""),
                    "region_str": row.get("region_str", ""),
                    "bins_total": row.get("bins_total", ""),
                    "bins_nonzero": row.get("bins_nonzero", ""),
                    "fraction_bins_nonzero": f_nz,
                    "spearman_sep_all": s_all,
                    "spearman_sep_nonzero": s_nz,
                    "delta_spearman_sep": d_s,
                    "abs_delta_spearman_sep": abs_d_s,
                    "jsd_sep_all": j_all,
                    "jsd_sep_nonzero": j_nz,
                    "delta_jsd_sep": d_j,
                    "abs_delta_jsd_sep": abs_d_j,
                    "max_abs_sep_shift": max_abs,
                })

            shift_rows.sort(
                key=lambda r: (
                    -float(r["max_abs_sep_shift"]) if np.isfinite(r["max_abs_sep_shift"]) else float("inf"),
                    -float(r["abs_delta_spearman_sep"]) if np.isfinite(r["abs_delta_spearman_sep"]) else float("inf"),
                    -float(r["abs_delta_jsd_sep"]) if np.isfinite(r["abs_delta_jsd_sep"]) else float("inf"),
                    str(r.get("region", "")),
                )
            )

            shift_tsv = region_metrics_dir / "region_nonzero_shift_ranking.tsv"
            shift_fields = [
                "region",
                "region_label",
                "region_str",
                "bins_total",
                "bins_nonzero",
                "fraction_bins_nonzero",
                "spearman_sep_all",
                "spearman_sep_nonzero",
                "delta_spearman_sep",
                "abs_delta_spearman_sep",
                "jsd_sep_all",
                "jsd_sep_nonzero",
                "delta_jsd_sep",
                "abs_delta_jsd_sep",
                "max_abs_sep_shift",
            ]
            with open(shift_tsv, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=shift_fields, delimiter="\t")
                writer.writeheader()
                for row in shift_rows:
                    writer.writerow(row)
            print(f"  Region nonzero shift ranking TSV: {shift_tsv}")

        print(f"  Region plots: {regions_dir}/")

    # 12. Write comparison INI
    ini_path = out_dir / "comparison_tracks.ini"
    write_compare_ini(
        ini_path=ini_path,
        out_dir=out_dir,
        samples=samples,
        conditions=conditions,
        qcat_paths=qcat_paths_for_ini,
        diff_pairs=diff_pairs_for_ini,
        cats_json_path=cats_json_path,
        num_states=num_states,
        genes_path=genes,
        beds=beds,
    )

    # 13. Optional Hi-C figure when both Hi-C files and a region are provided
    if hic_a and hic_b and region and len(set(conditions)) == 2:
        condA, condB = sorted(set(conditions))
        # choose first qcat belonging to each condition
        qcatA = next(s["qcat"] for s in samples if s["condition"] == condA)
        qcatB = next(s["qcat"] for s in samples if s["condition"] == condB)
        outfig = Path(hic_out) if hic_out else (out_dir / "hic_plot.pdf")
        print(f"\nGenerating Hi-C figure {outfig} for {condA} vs {condB}...")
        from bearing_hic_plot import make_figure
        make_figure(
            hic_a, hic_b,
            qcatA, qcatB,
            region, resolution,
            outfig,
            loops_path=loops, genes_path=genes,
            highlights_path=highlights,
            label_a=label_a or condA,
            label_b=label_b or condB,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare multiple qcat.bgz epilogos files across conditions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--sheet", required=True, metavar="TSV",
        help="Sample sheet TSV with columns: sample, condition, replicate, qcat",
    )
    parser.add_argument(
        "--out", required=True, metavar="DIR",
        help="Output directory (created if needed).",
    )
    parser.add_argument(
        "--chroms", nargs="+", metavar="CHR", default=None,
        help="Restrict analysis to specific chromosomes (e.g. chr1 chr2). "
             "Useful for quick testing.",
    )
    parser.add_argument(
        "--skip-diff", action="store_true",
        help="Skip differential epilogos qcat.bgz output.",
    )
    parser.add_argument(
        "--skip-pca", action="store_true",
        help="Skip PCA plot.",
    )
    parser.add_argument(
        "--diff-only", action="store_true",
        help=(
            "Generate only differential qcat outputs. "
            "Skips JSD/Spearman/PCA plots, INI generation, clipping, "
            "region plots, and optional Hi-C figure outputs."
        ),
    )
    parser.add_argument(
        "--diff-min-signal", type=float, default=0.0, metavar="FLOOR",
        help=(
            "Absolute-signal floor for differential bins. A bin is written to "
            "the diff qcat (and therefore tested downstream) only when its "
            "combined per-condition total raw signal clears FLOOR. Default 0.0 "
            "(off). Per-condition signal is the sum of the raw per-state vector, "
            "averaged across replicates. Setting this above the noise floor keeps "
            "near-empty bins (e.g. one antisense read at a pseudogene) out of the "
            "score set and shrinks the BH-FDR multiple-testing pool. A "
            "*.tested_bins.bed sidecar is written per comparison so the "
            "permutation null can be restricted to the identical bin set. "
            "Requires raw: fields in the qcat files."
        ),
    )
    parser.add_argument(
        "--diff-floor-combiner", choices=["max", "sum"], default="max",
        help=(
            "How to combine the two conditions' signal for --diff-min-signal: "
            "'max' keeps bins where max(sig_A, sig_B) >= FLOOR (retains genuine "
            "on/off events present in only one condition); 'sum' uses "
            "sig_A + sig_B. Default: max."
        ),
    )
    parser.add_argument(
        "--diff-signal-source", choices=["auto", "raw", "bigwig"], default="auto",
        help=(
            "Where --diff-min-signal reads per-bin magnitude. 'raw' uses the "
            "qcat raw: field (absent in canonical pyGenomeTracks-compatible "
            "qcat). 'bigwig' reads each sample's BigWigs named in the sheet 'bw' "
            "column, aligned to the comparison bins. 'auto' (default) uses raw "
            "if present, else bigwig."
        ),
    )
    parser.add_argument(
        "--diff-signal-track-aggregate", choices=["sum", "max"], default="sum",
        help=(
            "BigWig source: combine the per-assay signal within a sample. 'sum' "
            "= total activity across assays; 'max' = strongest single assay. "
            "Default: sum."
        ),
    )
    parser.add_argument(
        "--diff-signal-rep-aggregate", choices=["mean", "max", "min"], default="mean",
        help=(
            "BigWig source: combine replicate signal within a condition. 'mean' "
            "matches how condition scores are averaged; 'max' keeps a bin if any "
            "replicate has signal; 'min' requires all. Default: mean."
        ),
    )
    parser.add_argument(
        "--dump-signal-hist", action="store_true",
        help=(
            "Dump the genome-wide distribution of per-condition signal (and the "
            "max/sum combined signal per condition pair) to help pick "
            "--diff-min-signal from data rather than by guess. Writes "
            "signal_hist.tsv (percentiles + fixed-threshold counts per "
            "condition pair) and signal_per_bin.tsv (per-bin values, large) to "
            "the --out directory. Honors --diff-signal-source / "
            "--diff-signal-track-aggregate / --diff-signal-rep-aggregate / "
            "--diff-floor-combiner so the numbers match what a floor run would "
            "see. Works with --diff-min-signal 0."
        ),
    )
    # options relevant when the sheet contains BigWig paths instead of qcat
    parser.add_argument(
        "--bw-normalize", action="store_true",
        help="Quantile-normalize BigWig tracks when generating qcat files.",
    )
    parser.add_argument(
        "--bw-min-signal", type=float, default=None,
        help="Min-signal threshold passed to bigwig_to_qcat (bins below are zeroed).",
    )
    # optional Hi-C figure flags (only used when exactly two conditions exist)
    parser.add_argument("--hic-a", metavar="FILE",
                        help="Hi-C contact file for condition A (cool/.hic).")
    parser.add_argument("--hic-b", metavar="FILE",
                        help="Hi-C contact file for condition B (cool/.hic).")
    parser.add_argument("--region", metavar="CHR:START-END",
                        help="Genomic region to visualize in Hi-C plot.")
    parser.add_argument("--resolution", type=int, default=10000,
                        help="Hi-C bin resolution for figure (default 10000).")
    parser.add_argument("--hic-out", metavar="FILE", default=None,
                        help="Output PDF for Hi-C figure (defaults to out_dir/hic_plot.pdf)")
    parser.add_argument("--label-a", default=None,
                        help="Label for Hi-C condition A (overrides condition name)")
    parser.add_argument("--label-b", default=None,
                        help="Label for Hi-C condition B (overrides condition name)")
    parser.add_argument("--loops", metavar="FILE",
                        help="Loop anchors BEDPE for Hi-C figure.")
    parser.add_argument("--genes", metavar="FILE",
                        help="Gene annotation BED or GTF for Hi-C figure and optional INI BED track.")
    parser.add_argument(
        "--bed", "--beds", dest="beds", metavar="BED", action="append", default=None,
        help=(
            "Additional BED files to include as BED tracks in the generated INI. "
            "Can be provided multiple times: --bed file1.bed --bed file2.bed"
        ),
    )
    parser.add_argument("--highlights", metavar="FILE",
                        help="Highlights BED4 (col4=color) for Hi-C figure.")
    parser.add_argument(
        "--categories", metavar="YAML",
        help=(
            "YAML file defining category names and colors. "
            "Overrides the built-in 15-state mm10 definitions. "
            "Example files are in categories/."
        ),
    )
    parser.add_argument(
        "--qcat-max", type=float, default=5.0,
        help=(
            "Max absolute qcat score used for plotting in generated INI tracks. "
            "Because epilogos tracks do not support min/max fields, plotting-only "
            "clipped qcat copies are written under out_dir/clipped_for_ini/. "
            "Sample tracks use [0, qcat_max]; diff tracks use [-qcat_max, qcat_max]."
        ),
    )
    parser.add_argument(
        "--no-clip", action="store_true",
        help=(
            "Disable only the plotting-only clip_qcat_for_plot step used to build "
            "comparison INI tracks. When set, the INI uses original qcat values "
            "directly. This does NOT affect KL score computation/clamping in "
            "bigwig_to_qcat.py."
        ),
    )
    parser.add_argument(
        "--consensus-q", action="store_true",
        help=(
            "If set, compute element-wise median Q across all samples per assay "
            "track and rescore all bins before comparison. Requires either "
            "_cats.json files beside each qcat file, or raw: fields in the qcat "
            "files. Recommended when per-sample Q vectors diverge (high inter-sample JSD)."
        ),
    )
    parser.add_argument(
        "--compare-consensus-q", dest="consensus_q", action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--regions-file", metavar="TSV",
        help=(
            "TSV file (columns: name, region [label]) listing genomic regions "
            "for per-region matplotlib comparison plots. "
            "One PDF per region is written to out_dir/regions/. "
            "All sample tracks share y-limit = qcat_max; diff tracks \u00b1qcat_max. "
            "Per-region PCA is skipped when --skip-pca is set."
        ),
    )
    parser.add_argument(
        "--workers", type=int, default=1, metavar="N",
        help=(
            "Number of parallel workers for independent workflow steps "
            "(parse, diff generation, clipping, region plotting). "
            "Default: 1 (serial)."
        ),
    )
    parser.add_argument(
        "--spearman-nonzero-only", action="store_true",
        help=(
            "Also compute a second Spearman output using only bins where at least "
            "one/both samples in each pair are non-zero for the given state "
            "(controlled by --spearman-nonzero-mode). "
            "Writes spearman_nonzero_only.pdf and spearman_nonzero_diagnostics.tsv."
        ),
    )
    parser.add_argument(
        "--spearman-nonzero-mode", choices=["any", "both", "all"], default="any",
        help=(
            "Bin filter mode for --spearman-nonzero-only: "
            "'any' keeps bins where at least one sample in a pair is non-zero; "
            "'both' keeps bins where both samples are non-zero; "
            "'all' keeps all bins (useful for debugging/no-op checks). "
            "Default: any."
        ),
    )
    parser.add_argument(
        "--spearman-min-bins", type=int, default=2, metavar="N",
        help=(
            "Minimum number of bins required for pairwise/state Spearman in "
            "--spearman-nonzero-only mode. Default: 2."
        ),
    )
    parser.add_argument(
        "--jsd-nonzero-only", action="store_true",
        help=(
            "Also compute JSD heatmap on globally filtered nonzero bins. "
            "Writes jsd_heatmap_nonzero_only.pdf."
        ),
    )
    parser.add_argument(
        "--pca-nonzero-only", action="store_true",
        help=(
            "Also compute PCA on globally filtered nonzero bins. "
            "Writes pca_nonzero_only.pdf and pca_assay_loadings_nonzero.tsv."
        ),
    )
    parser.add_argument(
        "--global-nonzero-mode", choices=["any", "all"], default="any",
        help=(
            "Global nonzero filter mode for nonzero-only comparisons "
            "(--jsd-nonzero-only, --pca-nonzero-only): "
            "'any' keeps bins active in at least one sample; "
            "'all' keeps bins active in every sample. Default: any."
        ),
    )
    parser.add_argument(
        "--global-nonzero-min-bins", type=int, default=2, metavar="N",
        help=(
            "Minimum bins required after global nonzero filtering to run nonzero-only "
            "JSD/PCA. Default: 2."
        ),
    )
    parser.add_argument(
        "--skip-q-pair-jsd", action="store_true",
        help=(
            "Skip Q-pair JSD computation (per-sample Q-vector comparability check). "
            "Default: enabled."
        ),
    )
    parser.add_argument(
        "--per-track-jsd", action="store_true",
        help=(
            "Compute per-track JSD decomposition. Generates additional TSV and PDF "
            "showing how much each track contributes to the total Q-vector JSD between "
            "sample pairs. Useful for identifying which assays drive compatibility issues. "
            "Default: disabled (off by default to avoid extra computation)."
        ),
    )
    parser.add_argument(
        "--plots-only", action="store_true",
        help=(
            "Regenerate PDFs from cached region metrics pickle files only. "
            "Skips all expensive computation (parsing, diff, JSD, Spearman, PCA). "
            "Requires a previous full run to have generated the pickle cache. "
            "Useful for quick PDF regeneration with different plot settings."
        ),
    )
    args = parser.parse_args()

    # Load custom categories if provided
    cli_categories = None
    if args.categories:
        cli_categories, _ = load_categories_yaml(args.categories)
        print(f"Categories loaded from: {args.categories}  ({len(cli_categories)} states)")

    print(f"\n{'='*60}")
    print(f"  compare_qcat.py")
    print(f"  Sheet:  {args.sheet}")
    print(f"  Output: {args.out}")
    print(f"  Workers: {max(1, args.workers)}")
    if args.chroms:
        print(f"  Chroms: {', '.join(args.chroms)}")
    if args.skip_diff:
        print(f"  Differential epilogos: SKIPPED")
    if args.skip_pca:
        print(f"  PCA: SKIPPED")
    if args.diff_only:
        print("  Mode: DIFF-ONLY (no compare plots/INI/clipping)")
    if args.consensus_q:
        print("  Consensus Q rescoring: ENABLED")
    if args.diff_min_signal and args.diff_min_signal > 0:
        print(f"  Diff signal floor: {args.diff_floor_combiner}(sig_A, sig_B) "
              f">= {args.diff_min_signal}")
    if args.spearman_nonzero_only:
        print("  Spearman (nonzero-only): ENABLED")
        print(f"  Spearman nonzero mode: {args.spearman_nonzero_mode}")
        print(f"  Spearman min bins: {max(2, args.spearman_min_bins)}")
    if args.jsd_nonzero_only:
        print("  JSD (nonzero-only): ENABLED")
    if args.pca_nonzero_only:
        print("  PCA (nonzero-only): ENABLED")
    if args.jsd_nonzero_only or args.pca_nonzero_only:
        print(f"  Global nonzero mode: {args.global_nonzero_mode}")
        print(f"  Global nonzero min bins: {max(2, args.global_nonzero_min_bins)}")
    if args.skip_q_pair_jsd:
        print("  Q-pair JSD: SKIPPED")
    print(f"{'='*60}\n")

    run(
        sheet_path=args.sheet,
        out_dir=args.out,
        chroms=args.chroms,
        skip_diff=args.skip_diff,
        skip_pca=args.skip_pca,
        diff_only=args.diff_only,
        plots_only=args.plots_only,
        bw_normalize=args.bw_normalize,
        bw_min_signal=args.bw_min_signal,
        hic_a=args.hic_a,
        hic_b=args.hic_b,
        region=args.region,
        resolution=args.resolution,
        hic_out=args.hic_out,
        label_a=args.label_a,
        label_b=args.label_b,
        loops=args.loops,
        genes=args.genes,
        highlights=args.highlights,
        categories=cli_categories,
        qcat_max=args.qcat_max,
        regions_file=args.regions_file,
        workers=args.workers,
        clip_for_ini=(not args.no_clip),
        consensus_q=args.consensus_q,
        diff_min_signal=args.diff_min_signal,
        diff_floor_combiner=args.diff_floor_combiner,
        diff_signal_source=args.diff_signal_source,
        diff_signal_track_aggregate=args.diff_signal_track_aggregate,
        diff_signal_rep_aggregate=args.diff_signal_rep_aggregate,
        dump_signal_hist=args.dump_signal_hist,
        spearman_nonzero_only=args.spearman_nonzero_only,
        spearman_nonzero_mode=args.spearman_nonzero_mode,
        spearman_min_bins=args.spearman_min_bins,
        jsd_nonzero_only=args.jsd_nonzero_only,
        pca_nonzero_only=args.pca_nonzero_only,
        global_nonzero_mode=args.global_nonzero_mode,
        global_nonzero_min_bins=args.global_nonzero_min_bins,
        skip_q_pair_jsd=args.skip_q_pair_jsd,
        per_track_jsd=args.per_track_jsd,
    )


if __name__ == "__main__":
    main()
