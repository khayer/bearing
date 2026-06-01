#!/usr/bin/env python3
"""
rebin_qcat.py
=============
Rebin a Bearing qcat.bgz file from 200 bp bins to a coarser resolution
for cleaner visualisation at wide genomic regions.

Aggregation strategy
--------------------
For each merged window of W consecutive 200 bp bins:

  score_merged(i)  = mean( s(b,i) for b in window )   per track
  S_merged         = sum_i score_merged(i)             total

Using the MEAN preserves the per-bin interpretation: a merged-bin score of
1.5 means the same as a single-bin score of 1.5 (average enrichment per
200 bp), making scores directly comparable across zoom levels.

The dominant track colour is determined by whichever track has the largest
sum of s(b,i) across the window -- this represents the track most
consistently enriched across the region, not just the highest in one bin.

USAGE
-----
  # Merge to 2 kb bins (10 x 200 bp):
  python rebin_qcat.py --qcat sample.qcat.bgz --factor 10 --out sample_2kb.qcat.bgz

  # Merge to a specific target bin size:
  python rebin_qcat.py --qcat sample.qcat.bgz --bin-size 2000 --out sample_2kb.qcat.bgz

  # Multiple outputs for different zoom levels:
  python rebin_qcat.py \\
    --qcat sample.qcat.bgz \\
    --bin-size 200 2000 10000 \\
    --out-prefix sample

  # Also rebin the companion p-value BigWig:
  python rebin_qcat.py \\
    --qcat sample.qcat.bgz \\
    --pval sample.neglog10p.bw \\
    --bin-size 2000 10000 \\
    --out-prefix sample

REBINNING GUIDE
---------------
  Region width     Suggested bin size   --factor (from 200 bp)
  < 100 kb         200 bp               1   (no rebinning needed)
  100 kb - 500 kb  1 kb                 5
  500 kb - 2 Mb    2 kb                 10
  2 Mb - 10 Mb     5 kb                 25
  10 Mb - 50 Mb    10 kb                50
  > 50 Mb          25-50 kb             125-250

OUTPUT
------
  <out>.qcat.bgz      -- rebinned qcat, tabix-indexed
  <out>.qcat.bgz.tbi  -- tabix index
  <out>.neglog10p.bw  -- rebinned p-value BigWig (if --pval given)

DEPENDENCIES
------------
  pip install pysam numpy pyBigWig
"""

import argparse
import gzip
import json
import math
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np


BASE_BIN_SIZE = 200  # bp


# ---------------------------------------------------------------------------
# Parse qcat
# ---------------------------------------------------------------------------

def parse_qcat_region(path, chrom_filter=None):
    """
    Yield (chrom, start, end, per_track_dict) from a qcat.bgz.
    per_track_dict: {state_idx (1-based int): float score}
    """
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom_filter and chrom != chrom_filter:
                continue
            meta = parts[3]
            qi = meta.find("qcat:")
            if qi < 0:
                continue
            raw_start = meta.find(",raw:", qi)
            if raw_start >= 0:
                qcat_payload = meta[qi + 5:raw_start]
            else:
                qcat_payload = meta[qi + 5:]
            pairs = json.loads(qcat_payload)
            per_track = {int(state_idx): float(score)
                         for score, state_idx in pairs}
            yield chrom, start, end, per_track


# ---------------------------------------------------------------------------
# Rebinning
# ---------------------------------------------------------------------------

def rebin_qcat(in_path, out_path, target_bin_size):
    """
    Rebin a qcat.bgz to target_bin_size using mean aggregation.
    Writes a new qcat.bgz + .tbi index.
    """
    try:
        import pysam
    except ImportError:
        sys.exit("ERROR: pysam is required.  pip install pysam")

    factor = target_bin_size // BASE_BIN_SIZE
    if target_bin_size % BASE_BIN_SIZE != 0:
        sys.exit(
            f"ERROR: target bin size {target_bin_size} must be a multiple "
            f"of {BASE_BIN_SIZE} bp.")
    if factor < 1:
        sys.exit(f"ERROR: target bin size must be >= {BASE_BIN_SIZE} bp.")

    print(f"  Rebinning {Path(in_path).name}: "
          f"{BASE_BIN_SIZE} bp -> {target_bin_size} bp "
          f"(factor {factor}x, mean aggregation)")

    # ── Collect all bins, grouped by chromosome ───────────────────────────
    # chrom -> sorted list of (start, end, per_track_dict)
    chrom_bins = defaultdict(list)
    num_states = 0

    for chrom, start, end, per_track in parse_qcat_region(in_path):
        chrom_bins[chrom].append((start, end, per_track))
        if per_track:
            num_states = max(num_states, max(per_track.keys()))

    if not chrom_bins:
        sys.exit("ERROR: no bins parsed from input qcat.")

    print(f"    {sum(len(v) for v in chrom_bins.values()):,} input bins  "
          f"{num_states} tracks  "
          f"{len(chrom_bins)} chromosomes")

    # ── Write rebinned rows to temp file ──────────────────────────────────
    tmp_path = str(out_path) + ".tmp.tsv"
    rows_written = 0
    bin_id = 1

    with open(tmp_path, "w") as fh:
        for chrom in sorted(chrom_bins.keys(),
                             key=lambda c: (
                                 0 if c.startswith("chr") else 1,
                                 int(c.replace("chr", "")
                                      .replace("X", "98")
                                      .replace("Y", "99"))
                                 if c.replace("chr", "")
                                      .replace("X", "98")
                                      .replace("Y", "99").isdigit()
                                 else 0, c)):
            bins = sorted(chrom_bins[chrom], key=lambda x: x[0])

            # Snap all bins to target_bin_size windows
            # Group by: window_idx = start // target_bin_size
            windows = defaultdict(list)
            for start, end, pt in bins:
                win_idx = start // target_bin_size
                windows[win_idx].append((start, end, pt))

            for win_idx in sorted(windows.keys()):
                win_bins = windows[win_idx]
                win_start = win_idx * target_bin_size
                win_end   = win_start + target_bin_size

                # Mean per-track score across bins in window
                track_sums = defaultdict(float)
                n = len(win_bins)
                for _, _, pt in win_bins:
                    for state_idx, score in pt.items():
                        track_sums[state_idx] += score

                # Mean = sum / n
                track_means = {k: v / n for k, v in track_sums.items()}

                # Sort pairs descending by score for qcat format
                pairs = sorted(
                    ([round(v, 6), k] for k, v in track_means.items()),
                    key=lambda x: -x[0]
                )

                qcat_col = "id:" + str(bin_id) + ",qcat:" + json.dumps(pairs, separators=(",", ":"))
                fh.write(f"{chrom}\t{win_start}\t{win_end}\t{qcat_col}\n")
                rows_written += 1
                bin_id += 1

    print(f"    -> {rows_written:,} rebinned bins")

    # ── Sort ──────────────────────────────────────────────────────────────
    import subprocess
    sorted_path = tmp_path + ".sorted"
    result = subprocess.run(
        ["sort", "-k1,1", "-k2,2n", tmp_path],
        capture_output=True, text=True)
    if result.returncode == 0:
        with open(sorted_path, "w") as f:
            f.write(result.stdout)
    else:
        # Python fallback
        with open(tmp_path) as f:
            lines = f.readlines()
        lines.sort(key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
        with open(sorted_path, "w") as f:
            f.writelines(lines)
    os.remove(tmp_path)

    # ── Bgzip + tabix ─────────────────────────────────────────────────────
    out_str = str(out_path)
    pysam.tabix_compress(sorted_path, out_str, force=True)
    os.remove(sorted_path)
    pysam.tabix_index(out_str, preset="bed", force=True)
    print(f"    Written: {out_str}")
    print(f"             {out_str}.tbi")


def rebin_bigwig(in_path, out_path, target_bin_size, chrom_sizes=None):
    """
    Rebin a BigWig to target_bin_size using mean aggregation.
    Reads the source at BASE_BIN_SIZE resolution and averages.
    """
    try:
        import pyBigWig
    except ImportError:
        sys.exit("ERROR: pyBigWig is required.  pip install pyBigWig")

    print(f"  Rebinning BigWig {Path(in_path).name}: "
          f"{BASE_BIN_SIZE} bp -> {target_bin_size} bp")

    with pyBigWig.open(str(in_path)) as bw_in:
        chroms = bw_in.chroms()
        if chrom_sizes:
            chroms = {c: s for c, s in chroms.items()
                      if c in chrom_sizes}

        header = [(c, s) for c, s in chroms.items()]

        with pyBigWig.open(str(out_path), "w") as bw_out:
            bw_out.addHeader(header)

            for chrom, chrom_len in chroms.items():
                n_output_bins = math.ceil(chrom_len / target_bin_size)
                starts, ends, vals = [], [], []

                for i in range(n_output_bins):
                    s = i * target_bin_size
                    e = min(s + target_bin_size, chrom_len)
                    try:
                        v = bw_in.stats(chrom, s, e, type="mean")[0]
                    except Exception:
                        v = None
                    if v is None or math.isnan(v):
                        v = 0.0
                    if v != 0.0:
                        starts.append(s)
                        ends.append(e)
                        vals.append(float(v))

                if starts:
                    bw_out.addEntries(
                        [chrom] * len(starts), starts,
                        ends=ends, values=vals)

    print(f"    Written: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Rebin a Bearing qcat.bgz to coarser resolution.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--qcat", required=True, metavar="FILE",
                    help="Input qcat.bgz file.")
    ap.add_argument("--pval", default=None, metavar="FILE",
                    help="Optional: companion -log10(p) BigWig to rebin "
                         "(output of bearing_pvalue.py).")
    ap.add_argument("--out", default=None, metavar="FILE",
                    help="Output path for a single bin size. "
                         "Mutually exclusive with --out-prefix.")
    ap.add_argument("--out-prefix", default=None, metavar="STR",
                    help="Output prefix when generating multiple bin sizes. "
                         "Files are named <prefix>_<binsize>bp.qcat.bgz.")
    ap.add_argument("--bin-size", nargs="+", type=int, default=None,
                    metavar="BP",
                    help="Target bin size(s) in bp. Must be a multiple of "
                         f"{BASE_BIN_SIZE}. E.g. --bin-size 2000 10000.")
    ap.add_argument("--factor", type=int, default=None, metavar="N",
                    help=f"Merge factor (target = factor x {BASE_BIN_SIZE} bp). "
                         "Alternative to --bin-size for a single output.")
    args = ap.parse_args()

    # Resolve target bin sizes
    if args.factor and args.bin_size:
        sys.exit("ERROR: --factor and --bin-size are mutually exclusive.")
    if args.factor:
        targets = [args.factor * BASE_BIN_SIZE]
    elif args.bin_size:
        targets = args.bin_size
    else:
        sys.exit("ERROR: specify --bin-size or --factor.")

    if args.out and len(targets) > 1:
        sys.exit("ERROR: --out can only be used with a single bin size. "
                 "Use --out-prefix for multiple.")
    if args.out and args.out_prefix:
        sys.exit("ERROR: --out and --out-prefix are mutually exclusive.")
    if not args.out and not args.out_prefix:
        # Default: derive prefix from input filename
        base = str(Path(args.qcat).name)
        for ext in (".qcat.bgz", ".bgz"):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        args.out_prefix = base
        print(f"No --out or --out-prefix given; using prefix: {args.out_prefix}")

    for target_bp in targets:
        if target_bp % BASE_BIN_SIZE != 0:
            print(f"WARNING: {target_bp} is not a multiple of "
                  f"{BASE_BIN_SIZE} bp, skipping.")
            continue

        if args.out:
            qcat_out = Path(args.out)
        else:
            qcat_out = Path(f"{args.out_prefix}_{target_bp}bp.qcat.bgz")

        rebin_qcat(args.qcat, qcat_out, target_bp)

        if args.pval:
            if args.out:
                pval_out = Path(str(args.out).replace(".qcat.bgz", "")
                                .replace(".bgz", "") + ".neglog10p.bw")
            else:
                pval_out = Path(
                    f"{args.out_prefix}_{target_bp}bp.neglog10p.bw")
            rebin_bigwig(args.pval, pval_out, target_bp)

    print("\nDone.")


if __name__ == "__main__":
    main()
