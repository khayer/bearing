#!/usr/bin/env python3
"""Build a COHORT-WIDE, per-track quantile-normalization reference.

For each track (state), pool that track's genome-wide nonzero signal across ALL
samples in the run and summarize it as a fixed-length quantile grid; the
reference for a track is the average of the per-sample quantile grids. Each
sample's copy of a track is later mapped onto this shared per-track reference
(see bigwig_to_qcat.py --normalize-method cohort-quantile), which equalizes
dynamic-range differences BETWEEN samples for the same assay (e.g. a
low-dynamic-range ULI ChIP vs a deep ChIP, or ChIP vs CUT&RUN). This is the
axis that the within-sample cross-track normalization does NOT address.

The signal transform here matches the pre-normalization step in scoring:
abs() on negative-strand tracks and blacklist zeroing, computed on the same
200 bp bins. Only the per-track reference distribution is written; no scores.

Reads the BEARING sample sheet (sample<TAB>bw with comma-separated BigWig
paths, in track/state order). Writes an .npz with: ref_probs (L,), ref
(num_states, L), state_names (num_states,), bin_size.

ASCII-only.
"""
import argparse
import os
import sys
import numpy as np

# Reuse the exact loaders scoring uses, so the reference is built on the same
# bins/signal as the scored data.
from bigwig_to_qcat import (
    BIN_SIZE,
    NEGATIVE_STRAND_STATES,
    bins_for_chrom,
    bins_overlapping_blacklist,
    load_blacklist,
    mean_signal_in_bins,
    load_categories_yaml,
)


def read_sheet(path):
    """Return list of (sample_name, [bw_paths]) from a sample sheet TSV."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        s_i = idx.get("sample", 0)
        bw_i = idx.get("bw", 1)
        sheet_dir = os.path.dirname(os.path.abspath(path))
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) <= bw_i:
                continue
            name = f[s_i]
            paths = []
            for p in f[bw_i].split(","):
                p = p.strip()
                if not p:
                    continue
                paths.append(p if os.path.isabs(p) else os.path.join(sheet_dir, p))
            rows.append((name, paths))
    return rows


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    return sizes


def state_names_from_categories(cat_path, n_states):
    """Best-effort track names for labeling the reference; falls back to S1..Sn."""
    names = []
    try:
        if cat_path and cat_path.endswith((".yaml", ".yml")):
            cats = load_categories_yaml(cat_path)
            names = [c[0] if isinstance(c, (list, tuple)) else str(c) for c in cats]
        elif cat_path and cat_path.endswith(".json"):
            import json
            d = json.load(open(cat_path))
            cats = d.get("categories", d)
            # keys "1".."n" (1-indexed) -> names
            names = [None] * n_states
            for k, v in cats.items():
                try:
                    i = int(k) - 1
                except (ValueError, TypeError):
                    continue
                if 0 <= i < n_states:
                    names[i] = v[0] if isinstance(v, (list, tuple)) else str(v)
    except Exception as e:
        print("[WARN] could not read category names: %s" % e, file=sys.stderr)
    names = [n if n else "S%d" % (i + 1) for i, n in enumerate(names or [None] * n_states)]
    return names[:n_states] if len(names) >= n_states else \
        names + ["S%d" % (i + 1) for i in range(len(names), n_states)]


def sample_track_quantiles(bw_paths, chrom_sizes, blacklist, neg_states,
                           ref_probs, bin_size):
    """For one sample, return (n_states, L) array of per-track nonzero quantile
    grids over the genome (abs on neg-strand tracks, blacklist zeroed)."""
    import pyBigWig
    n_states = len(bw_paths)
    L = len(ref_probs)
    # accumulate nonzero values per track across chroms
    per_track_vals = [list() for _ in range(n_states)]
    for chrom, clen in chrom_sizes.items():
        bins = bins_for_chrom(clen, bin_size)
        if not bins:
            continue
        bl_mask = (bins_overlapping_blacklist(bins, blacklist, chrom)
                   if blacklist else None)
        for si, bw_path in enumerate(bw_paths):
            try:
                with pyBigWig.open(str(bw_path)) as bw:
                    col = mean_signal_in_bins(bw, chrom, bins)
            except Exception as e:
                print("  WARNING: %s %s: %s" % (chrom, bw_path, e), file=sys.stderr)
                continue
            if (si + 1) in neg_states:
                col = np.abs(col)
            if bl_mask is not None:
                col = col.copy()
                col[bl_mask] = 0.0
            nz = col[col > 0]
            if nz.size:
                per_track_vals[si].append(nz.astype(np.float64))
    grids = np.zeros((n_states, L), dtype=np.float64)
    for si in range(n_states):
        if per_track_vals[si]:
            allnz = np.concatenate(per_track_vals[si])
            grids[si] = np.quantile(allnz, ref_probs)
        # else: leave zeros (track absent/empty for this sample)
    return grids


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sheet", required=True, help="BEARING sample sheet TSV.")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--blacklist", default=None)
    ap.add_argument("--categories", default=None,
                    help="Category file (json/yaml) for track names only.")
    ap.add_argument("--out", required=True, help="Output .npz reference.")
    ap.add_argument("--quantile-points", type=int, default=4096,
                    help="Length of the per-track quantile grid (default 4096).")
    ap.add_argument("--bin-size", type=int, default=BIN_SIZE)
    args = ap.parse_args()

    samples = read_sheet(args.sheet)
    if not samples:
        print("[ERROR] no samples in sheet", file=sys.stderr)
        sys.exit(1)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    blacklist = load_blacklist(args.blacklist) if args.blacklist else None
    neg_states = set(NEGATIVE_STRAND_STATES)
    ref_probs = np.linspace(0.0, 1.0, args.quantile_points)

    n_states = len(samples[0][1])
    # running mean of per-track quantile grids across samples
    ref_sum = np.zeros((n_states, args.quantile_points), dtype=np.float64)
    ref_count = np.zeros(n_states, dtype=np.int64)
    for name, bw_paths in samples:
        if len(bw_paths) != n_states:
            print("[WARN] %s has %d tracks (expected %d); skipping"
                  % (name, len(bw_paths), n_states), file=sys.stderr)
            continue
        print("  scanning %s (%d tracks)" % (name, len(bw_paths)), file=sys.stderr)
        grids = sample_track_quantiles(bw_paths, chrom_sizes, blacklist,
                                       neg_states, ref_probs, args.bin_size)
        for si in range(n_states):
            if grids[si].any():
                ref_sum[si] += grids[si]
                ref_count[si] += 1
    ref = np.zeros_like(ref_sum)
    for si in range(n_states):
        if ref_count[si] > 0:
            ref[si] = ref_sum[si] / ref_count[si]

    names = state_names_from_categories(args.categories, n_states)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez_compressed(args.out, ref_probs=ref_probs, ref=ref.astype(np.float32),
                        state_names=np.array(names, dtype=object),
                        bin_size=np.int64(args.bin_size))
    print("Wrote %s: per-track cohort reference, %d states x %d quantile points, "
          "from %d samples" % (args.out, n_states, args.quantile_points, len(samples)))


if __name__ == "__main__":
    main()
