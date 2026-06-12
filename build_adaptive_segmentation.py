#!/usr/bin/env python3
"""
build_adaptive_segmentation.py -- build ONE genome-wide consensus equal-coverage
segmentation (the production version of the single-locus prototype).

Phase 1 of opt-in adaptive binning. Produces a BED of variable-width bins in which
every bin carries comparable pooled coverage, derived ONCE from pooled signal
(all samples, all tracks). This single consensus grid is then used by all samples
and by the permutation null, so the differential A - B stays defined.

Method (mirrors prototype_adaptive_binning.py, applied per chromosome):
  - fine grid at --fine-bp; pooled coverage = sum over all samples and tracks
    (abs() on negative-strand tracks), with the same loaders as the scorer
  - drop blacklisted fine bins and bins below --min-signal in EVERY sample
    (dead everywhere); a bin alive in any sample is kept
  - global coverage quota = (total retained coverage) / --target-bins; each
    chromosome gets round(chrom_coverage / quota) equal-coverage bins, gap-aware
    so a bin never spans dropped/dead space

Output BED columns: chrom  start  end  (optionally width and pooled coverage with
--with-stats). Bins are within-chromosome, non-overlapping, sorted.

  python build_adaptive_segmentation.py \
    --sheet workflow/results/samples.bearing.tsv \
    --categories categories/mm10_6track_panel.yaml \
    --chrom-sizes workflow/resources/mm10.chrom.sizes \
    --blacklist workflow/resources/mm10-blacklist.v2.bed \
    --fine-bp 50 --min-signal 0.1 --target-bins 2000000 \
    --out workflow/results/adaptive_segmentation.bed

ASCII only.
"""
import argparse
import os
import sys

import numpy as np

from bigwig_to_qcat import (
    bins_for_chrom,
    mean_signal_in_bins,
    load_blacklist,
    bins_overlapping_blacklist,
    NEGATIVE_STRAND_STATES,
)


def read_sheet(path):
    """Return list of (sample_name, [bw_paths]) from a BEARING sample sheet TSV."""
    sheet_dir = os.path.dirname(os.path.abspath(path))
    out = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        si, bi = idx.get("sample", 0), idx.get("bw", 1)
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) <= bi:
                continue
            paths = [p if os.path.isabs(p) else os.path.join(sheet_dir, p)
                     for p in (q.strip() for q in f[bi].split(",")) if p]
            out.append((f[si], paths))
    return out


def equal_coverage_segments(coverage_fine, fine_edges, n_bins, keep_mask):
    """
    ~equal-coverage segments over RETAINED fine bins, breaking at gaps. Canonical
    implementation; prototype_adaptive_binning.py mirrors this.
    Returns list of (start_bp, end_bp).
    """
    cov = np.asarray(coverage_fine, dtype=np.float64)
    keep_idx = np.flatnonzero(np.asarray(keep_mask, dtype=bool))
    if keep_idx.size == 0 or n_bins < 1:
        return []
    r_cov = cov[keep_idx]
    total = r_cov.sum()
    if total <= 0:
        bnd = np.linspace(0, len(keep_idx), n_bins + 1).round().astype(int)
    else:
        cum = np.concatenate([[0.0], np.cumsum(r_cov)])
        targets = np.linspace(0.0, total, n_bins + 1)
        bnd = np.searchsorted(cum, targets, side="left")
        bnd[0] = 0
        bnd[-1] = len(keep_idx)
        for k in range(1, len(bnd)):
            if bnd[k] <= bnd[k - 1]:
                bnd[k] = min(bnd[k - 1] + 1, len(keep_idx))
        bnd = np.minimum(bnd, len(keep_idx))
    segs = []
    for k in range(len(bnd) - 1):
        a, b = int(bnd[k]), int(bnd[k + 1])
        if b <= a:
            continue
        block = keep_idx[a:b]
        cuts = np.flatnonzero(np.diff(block) != 1)
        starts = [0] + (cuts + 1).tolist()
        ends = (cuts + 1).tolist() + [len(block)]
        for s, e in zip(starts, ends):
            lo = int(block[s])
            hi = int(block[e - 1]) + 1
            segs.append((int(fine_edges[lo]), int(fine_edges[hi])))
    return segs


def chrom_profile(sheet, chrom, chrom_len, fine_bp, blacklist, min_signal):
    """Return (fine_edges, pooled_coverage, keep_mask) for one chromosome."""
    import pyBigWig
    neg_cols = sorted(s - 1 for s in NEGATIVE_STRAND_STATES)
    edges = list(range(0, chrom_len, fine_bp)) + [chrom_len]
    bins = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    fine_edges = np.array(edges, dtype=np.int64)
    n_fine = len(bins)

    pooled = np.zeros(n_fine, dtype=np.float64)
    per_sample_total = np.zeros((len(sheet), n_fine), dtype=np.float64)
    for si, (_name, bw_paths) in enumerate(sheet):
        samp_total = np.zeros(n_fine, dtype=np.float64)
        for ci, p in enumerate(bw_paths):
            bw = pyBigWig.open(p)
            v = np.asarray(mean_signal_in_bins(bw, chrom, bins), dtype=np.float64)
            bw.close()
            if ci in neg_cols:
                v = np.abs(v)
            samp_total += v
        per_sample_total[si] = samp_total
        pooled += samp_total

    keep = np.ones(n_fine, dtype=bool)
    if blacklist is not None:
        keep &= ~bins_overlapping_blacklist(bins, blacklist, chrom)
    if min_signal > 0:
        keep &= per_sample_total.max(axis=0) >= min_signal
    return fine_edges, pooled, keep


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--categories", required=False,
                    help="kept for interface parity; neg-strand from NEGATIVE_STRAND_STATES")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--blacklist", default=None)
    ap.add_argument("--fine-bp", type=int, default=200)
    ap.add_argument("--min-signal", type=float, default=0.0)
    ap.add_argument("--target-bins", type=int, default=2000000,
                    help="approximate genome-wide bin count (sets the coverage quota)")
    ap.add_argument("--chroms", default=None,
                    help="comma-separated subset; default = all in chrom-sizes")
    ap.add_argument("--with-stats", action="store_true",
                    help="emit width and pooled-coverage columns")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sheet = read_sheet(args.sheet)
    if not sheet:
        sys.exit("[ERROR] no samples in sheet")
    blacklist = load_blacklist(args.blacklist) if args.blacklist else None

    sizes = {}
    with open(args.chrom_sizes) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    chroms = (args.chroms.split(",") if args.chroms
              else [c for c in sizes if c in sizes])
    chroms = [c for c in chroms if c in sizes]

    # pass 1: per-chrom pooled coverage + keep mask; accumulate total coverage
    cache = {}
    total_cov = 0.0
    for chrom in chroms:
        fe, pooled, keep = chrom_profile(sheet, chrom, sizes[chrom],
                                         args.fine_bp, blacklist, args.min_signal)
        cache[chrom] = (fe, pooled, keep)
        total_cov += float(pooled[keep].sum())
        sys.stderr.write("[pass1] %s: %d fine bins, %d retained, cov=%.3g\n"
                         % (chrom, len(pooled), int(keep.sum()), float(pooled[keep].sum())))

    if total_cov <= 0:
        sys.exit("[ERROR] zero total retained coverage")
    quota = total_cov / max(args.target_bins, 1)

    # pass 2: segment each chromosome with the global quota
    n_total = 0
    widths = []
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        for chrom in chroms:
            fe, pooled, keep = cache[chrom]
            chrom_cov = float(pooled[keep].sum())
            n_bins = max(1, int(round(chrom_cov / quota))) if chrom_cov > 0 else 0
            if n_bins == 0:
                continue
            segs = equal_coverage_segments(pooled, fe, n_bins, keep)
            for (s, e) in segs:
                if args.with_stats:
                    lo = int((s - fe[0]) // args.fine_bp)
                    hi = int((e - fe[0]) // args.fine_bp)
                    cov = float(pooled[lo:hi].sum())
                    out.write("%s\t%d\t%d\t%d\t%.4f\n" % (chrom, s, e, e - s, cov))
                else:
                    out.write("%s\t%d\t%d\n" % (chrom, s, e))
                widths.append(e - s)
            n_total += len(segs)

    widths = np.array(widths) if widths else np.array([0])
    sys.stderr.write(
        "[done] %d consensus bins genome-wide; quota=%.3g coverage/bin; "
        "width median=%d, range %d-%d -> %s\n"
        % (n_total, quota, int(np.median(widths)), int(widths.min()),
           int(widths.max()), args.out))


if __name__ == "__main__":
    main()
