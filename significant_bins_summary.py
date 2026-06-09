#!/usr/bin/env python3
"""Summarize FDR-significant BEARING bins per region per sample.

Reads the per-sample p-value tables written by bearing_pvalue.py
(pvalue/<sample>.stats.tsv) and a regions file (name<TAB>region<TAB>...),
and writes one row per region with one column per sample giving the count of
FDR-significant bins (and the total bins) in that region.

This replaces the cross-comparison summary that the old single-batch figure
runner produced as a side effect; it now stands alone so it does not depend on
the (parallelized) figure jobs.

The stats TSV columns are: chrom, start, end, bearing_score, pval,
pval_adj_bh, significant_fdr<level>. Significance is taken from the
significant_fdr* column when present (truthy = significant), else from
pval_adj_bh < --fdr.

ASCII-only; no third-party dependencies.
"""
import argparse
import glob
import math
import os
import sys


def parse_region(region_str):
    # "chr6:40793981-41688054"
    chrom, span = region_str.split(":")
    start_s, end_s = span.replace(",", "").split("-")
    return chrom, int(start_s), int(end_s)


def load_regions(path):
    regions = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # tolerate either a header row or a bare data row
        cols = {h: i for i, h in enumerate(header)}
        has_header = "region" in cols and "name" in cols
        if not has_header:
            fh.seek(0)
        name_i = cols.get("name", 0)
        region_i = cols.get("region", 1)
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) <= region_i:
                continue
            name = f[name_i]
            region_str = f[region_i]
            try:
                chrom, start, end = parse_region(region_str)
            except Exception:
                print("[WARN] skipping malformed region: %s" % region_str,
                      file=sys.stderr)
                continue
            regions.append((name, region_str, chrom, start, end))
    return regions


def _norm_chrom(c):
    return c[3:] if c.startswith("chr") else c


def count_bins_in_region(stats_path, chrom, start, end, fdr):
    """Return (n_sig, n_total) for bins overlapping [start, end) on chrom."""
    n_sig = 0
    n_total = 0
    target = _norm_chrom(chrom)
    with open(stats_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        if "chrom" not in idx or "start" not in idx or "end" not in idx:
            return (0, 0)
        ci, si, ei = idx["chrom"], idx["start"], idx["end"]
        sig_col = next((h for h in header if h.lower().startswith("significant_fdr")),
                       None)
        sig_i = idx.get(sig_col) if sig_col else None
        padj_i = idx.get("pval_adj_bh")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ci, si, ei):
                continue
            if _norm_chrom(f[ci]) != target:
                continue
            try:
                bs, be = int(f[si]), int(f[ei])
            except ValueError:
                continue
            if be <= start or bs >= end:
                continue
            n_total += 1
            is_sig = False
            if sig_i is not None and sig_i < len(f):
                v = f[sig_i].strip().lower()
                is_sig = v in ("1", "true", "yes", "t")
            elif padj_i is not None and padj_i < len(f):
                try:
                    is_sig = float(f[padj_i]) < fdr
                except ValueError:
                    is_sig = False
            if is_sig:
                n_sig += 1
    return (n_sig, n_total)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pvalue-dir", required=True,
                    help="Directory with <sample>.stats.tsv files.")
    ap.add_argument("--regions", required=True,
                    help="Regions file (name<TAB>region<TAB>...).")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="FDR level for pval_adj_bh fallback (default 0.05).")
    ap.add_argument("--exclude-diff", action="store_true", default=True,
                    help="Skip diff_*.stats.tsv (per-comparison) files.")
    ap.add_argument("--include-diff", dest="exclude_diff", action="store_false",
                    help="Also summarize diff_*.stats.tsv files.")
    args = ap.parse_args()

    stats_files = sorted(glob.glob(os.path.join(args.pvalue_dir, "*.stats.tsv")))
    samples = []
    for p in stats_files:
        base = os.path.basename(p)[:-len(".stats.tsv")]
        if args.exclude_diff and base.startswith("diff_"):
            continue
        samples.append((base, p))
    if not samples:
        print("[ERROR] no per-sample .stats.tsv files in %s" % args.pvalue_dir,
              file=sys.stderr)
        sys.exit(1)

    regions = load_regions(args.regions)
    sample_names = [s for s, _ in samples]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    header = ["region_name", "region"]
    for s in sample_names:
        header.append("%s_sig" % s)
        header.append("%s_total" % s)
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for (name, region_str, chrom, start, end) in regions:
            row = [name, region_str]
            for (_s, path) in samples:
                n_sig, n_total = count_bins_in_region(path, chrom, start, end,
                                                      args.fdr)
                row.append(str(n_sig))
                row.append(str(n_total))
            out.write("\t".join(row) + "\n")

    print("Wrote %s (%d regions x %d samples; FDR significance)" % (
        args.out, len(regions), len(sample_names)))


if __name__ == "__main__":
    main()
