#!/usr/bin/env python3
"""Spatial autocorrelation length of the per-bin BEARING score, and the implied
effective number of independent tests.

Motivation: the fixed 200 bp grid runs millions of correlated tests. If the
per-bin score stays correlated over many adjacent bins, the 200 bp grid
oversamples the genome and the raw "significant bin count" is one signal smeared
across correlated bins. The decorrelation length L (lag at which the
autocorrelation falls to 1/e and to 0.5) is the spatial scale of an effectively
independent unit; genome_covered / L is an order-of-magnitude estimate of the
honest number of independent tests, to be compared against the 200 bp grid, the
adaptive segmentation (~855 bp mean width, ~2.76M bins) and the Hi-C super-bin
floor (5 kb).

Reads a BEARING stats.tsv (per-sample or diff). Uses the 'bearing_score' column.
ASCII only. Needs pandas + numpy (both in the bearing env).

Usage:
  python score_autocorrelation.py STATS.tsv [--bin-size 200] [--max-lag-bp 60000]
    [--score-col bearing_score] [--main-chroms-only]
    [--out per_chrom.tsv] [--summary-out summary.tsv]

--out          per-chromosome decorrelation lengths (one row per chromosome).
--summary-out  the single-line summary values quoted in Methods (median L(1/e),
               median L(0.5), effective independent-test counts). Writing this
               is what lets the Methods number trace to a file instead of stdout.
"""
import argparse
import sys
import numpy as np
import pandas as pd

MAIN = ["chr%d" % i for i in range(1, 20)] + ["chrX"]
REFS = [("200 bp grid", 200), ("adaptive mean (~855 bp)", 855),
        ("Hi-C super-bin floor (5 kb)", 5000)]


def acf_for_chrom(start, score, bin_size, max_lag):
    """ACF vs integer bin-lag for one chromosome, on the dense bin grid with
    NaN gaps; each lag uses only bin pairs where both ends are present."""
    idx = (start // bin_size).astype(np.int64)
    idx -= idx.min()
    n = idx.max() + 1
    if n < 10:
        return None
    dense = np.full(n, np.nan, dtype=np.float64)
    dense[idx] = score.astype(np.float64)
    present = ~np.isnan(dense)
    mu = np.nanmean(dense)
    x = dense - mu
    var = np.nanmean(x * x)
    if not np.isfinite(var) or var <= 0:
        return None
    acf = np.full(max_lag + 1, np.nan)
    acf[0] = 1.0
    for k in range(1, max_lag + 1):
        a = x[:-k]; b = x[k:]
        m = present[:-k] & present[k:]
        if m.sum() < 1000:
            break
        acf[k] = np.mean(a[m] * b[m]) / var
    return acf


def crossing(acf, level, bin_size):
    """First lag (bp) at which acf drops to <= level; None if it never does."""
    for k in range(1, len(acf)):
        if not np.isfinite(acf[k]):
            return None
        if acf[k] <= level:
            # linear interp between k-1 and k
            a0, a1 = acf[k - 1], acf[k]
            frac = (a0 - level) / (a0 - a1) if a1 != a0 else 0.0
            return (k - 1 + frac) * bin_size
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("--bin-size", type=int, default=200)
    ap.add_argument("--max-lag-bp", type=int, default=60000)
    ap.add_argument("--score-col", default="bearing_score")
    ap.add_argument("--main-chroms-only", action="store_true", default=True)
    ap.add_argument("--out", default=None,
                    help="Per-chromosome decorrelation-length table (TSV).")
    ap.add_argument("--summary-out", default=None,
                    help="One-row summary of the values quoted in Methods (TSV).")
    args = ap.parse_args()

    max_lag = args.max_lag_bp // args.bin_size
    df = pd.read_csv(args.stats, sep="\t", usecols=lambda c: c in
                     ("chrom", "start", args.score_col))
    df = df.rename(columns={args.score_col: "score"})
    df = df[np.isfinite(df["score"])]
    if args.main_chroms_only:
        df = df[df["chrom"].isin(MAIN)]

    rows = []
    genome_bp = 0
    eff_e = 0.0
    eff_half = 0.0
    for chrom, g in df.groupby("chrom"):
        g = g.sort_values("start")
        acf = acf_for_chrom(g["start"].to_numpy(), g["score"].to_numpy(),
                            args.bin_size, max_lag)
        if acf is None:
            continue
        L_e = crossing(acf, 1.0 / np.e, args.bin_size)
        L_h = crossing(acf, 0.5, args.bin_size)
        span = int(g["start"].max() - g["start"].min()) + args.bin_size
        rows.append((chrom, L_h, L_e, len(g), span))
        genome_bp += span
        if L_e:
            eff_e += span / L_e
        if L_h:
            eff_half += span / L_h

    print("=" * 78)
    print("BEARING score spatial autocorrelation  |  %s" % args.stats)
    print("-" * 78)
    print("%-8s %10s %10s %12s" % ("chrom", "L(0.5)bp", "L(1/e)bp", "n_bins"))
    Lh_all, Le_all = [], []
    for chrom, L_h, L_e, n, span in rows:
        print("%-8s %10s %10s %12d" % (
            chrom,
            "%.0f" % L_h if L_h else "NA",
            "%.0f" % L_e if L_e else "NA", n))
        if L_h:
            Lh_all.append(L_h)
        if L_e:
            Le_all.append(L_e)
    print("-" * 78)
    med_h = float(np.median(Lh_all)) if Lh_all else float("nan")
    med_e = float(np.median(Le_all)) if Le_all else float("nan")
    n_bins_total = len(df)
    print("median decorrelation length:  L(0.5) = %.0f bp   L(1/e) = %.0f bp"
          % (med_h, med_e))
    print("genome covered (main chroms): %.1f Mb across %d present bins"
          % (genome_bp / 1e6, n_bins_total))
    print()
    print("Effective independent tests (genome_covered / L):")
    print("  using L(1/e): ~%.0f  (= %.2fx fewer than %d bins on the 200bp grid)"
          % (eff_e, n_bins_total / eff_e if eff_e else float("nan"), n_bins_total))
    print("  using L(0.5): ~%.0f  (= %.2fx fewer)"
          % (eff_half, n_bins_total / eff_half if eff_half else float("nan")))
    print()
    print("Decorrelation length vs candidate testing units:")
    for name, w in REFS:
        if np.isfinite(med_e):
            rel = med_e / w
            verdict = ("grid >> independence scale (oversampled %.1fx)" % rel
                       if rel > 1.5 else
                       "grid ~ independence scale" if rel > 0.67 else
                       "grid coarser than independence scale")
            print("  %-28s : L(1/e)/width = %.2f  -> %s" % (name, rel, verdict))
    print("=" * 78)

    if args.out:
        import os
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("# score_autocorrelation.py per-chromosome table\n")
            fh.write("# stats=%s score_col=%s bin_size=%d max_lag_bp=%d\n"
                     % (args.stats, args.score_col, args.bin_size, args.max_lag_bp))
            fh.write("chrom\tL_half_bp\tL_1e_bp\tn_bins\tspan_bp\n")
            for chrom, L_h, L_e, n, span in rows:
                fh.write("%s\t%s\t%s\t%d\t%d\n"
                         % (chrom,
                            "%.1f" % L_h if L_h else "NA",
                            "%.1f" % L_e if L_e else "NA", n, span))
        print("wrote %s (%d chromosomes)" % (args.out, len(rows)))

    if args.summary_out:
        import os
        os.makedirs(os.path.dirname(args.summary_out) or ".", exist_ok=True)
        with open(args.summary_out, "w") as fh:
            fh.write("# score_autocorrelation.py summary (values quoted in Methods)\n")
            fh.write("# stats=%s score_col=%s bin_size=%d\n"
                     % (args.stats, args.score_col, args.bin_size))
            fh.write("metric\tvalue\n")
            fh.write("median_L_1e_bp\t%.1f\n" % med_e)
            fh.write("median_L_half_bp\t%.1f\n" % med_h)
            fh.write("genome_covered_bp\t%d\n" % genome_bp)
            fh.write("present_bins\t%d\n" % n_bins_total)
            fh.write("eff_indep_tests_L1e\t%.0f\n" % eff_e)
            fh.write("eff_indep_tests_Lhalf\t%.0f\n" % eff_half)
            fh.write("oversampling_factor_L1e\t%.2f\n"
                     % (n_bins_total / eff_e if eff_e else float("nan")))
        print("wrote %s" % args.summary_out)


if __name__ == "__main__":
    main()
