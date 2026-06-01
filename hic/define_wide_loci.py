#!/usr/bin/env python3
"""
define_wide_loci_v2.py

Empirically define wide loci (e.g. Tcrb, Igh, Igk) as contiguous
chromosomal intervals where a contrast shows statistically significant
BEARING signal (FDR < threshold) above a small magnitude floor.

v2 changes from v1:
  - Primary filter is now `significant_fdr<thresh>` (uses the
    pval_adj_bh column from BEARING output) rather than raw |BES|
    genome-wide percentile. v1 missed Tcrb in DN-vs-EbKO because raw
    |BES| saturates at a ceiling (~3.8) across many noisy loci, so
    the top 1% by magnitude was dominated by spurious peaks. FDR
    correction is the appropriate filter for distinguishing real
    locus territory from background.
  - Optional secondary |BES| floor (default 1.0) drops bins that are
    statistically significant but biologically tiny.
  - Reports per-locus fraction of bins that are significant (a real
    AR locus should be near-uniformly significant across its width;
    a spurious locus will have low coverage).

Algorithm:
  1. Compute the genome-wide BES significance threshold from the
     pval_adj_bh column (default FDR 0.05).
  2. Take bins with significant_fdr<thresh> AND |bearing_score| >= floor.
  3. Merge bins within --merge-gap of each other (default 100 kb).
  4. Retain merged intervals with >= --min-bins above-threshold bins.
  5. For each retained interval, report n_bins, max |BES|, mean |BES|,
     and fraction of bins in the interval that were above-threshold
     (a quality metric — should be > 0.5 for real loci).

Usage:
  python define_wide_loci_v2.py \\
      --bes results_v6/diff_DN_vs_EbKO.stats.tsv \\
      --fdr 0.05 --bes-floor 1.0 --merge-gap 100000 --min-bins 50 \\
      --out wide_loci_DN_vs_EbKO.bed \\
      --out-stats wide_loci_DN_vs_EbKO.stats.tsv
"""
import argparse
import sys
import pandas as pd
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bes", required=True,
                    help="BEARING diff stats.tsv with bearing_score and "
                         "pval_adj_bh columns")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="FDR threshold for significance (default 0.05)")
    ap.add_argument("--bes-floor", type=float, default=1.0,
                    help="minimum |bearing_score| for a bin to be considered "
                         "locus-defining, even if FDR-significant. Drops "
                         "biologically negligible significant bins. "
                         "Default 1.0; set 0 to disable.")
    ap.add_argument("--merge-gap", type=int, default=100000,
                    help="merge above-threshold bins separated by < this "
                         "many bp (default 100000)")
    ap.add_argument("--min-bins", type=int, default=50,
                    help="minimum above-threshold bins per merged interval "
                         "(default 50)")
    ap.add_argument("--min-sig-fraction", type=float, default=0.0,
                    help="minimum fraction of bins in the merged interval "
                         "that must be significant. 0 disables (default). "
                         "Try 0.3-0.5 for stricter locus calling.")
    ap.add_argument("--out", required=True, help="output BED path")
    ap.add_argument("--out-stats", default=None,
                    help="optional per-locus stats TSV")
    ap.add_argument("--chrom-filter", default=None,
                    help="restrict output to comma-separated chrom list "
                         "(e.g. 'chr6,chr12,chr14'). For debugging.")
    args = ap.parse_args()

    print("Loading", args.bes, "...", flush=True)
    df = pd.read_csv(args.bes, sep="\t", low_memory=False)
    needed = {"chrom", "start", "end", "bearing_score", "pval_adj_bh"}
    missing = needed - set(df.columns)
    if missing:
        sys.exit("missing columns in " + args.bes + ": " + str(missing))

    df["abs_bes"] = df["bearing_score"].abs()
    n_total = len(df)

    # Significance + magnitude floor
    sig_mask = df["pval_adj_bh"] < args.fdr
    floor_mask = df["abs_bes"] >= args.bes_floor
    keep_mask = sig_mask & floor_mask

    n_sig = int(sig_mask.sum())
    n_kept = int(keep_mask.sum())
    print("  total bins: {:,}".format(n_total))
    print("  bins with pval_adj_bh < {}: {:,} ({:.2f}%)".format(
        args.fdr, n_sig, 100 * n_sig / n_total))
    print("  ... AND |BES| >= {}: {:,} ({:.2f}%)".format(
        args.bes_floor, n_kept, 100 * n_kept / n_total))

    if n_kept == 0:
        sys.exit("no bins pass filters; check thresholds or input file")

    # Keep all df (significant + non-significant) for fraction-significant
    # calculation later. Work on sig subset for merging.
    sig = df[keep_mask].copy().sort_values(
        ["chrom", "start"]).reset_index(drop=True)
    print("  {} chroms with at least one above-threshold bin".format(
        sig["chrom"].nunique()))

    # Merge within gap; compute per-locus stats
    loci = []
    for chrom, grp in sig.groupby("chrom"):
        grp = grp.sort_values("start").reset_index(drop=True)
        if len(grp) == 0:
            continue
        cur_start = int(grp["start"].iloc[0])
        cur_end = int(grp["end"].iloc[0])
        cur_n = 1
        cur_max = float(grp["abs_bes"].iloc[0])
        cur_sum = cur_max
        for i in range(1, len(grp)):
            s = int(grp["start"].iloc[i])
            e = int(grp["end"].iloc[i])
            b = float(grp["abs_bes"].iloc[i])
            if s - cur_end <= args.merge_gap:
                cur_end = e
                cur_n += 1
                cur_max = max(cur_max, b)
                cur_sum += b
            else:
                if cur_n >= args.min_bins:
                    loci.append((chrom, cur_start, cur_end, cur_n,
                                 cur_sum / cur_n, cur_max))
                cur_start, cur_end, cur_n = s, e, 1
                cur_max = cur_sum = b
        if cur_n >= args.min_bins:
            loci.append((chrom, cur_start, cur_end, cur_n,
                         cur_sum / cur_n, cur_max))

    # Compute fraction-significant for each merged locus (total bins of
    # any kind in the merged span, vs the n_above we already have).
    # Quick chrom-indexed lookup of all df coords.
    df_by_chrom = {c: g.sort_values("start").reset_index(drop=True)
                    for c, g in df.groupby("chrom")}

    enriched = []
    for c, s, e, n_above, mean_b, max_b in loci:
        g = df_by_chrom.get(c)
        if g is None:
            frac_sig = float("nan")
            n_total_in_locus = 0
        else:
            in_locus = g[(g["start"] >= s) & (g["end"] <= e)]
            n_total_in_locus = len(in_locus)
            if n_total_in_locus:
                frac_sig = n_above / n_total_in_locus
            else:
                frac_sig = float("nan")
        enriched.append((c, s, e, n_above, n_total_in_locus,
                         frac_sig, mean_b, max_b))

    # Optional fraction-significant gate
    if args.min_sig_fraction > 0:
        before = len(enriched)
        enriched = [r for r in enriched if r[5] >= args.min_sig_fraction]
        print("  fraction-significant >= {}: kept {}/{} loci".format(
            args.min_sig_fraction, len(enriched), before))

    # Optional chrom filter (debugging)
    if args.chrom_filter:
        keep_chroms = set(args.chrom_filter.split(","))
        enriched = [r for r in enriched if r[0] in keep_chroms]

    # Write BED
    with open(args.out, "w") as fh:
        for c, s, e, n_above, n_total_in, frac_sig, mean_b, max_b in enriched:
            name = "{}:{}-{}".format(c, s, e)
            fh.write("{}\t{}\t{}\t{}\n".format(c, s, e, name))
    print("Wrote", len(enriched), "loci to", args.out)

    if args.out_stats:
        stats = pd.DataFrame(enriched, columns=[
            "chrom", "start", "end", "n_bins_sig", "n_bins_total",
            "frac_sig", "mean_abs_BES", "max_abs_BES"])
        stats["width_bp"] = stats["end"] - stats["start"]
        stats = stats[["chrom", "start", "end", "width_bp", "n_bins_sig",
                       "n_bins_total", "frac_sig", "mean_abs_BES",
                       "max_abs_BES"]]
        # rank by max |BES| primarily, then n_bins_sig, then width
        stats = stats.sort_values(
            ["max_abs_BES", "n_bins_sig", "width_bp"],
            ascending=[False, False, False])
        stats.to_csv(args.out_stats, sep="\t", index=False)
        print("Wrote per-locus stats to", args.out_stats)

    if enriched:
        print("\nTop 10 loci by n_bins_sig (i.e. most extensive signal):")
        top_by_n = sorted(enriched, key=lambda r: -r[3])[:10]
        for c, s, e, n_a, n_t, f_s, m_b, mx_b in top_by_n:
            print("  {}:{:,}-{:,}  width={:>9,} bp  n_sig={:>4}  "
                  "frac_sig={:.2f}  max|BES|={:.2f}".format(
                      c, s, e, e - s, n_a, f_s, mx_b))


if __name__ == "__main__":
    main()
