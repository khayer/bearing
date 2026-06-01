#!/usr/bin/env python3
"""
per_locus_summary.py

For each interval in a BED file and each contrast in a list of BEARING
diff_*.stats.tsv files, compute a summary row with:

  - n_bins_total            (number of 200 bp bins in the interval)
  - n_bins_sig_fdr05        (FDR<0.05 by pval_adj_bh)
  - frac_sig_fdr05
  - mean_abs_bes_locus
  - mean_abs_bes_chrom_bg   (mean |BES| on the same chrom OUTSIDE locus)
  - elevation_ratio         (mean_abs_bes_locus / mean_abs_bes_chrom_bg)
  - max_abs_bes
  - top_track_contributor   (which kl_i contributes most |signal|)
  - per_track_fractions     (kl_1..kl_6 as fraction of summed |kl|)

Produces a long-format TSV (one row per locus x contrast) suitable for
the supplementary tables. Also a wide-format pivot for the main results
table if --pivot-out is given.

Usage:
  python per_locus_summary.py \\
      --loci-bed published_loci.bed \\
      --diffs DN_vs_DP=results_v6/diff_DN_vs_DP.stats.tsv \\
              DN_vs_EbKO=results_v6/diff_DN_vs_EbKO.stats.tsv \\
              DN_vs_ProB=results_v6/diff_DN_vs_ProB.stats.tsv \\
              DN_vs_S3T3=results_v6/diff_DN_vs_S3T3.stats.tsv \\
      --out per_locus_summary.tsv \\
      --pivot-out per_locus_summary.pivot.tsv

BED format (tab-separated):
  chrom  start  end  name  [optional extra columns ignored]
"""
import argparse
import sys
import pandas as pd
import numpy as np


def load_loci(bed_path):
    """Load 4-column BED. Tolerates extra columns."""
    rows = []
    with open(bed_path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#") or ln.startswith("track"):
                continue
            parts = ln.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            name = parts[3] if len(parts) >= 4 else "{}:{}-{}".format(
                chrom, s, e)
            rows.append({"chrom": chrom, "start": s, "end": e, "name": name})
    return pd.DataFrame(rows)


def summarize_one_locus(df_diff, chrom, start, end, kl_cols):
    """Compute summary for a single (locus, contrast) pair."""
    df_chrom = df_diff[df_diff["chrom"] == chrom]
    in_locus = df_chrom[
        (df_chrom["start"] >= start) & (df_chrom["end"] <= end)
    ]
    out_locus = df_chrom[
        (df_chrom["end"] < start) | (df_chrom["start"] > end)
    ]

    n_total = len(in_locus)
    if n_total == 0:
        return None

    abs_bes_locus = in_locus["bearing_score"].abs()
    abs_bes_bg = out_locus["bearing_score"].abs()

    sig = in_locus[in_locus["pval_adj_bh"] < 0.05]
    mean_locus = float(abs_bes_locus.mean())
    mean_bg = float(abs_bes_bg.mean()) if len(out_locus) > 0 else float("nan")
    ratio = mean_locus / mean_bg if (mean_bg and np.isfinite(mean_bg)) else float("nan")

    # Per-track contribution: sum |kl_i| over the locus
    track_sums = {col: float(in_locus[col].abs().sum()) for col in kl_cols}
    total_kl = sum(track_sums.values())
    if total_kl > 0:
        track_fracs = {col: track_sums[col] / total_kl for col in kl_cols}
    else:
        track_fracs = {col: 0.0 for col in kl_cols}

    if total_kl > 0:
        top_col = max(track_sums, key=track_sums.get)
    else:
        top_col = ""

    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "width_bp": end - start,
        "n_bins_total": n_total,
        "n_bins_sig_fdr05": int(len(sig)),
        "frac_sig_fdr05": len(sig) / n_total if n_total else 0.0,
        "mean_abs_bes_locus": mean_locus,
        "mean_abs_bes_chrom_bg": mean_bg,
        "elevation_ratio": ratio,
        "max_abs_bes": float(abs_bes_locus.max()) if n_total else float("nan"),
        "top_track_contributor": top_col,
        **{"frac_" + col: track_fracs[col] for col in kl_cols},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci-bed", required=True,
                    help="BED file with regions of interest")
    ap.add_argument("--diffs", required=True, nargs="+",
                    help="Space-separated CONTRAST=PATH entries, e.g. "
                         "DN_vs_DP=results_v6/diff_DN_vs_DP.stats.tsv")
    ap.add_argument("--out", required=True,
                    help="long-format per-locus output TSV")
    ap.add_argument("--pivot-out", default=None,
                    help="optional wide-format pivot TSV "
                         "(loci x contrasts x key metrics)")
    args = ap.parse_args()

    loci = load_loci(args.loci_bed)
    print("Loaded", len(loci), "loci from", args.loci_bed, flush=True)

    contrast_paths = []
    for spec in args.diffs:
        if "=" not in spec:
            sys.exit("invalid --diffs entry (need NAME=PATH): " + spec)
        name, path = spec.split("=", 1)
        contrast_paths.append((name.strip(), path.strip()))

    rows = []
    for cname, cpath in contrast_paths:
        print("\nLoading", cpath, "as", cname, "...", flush=True)
        df = pd.read_csv(cpath, sep="\t", low_memory=False)
        kl_cols = sorted([c for c in df.columns if c.startswith("kl_")],
                          key=lambda c: int(c.split("_")[1]))
        if not kl_cols:
            sys.exit(cpath + ": no kl_* columns found")

        for _, locus in loci.iterrows():
            summary = summarize_one_locus(
                df, locus["chrom"], int(locus["start"]),
                int(locus["end"]), kl_cols)
            if summary is None:
                print("  WARNING: no bins in {} for {}".format(
                    locus["name"], cname))
                continue
            summary["locus_name"] = locus["name"]
            summary["contrast"] = cname
            rows.append(summary)

    if not rows:
        sys.exit("no summaries produced; check inputs")

    out_df = pd.DataFrame(rows)
    # Reorder columns: locus_name, contrast, then summary metrics
    front = ["locus_name", "contrast", "chrom", "start", "end", "width_bp"]
    metric_cols = [c for c in out_df.columns if c not in front]
    out_df = out_df[front + metric_cols]
    out_df.to_csv(args.out, sep="\t", index=False, float_format="%.4f")
    print("\nWrote", len(out_df), "rows to", args.out)

    if args.pivot_out:
        # Compact pivot: locus rows, contrast columns, with the most
        # important metrics combined into one cell as "ratio (n_sig)".
        pivot_rows = []
        for locus_name, grp in out_df.groupby("locus_name", sort=False):
            row = {"locus_name": locus_name}
            for _, r in grp.iterrows():
                key = r["contrast"]
                row[key + "_ratio"] = r["elevation_ratio"]
                row[key + "_nsig"] = r["n_bins_sig_fdr05"]
                row[key + "_topkl"] = r["top_track_contributor"]
            pivot_rows.append(row)
        pivot_df = pd.DataFrame(pivot_rows)
        pivot_df.to_csv(args.pivot_out, sep="\t", index=False,
                         float_format="%.3f")
        print("Wrote pivot to", args.pivot_out)


if __name__ == "__main__":
    main()
