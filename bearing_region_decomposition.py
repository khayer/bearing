#!/usr/bin/env python3
"""
bearing_region_decomposition.py

Region-restricted BEARING decomposition report.

Given a directory of `diff_<A>_vs_<B>.stats.tsv` files (BEARING differential
output) and a BED file of named genomic regions, this script summarizes per
region:

  1. Per-region track contribution fraction (which tracks drive the
     compositional shift within each region)
  2. Top-K bins by |bearing_score| in each region, with full kl decomposition
  3. Sum and mean signed/absolute kl_i per region

This complements the genome-wide bearing_decomposition_report.py. Genome-wide
aggregation is dominated by bin-count-heavy comparisons (e.g., lymphoid vs
fibroblast) and by sparse, spatially-concentrated tracks (RNA-seq). For focal
loci such as antigen-receptor loci, region-restricted decomposition reveals
architectural (CTCF/Cohesin/H3K27ac) signal that is invisible in the
genome-wide aggregate.

Usage:
  python bearing_region_decomposition.py \\
      --diff-dir results_v6 \\
      --regions tcrb_regions_v5.bed \\
      --categories DN_rep1_cats.json \\
      --out tcrb_decomposition \\
      [--top-k 10] [--abs-bes-threshold 0.5]

Inputs:
  --diff-dir       Directory with diff_<A>_vs_<B>.stats.tsv files
  --regions        BED file with at least 3 columns (chrom, start, end);
                   optional 4th column = region name (auto-generated if absent)
  --categories     A *_cats.json file defining track names and colors
  --out            Output prefix
  --top-k          Number of top-|BES| bins to show per region (default 10)
  --abs-bes-threshold  Minimum |BES| for a bin to count in region fraction
                       statistics (default 0.5; helps exclude noise)

Outputs:
  <out>.pdf                One page per region with all comparisons
  <out>.per_region.tsv     Per-(region, comparison, track) contribution stats
  <out>.top_bins.tsv       Top-K bins per (region, comparison) with kl values
"""

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


def parse_categories(cats_json_path):
    """Read a *_cats.json file and return ordered list of (name, color)."""
    with open(cats_json_path) as f:
        d = json.load(f)
    cats = d["categories"]
    ordered = []
    for k in sorted(cats.keys(), key=int):
        name, color = cats[k]
        ordered.append((name, color))
    return ordered


def parse_comparison_name(filename):
    base = Path(filename).name
    base = re.sub(r"\.stats\.tsv(\.gz)?$", "", base)
    base = re.sub(r"^diff_", "", base)
    m = re.match(r"^(.+?)_vs_(.+?)$", base)
    if not m:
        raise ValueError(f"Cannot parse comparison name from: {filename}")
    return m.group(1), m.group(2)


def load_regions_bed(bed_path):
    """Read a BED file. Returns list of (chrom, start, end, name)."""
    regions = []
    with open(bed_path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"region_{i}"
            regions.append((chrom, start, end, name))
    return regions


def load_diff_tsv(path, n_tracks):
    df = pd.read_csv(path, sep="\t")
    kl_cols = [c for c in df.columns if c.startswith("kl_")]
    kl_cols.sort(key=lambda c: int(c.split("_")[1]))
    if len(kl_cols) != n_tracks:
        sys.exit(
            f"ERROR: {path}: expected {n_tracks} kl_* columns, got "
            f"{len(kl_cols)}: {kl_cols}"
        )
    return df, kl_cols


def subset_to_region(df, chrom, start, end):
    """Return bins overlapping the region."""
    mask = (df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)
    return df[mask].copy()


def compute_region_metrics(df_region, kl_cols, track_names,
                            abs_bes_threshold=0.5):
    """For a region's bins, compute summary metrics. Returns dict keyed by
    direction ('A', 'B', 'all').
    """
    out = {}
    if len(df_region) == 0:
        for d in ("A", "B", "all"):
            out[d] = {
                "n_bins": 0,
                "n_bins_above_threshold": 0,
                "frac_of_abs_bes": dict(zip(track_names,
                                              [0.0] * len(track_names))),
                "mean_kl": dict(zip(track_names,
                                      [0.0] * len(track_names))),
                "sum_abs_kl": dict(zip(track_names,
                                         [0.0] * len(track_names))),
            }
        return out

    # Filter to bins above threshold for contribution fraction
    thresh_mask = df_region["bearing_score_tested"] >= abs_bes_threshold

    for dir_label, dir_filter in [("A", df_region["direction"] == "+"),
                                    ("B", df_region["direction"] == "-"),
                                    ("all", pd.Series([True] *
                                                       len(df_region),
                                                       index=df_region.index))]:
        sub = df_region[dir_filter & thresh_mask]
        d = {
            "n_bins": int(dir_filter.sum()),
            "n_bins_above_threshold": len(sub),
        }
        if len(sub) > 0:
            kl_mat = sub[kl_cols].to_numpy()
            abs_kl = np.abs(kl_mat)
            total_abs = abs_kl.sum()
            if total_abs > 0:
                track_frac = abs_kl.sum(axis=0) / total_abs
            else:
                track_frac = np.zeros(len(kl_cols))
            d["frac_of_abs_bes"] = dict(zip(track_names,
                                              track_frac.tolist()))
            d["mean_kl"] = dict(zip(track_names,
                                      kl_mat.mean(axis=0).tolist()))
            d["sum_abs_kl"] = dict(zip(track_names,
                                         abs_kl.sum(axis=0).tolist()))
        else:
            d["frac_of_abs_bes"] = dict(zip(track_names,
                                              [0.0] * len(track_names)))
            d["mean_kl"] = dict(zip(track_names,
                                      [0.0] * len(track_names)))
            d["sum_abs_kl"] = dict(zip(track_names,
                                         [0.0] * len(track_names)))
        out[dir_label] = d
    return out


def get_top_bins(df_region, kl_cols, track_names, top_k=10):
    """Return top-K bins in region by |bearing_score| with full kl values."""
    if len(df_region) == 0:
        return pd.DataFrame()
    sorted_df = df_region.sort_values("bearing_score_tested",
                                       ascending=False).head(top_k)
    cols_out = ["chrom", "start", "end", "bearing_score",
                "significant_fdr0.05", "direction"]
    cols_out = [c for c in cols_out if c in sorted_df.columns]
    cols_out += kl_cols
    return sorted_df[cols_out].rename(
        columns={kl_cols[i]: track_names[i]
                 for i in range(len(track_names))}
    )


# ---- Plotting -----------------------------------------------------------

def plot_region_track_fraction(metrics_per_comparison, comparisons,
                                 track_names, track_colors, ax,
                                 region_name, region_label,
                                 direction_label="all"):
    """Stacked bar: one bar per comparison, segments = track contribution
    fraction (frac_of_abs_bes), for bins in this region above threshold."""
    n = len(comparisons)
    x = np.arange(n)
    bottoms = np.zeros(n)

    for j, t in enumerate(track_names):
        heights = np.zeros(n)
        for i, (a, b) in enumerate(comparisons):
            m = metrics_per_comparison[(a, b)]
            heights[i] = m[direction_label]["frac_of_abs_bes"].get(t, 0)
        ax.bar(x, heights, bottom=bottoms, color=track_colors[j],
               edgecolor="white", linewidth=0.5, label=t)
        bottoms += heights

    labels = []
    for a, b in comparisons:
        n_bins = metrics_per_comparison[(a, b)][direction_label][
            "n_bins_above_threshold"]
        labels.append(f"{a} vs {b}\n(n={n_bins})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Fraction of total |BES|")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0 / len(track_names), color="gray", linestyle="--",
                linewidth=0.6, alpha=0.5,
                label=f"equal share ({1.0/len(track_names):.2f})")
    ax.set_title(f"{region_name}: track contribution to compositional shift "
                  f"(bins above |BES| threshold)")


def plot_region_top_bins_heatmap(top_bins_per_comp, comparisons,
                                   track_names, track_colors, ax,
                                   region_name):
    """Heatmap: rows = bins (concatenated across comparisons), cols = tracks,
    values = signed kl. Below each row, comparison and |BES|.
    """
    rows = []
    row_labels = []
    for (a, b) in comparisons:
        df_top = top_bins_per_comp.get((a, b))
        if df_top is None or len(df_top) == 0:
            continue
        for _, r in df_top.iterrows():
            vals = [r[t] for t in track_names]
            rows.append(vals)
            sig = "*" if r.get("significant_fdr0.05", 0) == 1 else ""
            row_labels.append(
                f"{a}/{b}  {r['chrom']}:{int(r['start']):,}  "
                f"BES={r['bearing_score']:+.2f}{sig}"
            )
    if not rows:
        ax.text(0.5, 0.5, f"No top bins in {region_name}",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, alpha=0.6)
        ax.axis("off")
        return

    data = np.array(rows)
    vmax = max(0.5, np.nanmax(np.abs(data)))
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(track_names)))
    ax.set_xticklabels(track_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=6.5)
    ax.set_title(f"{region_name}: top bins by |BES| per comparison "
                  f"(* = FDR-significant)")
    # Annotate non-zero cells
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if abs(v) > 0.1:
                color = "black" if abs(v) < vmax * 0.6 else "white"
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=6, color=color)
    plt.colorbar(im, ax=ax, label="Signed kl_i")


# ---- Main ---------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--diff-dir", required=True)
    ap.add_argument("--regions", required=True,
                    help="BED file with regions (cols: chrom start end [name])")
    ap.add_argument("--categories", required=True,
                    help="A *_cats.json file (any sample's)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-k", type=int, default=10,
                    help="Top-K bins per region by |BES| (default 10)")
    ap.add_argument("--abs-bes-threshold", type=float, default=0.5,
                    help="Minimum |BES| for region fraction stats "
                         "(default 0.5)")
    ap.add_argument("--filename-pattern", default="diff_*.stats.tsv")
    args = ap.parse_args()

    # Load category definitions
    cats = parse_categories(args.categories)
    track_names = [c[0] for c in cats]
    track_colors = [c[1] for c in cats]
    n_tracks = len(track_names)
    print(f"Tracks ({n_tracks}): {track_names}", flush=True)

    # Load regions
    regions = load_regions_bed(args.regions)
    print(f"Loaded {len(regions)} regions from {args.regions}", flush=True)
    for r in regions:
        print(f"  {r[3]}: {r[0]}:{r[1]:,}-{r[2]:,}", flush=True)

    # Enumerate comparisons
    diff_files = sorted(Path(args.diff_dir).glob(args.filename_pattern))
    if not diff_files:
        sys.exit(f"ERROR: no diff files in {args.diff_dir}")

    comparisons = []
    diff_dfs = {}
    for f in diff_files:
        try:
            a, b = parse_comparison_name(str(f))
        except ValueError as ex:
            print(f"  WARN: skipping {f}: {ex}", file=sys.stderr)
            continue
        print(f"  Loading {a} vs {b}...", flush=True)
        df, kl_cols = load_diff_tsv(str(f), n_tracks)
        diff_dfs[(a, b)] = (df, kl_cols)
        comparisons.append((a, b))

    print(f"\nProcessing {len(regions)} regions x "
          f"{len(comparisons)} comparisons", flush=True)

    # For each (region, comparison): compute metrics + top-K bins
    all_metrics = {}     # (region_name, A, B) -> metrics dict
    all_top_bins = {}    # (region_name, A, B) -> DataFrame
    rows_tsv = []
    top_rows = []

    for chrom, start, end, region_name in regions:
        print(f"\n  Region {region_name}: {chrom}:{start:,}-{end:,}",
              flush=True)
        for (a, b) in comparisons:
            df, kl_cols = diff_dfs[(a, b)]
            df_region = subset_to_region(df, chrom, start, end)
            metrics = compute_region_metrics(df_region, kl_cols,
                                              track_names,
                                              abs_bes_threshold=
                                              args.abs_bes_threshold)
            all_metrics[(region_name, a, b)] = metrics
            top_bins = get_top_bins(df_region, kl_cols, track_names,
                                     top_k=args.top_k)
            all_top_bins[(region_name, a, b)] = top_bins

            # TSV summary
            for dir_label in ("A", "B", "all"):
                d = metrics[dir_label]
                for t in track_names:
                    rows_tsv.append({
                        "region": region_name,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "comparison_A": a,
                        "comparison_B": b,
                        "direction": dir_label,
                        "n_bins_in_direction": d["n_bins"],
                        "n_bins_above_threshold": d["n_bins_above_threshold"],
                        "track": t,
                        "frac_of_abs_bes": d["frac_of_abs_bes"][t],
                        "mean_kl": d["mean_kl"][t],
                        "sum_abs_kl": d["sum_abs_kl"][t],
                    })

            # Top-bins TSV
            if len(top_bins):
                for _, r in top_bins.iterrows():
                    row_out = {"region": region_name,
                                "comparison_A": a, "comparison_B": b,
                                "chrom": r["chrom"], "start": r["start"],
                                "end": r["end"],
                                "bearing_score": r["bearing_score"],
                                "FDR_sig": r.get("significant_fdr0.05", 0),
                                "direction": r.get("direction", "")}
                    for t in track_names:
                        row_out[f"kl_{t}"] = r[t]
                    top_rows.append(row_out)

    # Write TSVs
    pd.DataFrame(rows_tsv).to_csv(f"{args.out}.per_region.tsv",
                                    sep="\t", index=False)
    print(f"\nWrote {args.out}.per_region.tsv", flush=True)

    pd.DataFrame(top_rows).to_csv(f"{args.out}.top_bins.tsv",
                                    sep="\t", index=False)
    print(f"Wrote {args.out}.top_bins.tsv", flush=True)

    # Build PDF: summary page first, then one page per region
    pdf_path = f"{args.out}.pdf"
    with PdfPages(pdf_path) as pdf:

        # SUMMARY PAGE: track fraction heatmap across all regions x comparisons
        # Rows = (region, comparison), cols = tracks; one heatmap.
        n_rows = len(regions) * len(comparisons)
        summary_data = np.zeros((n_rows, len(track_names)))
        row_labels = []
        idx = 0
        for chrom, start, end, region_name in regions:
            for (a, b) in comparisons:
                m = all_metrics[(region_name, a, b)]["all"]
                for j, t in enumerate(track_names):
                    summary_data[idx, j] = m["frac_of_abs_bes"].get(t, 0)
                n = m["n_bins_above_threshold"]
                row_labels.append(f"{region_name}  |  {a} vs {b}  (n={n})")
                idx += 1

        fig, ax = plt.subplots(figsize=(11, max(4, 0.3 * n_rows)))
        vmax = max(0.5, np.nanmax(summary_data) if summary_data.size else 0.5)
        im = ax.imshow(summary_data, aspect="auto", cmap="YlOrRd",
                       vmin=0, vmax=vmax)
        ax.set_xticks(range(len(track_names)))
        ax.set_xticklabels(track_names, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(row_labels, fontsize=7)
        ax.set_title(f"Per-region track contribution to compositional shift\n"
                      f"(bins with |BES| >= {args.abs_bes_threshold}; "
                      f"each row sums to ~1.0)",
                      fontsize=11)
        equal_share = 1.0 / len(track_names)
        for i in range(n_rows):
            for j in range(len(track_names)):
                v = summary_data[i, j]
                if v > 0.02:
                    color = "black" if v < vmax * 0.6 else "white"
                    marker = "*" if v > equal_share * 1.5 else ""
                    ax.text(j, i, f"{v:.2f}{marker}", ha="center",
                            va="center", fontsize=6, color=color)
        plt.colorbar(im, ax=ax,
                      label=f"Fraction of |BES|  "
                            f"(equal share={equal_share:.2f}; * = >1.5× equal)")
        # Add horizontal lines separating regions
        for i in range(1, len(regions)):
            ax.axhline(i * len(comparisons) - 0.5, color="black",
                        linewidth=0.6, alpha=0.4)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Per-region detail pages
        for chrom, start, end, region_name in regions:
            # Build metrics dict for THIS region only
            region_metrics = {(a, b): all_metrics[(region_name, a, b)]
                              for (a, b) in comparisons}
            region_top = {(a, b): all_top_bins[(region_name, a, b)]
                          for (a, b) in comparisons}

            # Compute layout: top half = fraction bar, bottom half = heatmap
            # Heatmap height depends on number of rows (top_k * n_comparisons)
            n_heatmap_rows = sum(len(region_top[(a, b)])
                                  for (a, b) in comparisons)
            heatmap_height = max(4, 0.18 * n_heatmap_rows)

            fig = plt.figure(figsize=(13, 6 + heatmap_height))
            gs = fig.add_gridspec(2, 1, height_ratios=[5, heatmap_height])

            ax_top = fig.add_subplot(gs[0])
            plot_region_track_fraction(region_metrics, comparisons,
                                          track_names, track_colors,
                                          ax_top, region_name,
                                          f"{chrom}:{start:,}-{end:,}",
                                          direction_label="all")
            ax_top.legend(loc="upper right", fontsize=7, ncol=2,
                           framealpha=0.9)

            ax_bot = fig.add_subplot(gs[1])
            plot_region_top_bins_heatmap(region_top, comparisons,
                                            track_names, track_colors,
                                            ax_bot, region_name)

            fig.suptitle(f"{region_name}  ({chrom}:{start:,}-{end:,})",
                          fontsize=12, y=0.995)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"Wrote {pdf_path}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
