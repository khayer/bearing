#!/usr/bin/env python3
"""
bearing_specificity_profile.py

Compute and plot a sliding-window 1D "specificity profile" of BEARING
compositional shift magnitude across a chromosomal region.

This is the quantitative companion to the wide-locus claim that BEARING
signal is concentrated at antigen-receptor loci and falls off in flanking
regions. It is a positive control for the BEARING methodology: if the
method were producing global noise, |BES| would be similar throughout the
window; instead, |BES| is high at biologically meaningful loci (Tcrb,
Tcra) and near baseline elsewhere.

Three optional metrics, each toggleable:

  1. Mean |BES| per window  -  compositional shift magnitude
  2. FDR-significant bin fraction per window  -  power-adjusted signal
  3. Per-track fraction of |BES|  -  which assays dominate in each window

Usage:
  python bearing_specificity_profile.py \\
      --diff results_v6/diff_DN_vs_EbKO.stats.tsv \\
      --categories DN_rep1_cats.json \\
      --region chr6:40000000-47400000 \\
      --window-size 100000 --step 25000 \\
      --annotations tcrb_regions_v5.bed \\
      --out tcrb_wide_specificity_DN_vs_EbKO

  # Or for all diff files in a directory:
  python bearing_specificity_profile.py \\
      --diff-dir results_v6 \\
      --categories DN_rep1_cats.json \\
      --region chr6:40000000-47400000 \\
      --annotations tcrb_regions_v5.bed \\
      --out tcrb_wide_specificity_all

Inputs:
  --diff FILE         A single diff_*_vs_*.stats.tsv file
  --diff-dir DIR      Directory containing multiple diff_*.stats.tsv files
                      (mutually exclusive with --diff)
  --categories JSON   A *_cats.json file (any sample's)
  --region CHR:S-E    Genomic window to profile
  --window-size BP    Sliding window size (default 100,000)
  --step BP           Step size for sliding (default 25,000)
  --annotations BED   Optional BED file of feature regions to overlay
  --out PREFIX        Output file prefix (creates .pdf, .tsv, .bedgraph)

Outputs:
  <out>.pdf           Multi-panel specificity profile plot
  <out>.tsv           Per-window statistics (machine-readable)
  <out>.bedgraph      Mean |BES| per window in bedGraph format (for IGV)
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


def parse_categories(path):
    with open(path) as f:
        d = json.load(f)
    cats = d["categories"]
    ordered = []
    for k in sorted(cats.keys(), key=int):
        name, color = cats[k]
        ordered.append((name, color))
    return ordered


def parse_region_str(s):
    m = re.match(r"^([^:]+):([0-9,]+)-([0-9,]+)$", s)
    if not m:
        raise argparse.ArgumentTypeError(f"Cannot parse region: {s}")
    return m.group(1), int(m.group(2).replace(",", "")), \
        int(m.group(3).replace(",", ""))


def parse_comparison_name(filename):
    base = Path(filename).name
    base = re.sub(r"\.stats\.tsv(\.gz)?$", "", base)
    base = re.sub(r"^diff_", "", base)
    m = re.match(r"^(.+?)_vs_(.+?)$", base)
    if not m:
        raise ValueError(f"Cannot parse: {filename}")
    return m.group(1), m.group(2)


def load_diff_tsv(path, n_tracks):
    df = pd.read_csv(path, sep="\t")
    kl_cols = sorted([c for c in df.columns if c.startswith("kl_")],
                     key=lambda c: int(c.split("_")[1]))
    if len(kl_cols) != n_tracks:
        sys.exit(f"ERROR: {path}: expected {n_tracks} kl_* cols, "
                 f"got {len(kl_cols)}")
    return df, kl_cols


def load_annotations(bed_path):
    if not bed_path:
        return []
    regions = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"region_{len(regions)+1}"
            regions.append((chrom, start, end, name))
    return regions


def compute_window_stats(df, kl_cols, track_names, chrom, region_start,
                          region_end, window_size, step,
                          fdr_col="significant_fdr0.05"):
    """For each sliding window, compute mean |BES|, FDR fraction, per-track
    contribution fraction.
    """
    # Subset to region
    mask = ((df["chrom"] == chrom) &
            (df["end"] > region_start) &
            (df["start"] < region_end))
    region_df = df[mask].copy()
    n_region_bins = len(region_df)

    # Define windows
    win_starts = np.arange(region_start, region_end, step)
    rows = []

    for ws in win_starts:
        we = ws + window_size
        bin_mask = ((region_df["end"] > ws) &
                    (region_df["start"] < we))
        sub = region_df[bin_mask]
        n_bins = len(sub)

        if n_bins == 0:
            row = {
                "chrom": chrom,
                "win_start": ws,
                "win_end": we,
                "win_mid": ws + window_size // 2,
                "n_bins": 0,
                "mean_abs_bes": 0.0,
                "max_abs_bes": 0.0,
                "n_fdr_sig": 0,
                "frac_fdr_sig": 0.0,
            }
            for t in track_names:
                row[f"frac_{t}"] = 0.0
                row[f"mean_abs_kl_{t}"] = 0.0
            rows.append(row)
            continue

        abs_bes = sub["bearing_score_tested"].to_numpy()
        kl_mat = sub[kl_cols].to_numpy()
        abs_kl = np.abs(kl_mat)
        total_abs = abs_kl.sum()

        if fdr_col in sub.columns:
            n_fdr = int(sub[fdr_col].astype(int).sum())
        else:
            n_fdr = 0

        row = {
            "chrom": chrom,
            "win_start": ws,
            "win_end": we,
            "win_mid": ws + window_size // 2,
            "n_bins": n_bins,
            "mean_abs_bes": float(abs_bes.mean()),
            "max_abs_bes": float(abs_bes.max()),
            "n_fdr_sig": n_fdr,
            "frac_fdr_sig": n_fdr / n_bins,
        }
        if total_abs > 0:
            track_frac = abs_kl.sum(axis=0) / total_abs
        else:
            track_frac = np.zeros(len(kl_cols))
        for j, t in enumerate(track_names):
            row[f"frac_{t}"] = float(track_frac[j])
            row[f"mean_abs_kl_{t}"] = float(abs_kl[:, j].mean())
        rows.append(row)

    return pd.DataFrame(rows), n_region_bins


def plot_summary_heatmap(all_stats_dfs, annotations, region_start,
                          region_end, chrom, window_size, step):
    """Build a heatmap: rows = comparisons, columns = windows, color = mean |BES|.
    This is the summary "are signals specific to AR loci across all comparisons?" view.
    """
    comparisons = sorted(all_stats_dfs.keys())
    if not comparisons:
        return None
    first_df = all_stats_dfs[comparisons[0]]
    win_mid = first_df["win_mid"].to_numpy()
    n_win = len(win_mid)
    n_comp = len(comparisons)

    # Build matrix
    data = np.zeros((n_comp, n_win))
    for i, comp in enumerate(comparisons):
        data[i] = all_stats_dfs[comp]["mean_abs_bes"].to_numpy()

    # Figure: heatmap with annotation strip above
    fig = plt.figure(figsize=(14, max(4, 0.4 * n_comp + 1.5)))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.5, n_comp * 0.5])

    # Annotation panel
    ax_anno = fig.add_subplot(gs[0])
    ax_anno.set_ylim(0, 1.5)
    for c2, s2, e2, n2 in annotations:
        if c2 != chrom:
            continue
        if e2 < region_start or s2 > region_end:
            continue
        width = e2 - s2
        width_frac = width / (region_end - region_start)
        ax_anno.add_patch(plt.Rectangle((s2, 0.2), width, 0.6,
                                          facecolor="#d0e8ff",
                                          edgecolor="navy", linewidth=1.0))
        if width_frac < 0.04:
            ax_anno.annotate(n2, xy=((s2 + e2) / 2, 0.8),
                              xytext=((s2 + e2) / 2, 1.1),
                              ha="center", va="bottom", fontsize=7,
                              fontweight="bold",
                              arrowprops=dict(arrowstyle="-", color="navy",
                                               linewidth=0.5))
        else:
            ax_anno.text((s2 + e2) / 2, 0.5, n2,
                          ha="center", va="center",
                          fontsize=8, fontweight="bold")
    ax_anno.set_yticks([])
    ax_anno.set_ylabel("Features", fontsize=8, rotation=0,
                        ha="right", va="center")
    ax_anno.spines["top"].set_visible(False)
    ax_anno.spines["right"].set_visible(False)
    ax_anno.spines["left"].set_visible(False)
    ax_anno.set_xlim(region_start, region_end)
    ax_anno.tick_params(axis="x", labelbottom=False)

    # Heatmap
    ax = fig.add_subplot(gs[1], sharex=ax_anno)
    vmax = max(data.max(), 0.5) if data.size else 1.0
    im = ax.imshow(
        data, aspect="auto", cmap="magma_r",
        extent=[region_start, region_end, n_comp, 0],
        interpolation="nearest", vmin=0, vmax=vmax,
    )
    ax.set_yticks(np.arange(n_comp) + 0.5)
    ax.set_yticklabels([f"{a} vs {b}" for (a, b) in comparisons],
                        fontsize=9)
    ax.set_xlim(region_start, region_end)
    ax.set_xlabel(f"{chrom} position (bp)", fontsize=10)
    ax.ticklabel_format(useOffset=False, style="plain", axis="x")
    plt.colorbar(im, ax=ax, label="Mean |BES| per window",
                  fraction=0.02, pad=0.01)

    fig.suptitle(
        f"BEARING wide-locus specificity profile across comparisons\n"
        f"{chrom}:{region_start:,}-{region_end:,}  "
        f"(window={window_size:,} bp, step={step:,} bp)",
        fontsize=11, y=0.995
    )
    plt.tight_layout()
    return fig


def plot_specificity_profile(stats_df, annotations, track_names,
                              track_colors, region_start, region_end,
                              chrom, comparison_label, panels_to_show):
    """Build the multi-panel figure. panels_to_show is a list:
    ['mean_bes', 'fdr', 'per_track', 'annotations'] in display order.
    """
    n_panels = len(panels_to_show)
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(13, 1.6 * n_panels + 1),
        sharex=True,
        gridspec_kw={"hspace": 0.2,
                     "height_ratios": [
                         0.4 if p == "annotations" else 1.0
                         for p in panels_to_show
                     ]}
    )
    if n_panels == 1:
        axes = [axes]

    win_mid = stats_df["win_mid"].to_numpy()

    for ax_idx, panel in enumerate(panels_to_show):
        ax = axes[ax_idx]

        if panel == "annotations":
            ax.set_ylim(0, 1)
            for c2, s2, e2, n2 in annotations:
                if c2 != chrom:
                    continue
                if e2 < region_start or s2 > region_end:
                    continue
                width = e2 - s2
                width_frac = width / (region_end - region_start)
                ax.add_patch(plt.Rectangle((s2, 0.2), width, 0.6,
                                            facecolor="#d0e8ff",
                                            edgecolor="navy",
                                            linewidth=1.0))
                if width_frac < 0.04:
                    ax.annotate(n2, xy=((s2 + e2) / 2, 0.8),
                                 xytext=((s2 + e2) / 2, 1.05),
                                 ha="center", va="bottom", fontsize=7,
                                 fontweight="bold",
                                 arrowprops=dict(arrowstyle="-",
                                                  color="navy",
                                                  linewidth=0.5))
                else:
                    ax.text((s2 + e2) / 2, 0.5, n2,
                             ha="center", va="center",
                             fontsize=8, fontweight="bold")
            ax.set_yticks([])
            ax.set_ylabel("Features", fontsize=8, rotation=0,
                           ha="right", va="center")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.set_ylim(0, 1.5)

        elif panel == "mean_bes":
            y = stats_df["mean_abs_bes"].to_numpy()
            ax.fill_between(win_mid, 0, y, color="#404040", alpha=0.7,
                             linewidth=0)
            ax.plot(win_mid, y, color="black", linewidth=0.8)
            ax.set_ylabel("Mean |BES|\nper window",
                           fontsize=9, rotation=0, ha="right", va="center")
            # Annotate the top windows
            top_idx = np.argsort(y)[-3:][::-1]
            for ti in top_idx:
                if y[ti] > 0.1:
                    ax.annotate(f"{y[ti]:.2f}",
                                 xy=(win_mid[ti], y[ti]),
                                 xytext=(0, 5),
                                 textcoords="offset points",
                                 fontsize=7, ha="center",
                                 color="darkred")
            # Baseline marker - the median of the lower 50% of windows
            baseline = np.median(np.sort(y)[:len(y) // 2])
            ax.axhline(baseline, color="gray", linestyle="--",
                        linewidth=0.5, alpha=0.7)
            ax.text(region_end, baseline,
                     f" baseline\n {baseline:.3f}",
                     fontsize=6, color="gray", va="center")
            ax.set_ylim(0, max(y.max() * 1.1, 0.01))

        elif panel == "fdr":
            y = stats_df["frac_fdr_sig"].to_numpy()
            if y.max() == 0:
                ax.text(0.5, 0.5,
                         "No FDR-significant bins in this window\n"
                         "(typical for N=10 perms on lymphoid-vs-lymphoid)",
                         ha="center", va="center",
                         transform=ax.transAxes, fontsize=9, color="gray")
                ax.set_yticks([])
            else:
                ax.fill_between(win_mid, 0, y, color="#cc6600", alpha=0.7,
                                 linewidth=0)
                ax.plot(win_mid, y, color="#aa4400", linewidth=0.8)
            ax.set_ylabel("Frac FDR-sig\nper window",
                           fontsize=9, rotation=0, ha="right", va="center")

        elif panel == "per_track":
            # Stacked area plot of per-track fraction of |BES|
            stack_data = []
            for t in track_names:
                stack_data.append(stats_df[f"frac_{t}"].to_numpy())
            stack_data = np.array(stack_data)
            # Use ax.stackplot
            ax.stackplot(win_mid, stack_data,
                          colors=track_colors, edgecolor="white",
                          linewidth=0.1, alpha=0.9)
            ax.set_ylabel("Track\ncontribution",
                           fontsize=9, rotation=0, ha="right", va="center")
            ax.set_ylim(0, 1)
            # Legend
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=c, label=n)
                        for c, n in zip(track_colors, track_names)]
            ax.legend(handles=handles, loc="upper left", fontsize=6,
                       ncol=len(track_names), framealpha=0.9,
                       bbox_to_anchor=(0.02, 1.18))

    axes[-1].set_xlim(region_start, region_end)
    axes[-1].set_xlabel(f"{chrom} position (bp)", fontsize=10)
    axes[-1].ticklabel_format(useOffset=False, style="plain", axis="x")

    fig.suptitle(f"BEARING wide-locus specificity profile: {comparison_label}\n"
                  f"{chrom}:{region_start:,}-{region_end:,}  "
                  f"(window={stats_df.iloc[0]['win_end'] - stats_df.iloc[0]['win_start']:,} bp, "
                  f"step={stats_df.iloc[1]['win_start'] - stats_df.iloc[0]['win_start']:,} bp)",
                  fontsize=11, y=0.995)
    plt.tight_layout()
    return fig


def write_bedgraph(stats_df, chrom, out_path, value_col="mean_abs_bes"):
    """Write a bedGraph file of the chosen value column."""
    with open(out_path, "w") as f:
        f.write(f"track type=bedGraph name=\"{value_col}\" "
                f"description=\"BEARING specificity profile\" visibility=full\n")
        for _, r in stats_df.iterrows():
            f.write(f"{chrom}\t{int(r['win_start'])}\t{int(r['win_end'])}\t"
                    f"{r[value_col]:.6f}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--diff", help="Single diff_*.stats.tsv file")
    grp.add_argument("--diff-dir",
                     help="Directory containing diff_*.stats.tsv files")
    ap.add_argument("--filename-pattern", default="diff_*.stats.tsv",
                    help="Glob pattern when using --diff-dir")
    ap.add_argument("--categories", required=True,
                    help="A *_cats.json file")
    ap.add_argument("--region", required=True, type=parse_region_str,
                    help="Region in format chrom:start-end")
    ap.add_argument("--window-size", type=int, default=50000,
                    help="Sliding window size (default 50,000; "
                         "matches Hi-C resolution)")
    ap.add_argument("--step", type=int, default=50000,
                    help="Step size for sliding (default 50,000; "
                         "non-overlapping at Hi-C resolution)")
    ap.add_argument("--annotations", default=None,
                    help="BED file with feature regions for annotation panel")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("--panels", nargs="+",
                    default=["annotations", "mean_bes", "fdr", "per_track"],
                    choices=["annotations", "mean_bes", "fdr", "per_track"],
                    help="Which panels to show, in order")
    args = ap.parse_args()

    cats = parse_categories(args.categories)
    track_names = [c[0] for c in cats]
    track_colors = [c[1] for c in cats]
    n_tracks = len(track_names)
    print(f"Tracks ({n_tracks}): {track_names}", flush=True)

    chrom, region_start, region_end = args.region
    print(f"Region: {chrom}:{region_start:,}-{region_end:,}", flush=True)
    print(f"Window: {args.window_size:,} bp, step: {args.step:,} bp",
          flush=True)
    n_windows = (region_end - region_start) // args.step
    print(f"Approximate window count: {n_windows}", flush=True)

    annotations = load_annotations(args.annotations)
    print(f"Loaded {len(annotations)} feature regions", flush=True)

    # Enumerate diff files
    if args.diff:
        diff_files = [Path(args.diff)]
    else:
        diff_files = sorted(Path(args.diff_dir).glob(args.filename_pattern))
        if not diff_files:
            sys.exit(f"ERROR: no diff files in {args.diff_dir}")

    print(f"\nProcessing {len(diff_files)} diff file(s)", flush=True)

    pdf_path = f"{args.out}.pdf"
    all_stats_rows = []
    all_stats_dfs = {}  # (A, B) -> stats_df for summary heatmap

    with PdfPages(pdf_path) as pdf:
        # First pass: compute stats for all comparisons (so we can build
        # the summary heatmap before the detail pages)
        for diff_file in diff_files:
            try:
                a, b = parse_comparison_name(str(diff_file))
            except ValueError as ex:
                print(f"  WARN: skipping {diff_file}: {ex}",
                      file=sys.stderr)
                continue
            comparison_label = f"{a} vs {b}"
            print(f"  {comparison_label} ...", flush=True)

            df, kl_cols = load_diff_tsv(str(diff_file), n_tracks)
            stats_df, n_region_bins = compute_window_stats(
                df, kl_cols, track_names, chrom, region_start, region_end,
                args.window_size, args.step
            )
            print(f"    {n_region_bins} bins in region, "
                  f"{len(stats_df)} windows", flush=True)
            if stats_df["mean_abs_bes"].max() == 0:
                print(f"    WARN: no signal in this region for "
                      f"{comparison_label}", file=sys.stderr)

            stats_df["comparison_A"] = a
            stats_df["comparison_B"] = b
            all_stats_rows.append(stats_df)
            all_stats_dfs[(a, b)] = stats_df

            # Per-comparison bedGraph
            if len(diff_files) > 1:
                bg_path = f"{args.out}_{a}_vs_{b}.bedgraph"
            else:
                bg_path = f"{args.out}.bedgraph"
            write_bedgraph(stats_df, chrom, bg_path)
            print(f"    wrote {bg_path}", flush=True)

        # Page 1: summary heatmap (if multiple comparisons)
        if len(all_stats_dfs) > 1:
            print("\nBuilding summary heatmap (page 1)...", flush=True)
            fig = plot_summary_heatmap(
                all_stats_dfs, annotations, region_start, region_end,
                chrom, args.window_size, args.step
            )
            if fig is not None:
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        # Subsequent pages: per-comparison detail
        for (a, b), stats_df in all_stats_dfs.items():
            comparison_label = f"{a} vs {b}"
            fig = plot_specificity_profile(
                stats_df, annotations, track_names, track_colors,
                region_start, region_end, chrom,
                comparison_label, args.panels
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    # Combined TSV
    combined = pd.concat(all_stats_rows, ignore_index=True)
    tsv_path = f"{args.out}.tsv"
    combined.to_csv(tsv_path, sep="\t", index=False)
    print(f"\nWrote {pdf_path}", flush=True)
    print(f"Wrote {tsv_path}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
