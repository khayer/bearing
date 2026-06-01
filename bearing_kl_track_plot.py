#!/usr/bin/env python3
"""
bearing_kl_track_plot.py

Plot BEARING per-track KL contributions (kl_1...kl_K) as separate genome-browser-
style tracks along a genomic interval, with optional region annotations.

Use this to visualize WHICH tracks contribute to the BEARING compositional shift
at a focal locus. Each track is plotted as a signed bar/line: positive = condition
A favored, negative = condition B favored, height = |kl_i|.

Designed for the Tcrb/Igh focal-locus figures in the BEARING paper.

Usage:
  python bearing_kl_track_plot.py \\
      --diff results_v6/diff_DN_vs_EbKO.stats.tsv \\
      --categories DN_rep1_cats.json \\
      --region chr6:40790000-41690000 \\
      --regions-bed tcrb_regions_v5.bed \\
      --out tcrb_DN_vs_EbKO_kl_tracks.pdf

Optional:
  --tracks ATAC CTCF Cohesin H3K27ac        Only show these tracks
  --fdr-only                                 Only plot FDR-significant bins
  --label-a DN  --label-b EbKO              Override condition labels
  --highlight-threshold 1.5                  Highlight bins with |kl| > this
"""

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_categories(cats_json_path):
    with open(cats_json_path) as f:
        d = json.load(f)
    cats = d["categories"]
    ordered = []
    for k in sorted(cats.keys(), key=int):
        name, color = cats[k]
        ordered.append((name, color))
    return ordered


def parse_region_str(s):
    """Parse chrom:start-end string."""
    m = re.match(r"^([^:]+):([0-9,]+)-([0-9,]+)$", s)
    if not m:
        raise argparse.ArgumentTypeError(f"Cannot parse region: {s}")
    return m.group(1), int(m.group(2).replace(",", "")), \
        int(m.group(3).replace(",", ""))


def load_regions_bed(bed_path):
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
            name = parts[3] if len(parts) >= 4 else f"r{i}"
            regions.append((chrom, start, end, name))
    return regions


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--diff", required=True,
                    help="diff_<A>_vs_<B>.stats.tsv file")
    ap.add_argument("--categories", required=True,
                    help="A *_cats.json file defining track names/colors")
    ap.add_argument("--region", required=True, type=parse_region_str,
                    help="Region to plot, format chrom:start-end")
    ap.add_argument("--regions-bed", default=None,
                    help="Optional BED file of region annotations to overlay")
    ap.add_argument("--out", required=True, help="Output PDF or PNG path")
    ap.add_argument("--tracks", nargs="+", default=None,
                    help="Subset of tracks to plot (default: all)")
    ap.add_argument("--fdr-only", action="store_true",
                    help="Only plot FDR-significant bins")
    ap.add_argument("--label-a", default="A", help="Label for condition A")
    ap.add_argument("--label-b", default="B", help="Label for condition B")
    ap.add_argument("--highlight-threshold", type=float, default=1.5,
                    help="Highlight bins with |kl| above this (default 1.5)")
    ap.add_argument("--bin-width", type=int, default=200,
                    help="Bin width in bp (default 200; used for bar width)")
    args = ap.parse_args()

    cats = parse_categories(args.categories)
    all_track_names = [c[0] for c in cats]
    all_track_colors = [c[1] for c in cats]

    if args.tracks:
        plot_tracks = [t for t in all_track_names if t in args.tracks]
        if not plot_tracks:
            sys.exit(f"ERROR: none of --tracks {args.tracks} match known: "
                     f"{all_track_names}")
    else:
        plot_tracks = all_track_names
    track_indices = [all_track_names.index(t) for t in plot_tracks]
    track_colors = [all_track_colors[i] for i in track_indices]

    # Load diff
    print(f"Loading {args.diff}...", flush=True)
    df = pd.read_csv(args.diff, sep="\t")
    kl_cols_all = sorted([c for c in df.columns if c.startswith("kl_")],
                          key=lambda c: int(c.split("_")[1]))
    if len(kl_cols_all) != len(all_track_names):
        sys.exit(f"ERROR: expected {len(all_track_names)} kl_* cols, "
                 f"got {len(kl_cols_all)}")
    kl_cols_plot = [kl_cols_all[i] for i in track_indices]

    # Subset
    chrom, start, end = args.region
    mask = (df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)
    df_r = df[mask].copy().reset_index(drop=True)
    if args.fdr_only:
        df_r = df_r[df_r["significant_fdr0.05"].astype(int) == 1]
    print(f"  {len(df_r)} bins in region", flush=True)

    if len(df_r) == 0:
        sys.exit("No bins to plot.")

    # Load region annotations
    anno_regions = []
    if args.regions_bed:
        for c2, s2, e2, n2 in load_regions_bed(args.regions_bed):
            if c2 == chrom and e2 > start and s2 < end:
                anno_regions.append((c2, s2, e2, n2))

    # Figure: one panel per track + 1 for the total BES + 1 for annotation
    n_panels = len(plot_tracks) + 1 + (1 if anno_regions else 0)
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(13, max(5, 1.2 * n_panels)),
        sharex=True,
        gridspec_kw={"hspace": 0.15,
                     "height_ratios":
                     [0.4 if anno_regions else 0] + [0.6] +
                     [1.0] * len(plot_tracks)}
        if anno_regions else
        {"hspace": 0.15,
         "height_ratios": [0.6] + [1.0] * len(plot_tracks)}
    )
    if n_panels == 1:
        axes = [axes]

    panel_idx = 0

    # Panel: region annotations
    if anno_regions:
        ax = axes[panel_idx]
        panel_idx += 1
        ax.set_ylim(0, 1)
        # Stagger labels for overlapping regions
        for i, (_, s2, e2, n2) in enumerate(anno_regions):
            y_text = 0.5 if (i % 2 == 0) else 0.85
            y_box = 0.2
            ax.add_patch(plt.Rectangle((s2, y_box), e2 - s2, 0.6,
                                         facecolor="#d0e8ff",
                                         edgecolor="navy", linewidth=1.0))
            # If region is very narrow, put label above
            region_width_frac = (e2 - s2) / (end - start)
            if region_width_frac < 0.05:
                ax.annotate(n2, xy=((s2 + e2) / 2, 0.8),
                             xytext=((s2 + e2) / 2, 1.05),
                             ha="center", va="bottom", fontsize=7,
                             fontweight="bold",
                             arrowprops=dict(arrowstyle="-", color="navy",
                                              linewidth=0.5))
            else:
                ax.text((s2 + e2) / 2, 0.5, n2, ha="center", va="center",
                        fontsize=8, fontweight="bold")
        ax.set_yticks([])
        ax.set_ylabel("Regions", fontsize=8, rotation=0,
                       ha="right", va="center")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_ylim(0, 1.5)

    # Panel: total BES (bearing_score)
    ax = axes[panel_idx]
    panel_idx += 1
    centers = (df_r["start"].to_numpy() + df_r["end"].to_numpy()) / 2
    bes = df_r["bearing_score"].to_numpy()
    pos = bes >= 0
    ax.bar(centers[pos], bes[pos], width=args.bin_width,
           color="#404040", edgecolor="none", alpha=0.9)
    ax.bar(centers[~pos], bes[~pos], width=args.bin_width,
           color="#bbbbbb", edgecolor="none", alpha=0.9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylabel(f"BES\n({args.label_a}↑ vs {args.label_b}↓)",
                   fontsize=8, rotation=0, ha="right", va="center")

    # FDR-significant bin markers
    sig_mask = df_r["significant_fdr0.05"].astype(int) == 1
    if sig_mask.any():
        ax.scatter(centers[sig_mask],
                    np.sign(bes[sig_mask]) * (np.abs(bes[sig_mask]) + 0.3),
                    marker="*", color="gold", s=40, edgecolor="black",
                    linewidth=0.5, zorder=10, label="FDR<0.05")
        ax.legend(loc="upper right", fontsize=7)

    # Panel per track
    for j, (t, col, color, idx) in enumerate(zip(
            plot_tracks, kl_cols_plot, track_colors, track_indices)):
        ax = axes[panel_idx]
        panel_idx += 1
        vals = df_r[col].to_numpy()
        nonzero = vals != 0
        pos = nonzero & (vals > 0)
        neg = nonzero & (vals < 0)

        # Color: track color for positive (A-direction), darker shade for B
        ax.bar(centers[pos], vals[pos], width=args.bin_width,
               color=color, edgecolor="none", alpha=0.9)
        ax.bar(centers[neg], vals[neg], width=args.bin_width,
               color=color, edgecolor="black", linewidth=0.3, alpha=0.5)
        ax.axhline(0, color="black", linewidth=0.5)

        # Highlight strong bins
        strong = np.abs(vals) >= args.highlight_threshold
        if strong.any():
            ax.scatter(centers[strong],
                        np.sign(vals[strong]) * (np.abs(vals[strong]) + 0.2),
                        marker="v", color="red", s=20, zorder=10)

        ax.set_ylabel(f"{t}\nkl", fontsize=8, rotation=0,
                        ha="right", va="center", color=color, fontweight="bold")

        # Hide x-axis on all but last panel
        if panel_idx < n_panels:
            ax.tick_params(axis="x", labelbottom=False)

    # X-axis on bottom panel
    axes[-1].set_xlim(start, end)
    axes[-1].set_xlabel(f"{chrom} position (bp)", fontsize=9)
    axes[-1].ticklabel_format(useOffset=False, style="plain", axis="x")

    fig.suptitle(f"BEARING KL decomposition along {chrom}:{start:,}-{end:,}\n"
                  f"{args.label_a} vs {args.label_b}",
                  fontsize=11, y=0.995)
    plt.tight_layout()
    plt.savefig(args.out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
