#!/usr/bin/env python3
"""
bearing_tad_decomposition.py

Compute and plot BEARING compositional decomposition aggregated by TAD
(Topologically Associating Domain) boundaries from Hi-C data.

For each TAD, computes mean |BES|, peak |BES|, FDR rate, and per-track
contribution fractions. Plots TAD bars + insulation profile + per-comparison
BES and FDR heatmaps with one cell per TAD.

Usage:
  python bearing_tad_decomposition.py \\
      --diff-dir results_v6 \\
      --categories DN_rep1_cats.json \\
      --tads tads/merged_corrected_KR_DN_bs_25000_tads.bed \\
      --insulation tads/merged_corrected_KR_DN_bs_25000_tad_score.bm \\
      --region chr6:40000000-47400000 \\
      --annotations tcrb_wide_annotations.bed \\
      --out tcrb_wide_tad_DN

Outputs:
  <out>.pdf   Multi-panel TAD decomposition plot
  <out>.tsv   Per-TAD per-comparison stats (machine-readable)
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


def load_tads(bed_path, chrom_filter=None, start_filter=None,
              end_filter=None):
    tads = []
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
            tad_id = parts[3] if len(parts) >= 4 else f"tad_{len(tads)+1}"
            if chrom_filter and chrom != chrom_filter:
                continue
            if start_filter is not None and end < start_filter:
                continue
            if end_filter is not None and start > end_filter:
                continue
            tads.append((chrom, start, end, tad_id))
    return tads


def load_insulation(bm_path, chrom_filter=None):
    """Load insulation scores from a HiCExplorer .bm file.
    Returns dict: chrom -> list of (start, end, score).
    Uses the second score column if available (wider depth).
    """
    insulation = {}
    with open(bm_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            if chrom_filter and chrom != chrom_filter:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                score = float(parts[4]) if len(parts) >= 5 \
                    else float(parts[3])
            except ValueError:
                continue
            insulation.setdefault(chrom, []).append((start, end, score))
    return insulation


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
            name = parts[3] if len(parts) >= 4 else f"r_{len(regions)+1}"
            regions.append((chrom, start, end, name))
    return regions


def load_diff_tsv(path, n_tracks):
    df = pd.read_csv(path, sep="\t")
    kl_cols = sorted([c for c in df.columns if c.startswith("kl_")],
                     key=lambda c: int(c.split("_")[1]))
    if len(kl_cols) != n_tracks:
        sys.exit(f"ERROR: {path}: expected {n_tracks} kl_* cols, "
                 f"got {len(kl_cols)}")
    return df, kl_cols


def compute_tad_stats(df, kl_cols, track_names, tads,
                       fdr_col="significant_fdr0.05"):
    rows = []
    for chrom, t_start, t_end, tad_id in tads:
        mask = ((df["chrom"] == chrom) &
                (df["end"] > t_start) & (df["start"] < t_end))
        sub = df[mask]
        n_bins = len(sub)

        row = {
            "chrom": chrom,
            "tad_start": t_start,
            "tad_end": t_end,
            "tad_id": tad_id,
            "tad_size_bp": t_end - t_start,
            "n_bins": n_bins,
        }
        if n_bins == 0:
            row["mean_abs_bes"] = 0.0
            row["max_abs_bes"] = 0.0
            row["n_fdr_sig"] = 0
            row["frac_fdr_sig"] = 0.0
            for t in track_names:
                row[f"frac_{t}"] = 0.0
                row[f"mean_abs_kl_{t}"] = 0.0
            rows.append(row)
            continue

        abs_bes = sub["bearing_score_tested"].to_numpy()
        kl_mat = sub[kl_cols].to_numpy()
        abs_kl = np.abs(kl_mat)
        total_abs = abs_kl.sum()

        n_fdr = (int(sub[fdr_col].astype(int).sum())
                  if fdr_col in sub.columns else 0)

        row["mean_abs_bes"] = float(abs_bes.mean())
        row["max_abs_bes"] = float(abs_bes.max())
        row["n_fdr_sig"] = n_fdr
        row["frac_fdr_sig"] = n_fdr / n_bins
        track_frac = (abs_kl.sum(axis=0) / total_abs
                       if total_abs > 0 else np.zeros(len(kl_cols)))
        for j, t in enumerate(track_names):
            row[f"frac_{t}"] = float(track_frac[j])
            row[f"mean_abs_kl_{t}"] = float(abs_kl[:, j].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def plot_tad_overview(stats_dfs_by_comp, tads, insulation, annotations,
                      track_names, track_colors,
                      region_start, region_end, chrom):
    """Multi-panel figure: annotations + TAD bars + insulation +
    per-TAD BES heatmap + per-TAD FDR heatmap.
    """
    comparisons = sorted(stats_dfs_by_comp.keys())
    n_comp = len(comparisons)

    has_insulation = bool(insulation) and chrom in insulation
    n_panels = 4 + (1 if has_insulation else 0)
    height_ratios = [0.5, 0.6]  # annotations, TAD bars
    if has_insulation:
        height_ratios.append(1.0)
    height_ratios.extend([n_comp * 0.35, n_comp * 0.35])

    fig = plt.figure(figsize=(14, max(6, sum(height_ratios) * 0.8)))
    gs = fig.add_gridspec(n_panels, 1, height_ratios=height_ratios,
                           hspace=0.15)

    ax_anno = fig.add_subplot(gs[0])
    ax_tad = fig.add_subplot(gs[1], sharex=ax_anno)
    if has_insulation:
        ax_ins = fig.add_subplot(gs[2], sharex=ax_anno)
        next_idx = 3
    else:
        ax_ins = None
        next_idx = 2
    ax_bes = fig.add_subplot(gs[next_idx], sharex=ax_anno)
    ax_fdr = fig.add_subplot(gs[next_idx + 1], sharex=ax_anno)

    # annotations
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
                              xytext=((s2 + e2) / 2, 1.05),
                              ha="center", va="bottom", fontsize=7,
                              fontweight="bold",
                              arrowprops=dict(arrowstyle="-",
                                               color="navy",
                                               linewidth=0.5))
        else:
            ax_anno.text((s2 + e2) / 2, 0.5, n2,
                          ha="center", va="center",
                          fontsize=8, fontweight="bold")
    ax_anno.set_yticks([])
    ax_anno.set_ylabel("Features", fontsize=8, rotation=0,
                        ha="right", va="center")
    for spine in ["top", "right", "left"]:
        ax_anno.spines[spine].set_visible(False)
    ax_anno.set_ylim(0, 1.5)
    ax_anno.tick_params(axis="x", labelbottom=False)

    # TAD bars
    in_window_tads = [t for t in tads
                       if t[0] == chrom and t[2] > region_start
                       and t[1] < region_end]
    ax_tad.set_ylim(0, 1)
    for i, (_, t_start, t_end, _) in enumerate(in_window_tads):
        color = "#7799cc" if i % 2 == 0 else "#88cc88"
        ax_tad.add_patch(plt.Rectangle(
            (t_start, 0.3), t_end - t_start, 0.4,
            facecolor=color, edgecolor="black", linewidth=0.4, alpha=0.75
        ))
    ax_tad.set_yticks([])
    ax_tad.set_ylabel(f"TADs (n={len(in_window_tads)})", fontsize=8,
                        rotation=0, ha="right", va="center")
    ax_tad.tick_params(axis="x", labelbottom=False)

    # insulation profile
    if ax_ins is not None:
        ins_data = insulation[chrom]
        xs, ys = [], []
        for s, e, sc in ins_data:
            if e < region_start or s > region_end:
                continue
            xs.append((s + e) / 2)
            ys.append(sc)
        if xs:
            xs_a, ys_a = np.array(xs), np.array(ys)
            ax_ins.plot(xs_a, ys_a, color="#444", linewidth=0.8)
            ax_ins.fill_between(xs_a, 0, ys_a,
                                  where=(ys_a < 0),
                                  color="#cc6600", alpha=0.3)
            ax_ins.axhline(0, color="gray", linestyle="--",
                            linewidth=0.5, alpha=0.5)
        ax_ins.set_ylabel("Insulation", fontsize=8, rotation=0,
                           ha="right", va="center")
        ax_ins.tick_params(axis="x", labelbottom=False)

    # Per-comparison heatmaps
    if not in_window_tads:
        ax_bes.text(0.5, 0.5, "No TADs in window",
                     ha="center", va="center", transform=ax_bes.transAxes)
        ax_fdr.text(0.5, 0.5, "No TADs in window",
                     ha="center", va="center", transform=ax_fdr.transAxes)
    else:
        data_bes = np.zeros((n_comp, len(in_window_tads)))
        data_fdr = np.zeros((n_comp, len(in_window_tads)))
        for i, comp in enumerate(comparisons):
            stats_df = stats_dfs_by_comp[comp]
            for j, (_, t_start, t_end, _) in enumerate(in_window_tads):
                row = stats_df[
                    (stats_df["tad_start"] == t_start) &
                    (stats_df["tad_end"] == t_end)
                ]
                if len(row) > 0:
                    data_bes[i, j] = row.iloc[0]["mean_abs_bes"]
                    data_fdr[i, j] = row.iloc[0]["frac_fdr_sig"]

        # BES heatmap
        vmax_bes = max(data_bes.max(), 0.5)
        for j, (_, t_start, t_end, _) in enumerate(in_window_tads):
            for i in range(n_comp):
                ax_bes.add_patch(plt.Rectangle(
                    (t_start, n_comp - i - 1),
                    t_end - t_start, 1,
                    facecolor=plt.cm.magma_r(data_bes[i, j] / vmax_bes),
                    edgecolor="white", linewidth=0.3
                ))
        ax_bes.set_xlim(region_start, region_end)
        ax_bes.set_ylim(0, n_comp)
        ax_bes.set_yticks(np.arange(n_comp) + 0.5)
        ax_bes.set_yticklabels(
            [f"{a} vs {b}" for (a, b) in reversed(comparisons)],
            fontsize=8
        )
        ax_bes.set_ylabel("Per-TAD\nmean |BES|", fontsize=9)
        ax_bes.tick_params(axis="x", labelbottom=False)

        # FDR heatmap
        vmax_fdr = max(data_fdr.max(), 0.01)
        for j, (_, t_start, t_end, _) in enumerate(in_window_tads):
            for i in range(n_comp):
                ax_fdr.add_patch(plt.Rectangle(
                    (t_start, n_comp - i - 1),
                    t_end - t_start, 1,
                    facecolor=plt.cm.Oranges(data_fdr[i, j] / vmax_fdr),
                    edgecolor="white", linewidth=0.3
                ))
        ax_fdr.set_xlim(region_start, region_end)
        ax_fdr.set_ylim(0, n_comp)
        ax_fdr.set_yticks(np.arange(n_comp) + 0.5)
        ax_fdr.set_yticklabels(
            [f"{a} vs {b}" for (a, b) in reversed(comparisons)],
            fontsize=8
        )
        ax_fdr.set_ylabel("Per-TAD\nFDR rate", fontsize=9)

    ax_fdr.set_xlim(region_start, region_end)
    ax_fdr.set_xlabel(f"{chrom} position (bp)", fontsize=10)
    ax_fdr.ticklabel_format(useOffset=False, style="plain", axis="x")

    fig.suptitle(
        f"BEARING per-TAD decomposition: {chrom}:{region_start:,}-{region_end:,}\n"
        f"n_TADs in window = {len(in_window_tads)}",
        fontsize=11, y=0.995
    )
    return fig


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--diff-dir", required=True,
                    help="Directory with diff_*.stats.tsv files")
    ap.add_argument("--categories", required=True,
                    help="A *_cats.json file")
    ap.add_argument("--tads", required=True,
                    help="TAD BED file from HiCExplorer")
    ap.add_argument("--insulation", default=None,
                    help="HiCExplorer insulation .bm file")
    ap.add_argument("--region", required=True, type=parse_region_str,
                    help="Region in format chrom:start-end")
    ap.add_argument("--annotations", default=None,
                    help="Feature annotation BED")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("--filename-pattern", default="diff_*.stats.tsv")
    args = ap.parse_args()

    cats = parse_categories(args.categories)
    track_names = [c[0] for c in cats]
    track_colors = [c[1] for c in cats]
    n_tracks = len(track_names)
    print(f"Tracks ({n_tracks}): {track_names}", flush=True)

    chrom, region_start, region_end = args.region
    print(f"Region: {chrom}:{region_start:,}-{region_end:,}", flush=True)

    tads = load_tads(args.tads, chrom_filter=chrom,
                      start_filter=region_start, end_filter=region_end)
    print(f"Loaded {len(tads)} TADs in region", flush=True)

    insulation = {}
    if args.insulation:
        insulation = load_insulation(args.insulation, chrom_filter=chrom)
        n_ins = len(insulation.get(chrom, []))
        print(f"Loaded insulation scores: {n_ins} bins on {chrom}",
              flush=True)

    annotations = load_annotations(args.annotations)
    print(f"Loaded {len(annotations)} feature annotations", flush=True)

    diff_files = sorted(Path(args.diff_dir).glob(args.filename_pattern))
    if not diff_files:
        sys.exit(f"ERROR: no diff files in {args.diff_dir}")
    print(f"\nProcessing {len(diff_files)} diff files", flush=True)

    stats_dfs_by_comp = {}
    all_rows = []

    for diff_file in diff_files:
        try:
            a, b = parse_comparison_name(str(diff_file))
        except ValueError as ex:
            print(f"  WARN skip {diff_file}: {ex}", file=sys.stderr)
            continue
        comparison_label = f"{a} vs {b}"
        print(f"  {comparison_label} ...", flush=True)

        df, kl_cols = load_diff_tsv(str(diff_file), n_tracks)
        df_region = df[
            (df["chrom"] == chrom) &
            (df["end"] > region_start) &
            (df["start"] < region_end)
        ]
        print(f"    {len(df_region)} bins in region", flush=True)

        stats_df = compute_tad_stats(df_region, kl_cols, track_names, tads)
        stats_df["comparison_A"] = a
        stats_df["comparison_B"] = b
        stats_dfs_by_comp[(a, b)] = stats_df
        all_rows.append(stats_df)

    combined = pd.concat(all_rows, ignore_index=True)
    tsv_path = f"{args.out}.tsv"
    combined.to_csv(tsv_path, sep="\t", index=False)
    print(f"\nWrote {tsv_path}", flush=True)

    pdf_path = f"{args.out}.pdf"
    with PdfPages(pdf_path) as pdf:
        fig = plot_tad_overview(
            stats_dfs_by_comp, tads, insulation, annotations,
            track_names, track_colors,
            region_start, region_end, chrom
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
