#!/usr/bin/env python3
"""
bearing_decomposition_report.py

Cross-comparison decomposition report for BEARING differential analyses.

Given a directory of `diff_<A>_vs_<B>.stats.tsv` files (one per pairwise
condition comparison, output by bearing_pvalue.py --diff after a
compare_qcat.py + permutation null run), this script summarizes:

  1. For each comparison: how many bins are "gained in A" vs "gained in B"
     at two thresholds:
       (a) FDR < 0.05 (significant_fdr0.05 column == 1)
       (b) Top-N by |bearing_score| (default N=500 per direction)

  2. For each bin: which of the K input tracks dominated the per-bin score
     (largest |kl_i| where i = track index). This identifies "driver tracks"
     for each gained-direction × comparison cell.

  3. For each comparison: signed per-track contributions (mean kl_i across
     gained-A bins; mean kl_i across gained-B bins). Reveals systematic
     compositional shifts, e.g., "DN vs DP is driven by Cohesin gain in DP".

Outputs:
  - <out_prefix>.pdf            multi-page summary report
  - <out_prefix>.bin_counts.tsv per-comparison bin counts (machine-readable)
  - <out_prefix>.driver_counts.tsv driver-track tally per comparison
  - <out_prefix>.track_contrib.tsv per-comparison per-track mean kl_i

Usage:
  python bearing_decomposition_report.py \\
      --diff-dir results_v6 \\
      --categories DN_rep1_cats.json \\
      --out cross_comparison_report \\
      [--top-n 500] [--fdr-col significant_fdr0.05]

The --categories argument is a per-sample _cats.json file (any sample
suffices since all samples share the same category list in BEARING).
The categories define track names and colors for plotting.
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


# -----------------------------------------------------------------------
# Inputs / parsing
# -----------------------------------------------------------------------

def parse_categories(cats_json_path):
    """Read a *_cats.json file and return ordered list of (name, color)."""
    with open(cats_json_path) as f:
        d = json.load(f)
    cats = d["categories"]
    # Keys are 1-indexed string integers
    ordered = []
    for k in sorted(cats.keys(), key=int):
        name, color = cats[k]
        ordered.append((name, color))
    return ordered


def parse_comparison_name(filename):
    """Extract (A, B) from 'diff_<A>_vs_<B>.stats.tsv'."""
    base = Path(filename).name
    base = re.sub(r"\.stats\.tsv(\.gz)?$", "", base)
    base = re.sub(r"^diff_", "", base)
    m = re.match(r"^(.+?)_vs_(.+?)$", base)
    if not m:
        raise ValueError(f"Cannot parse comparison name from: {filename}")
    return m.group(1), m.group(2)


def load_diff_tsv(path, n_tracks):
    """Load one diff TSV. Returns DataFrame with parsed kl_* columns."""
    df = pd.read_csv(path, sep="\t")
    # Find kl_* columns; ensure we get exactly n_tracks
    kl_cols = [c for c in df.columns if c.startswith("kl_")]
    kl_cols.sort(key=lambda c: int(c.split("_")[1]))
    if len(kl_cols) != n_tracks:
        sys.exit(
            f"ERROR: {path}: expected {n_tracks} kl_* columns, got "
            f"{len(kl_cols)}: {kl_cols}"
        )
    return df, kl_cols


# -----------------------------------------------------------------------
# Per-comparison metrics
# -----------------------------------------------------------------------

def compute_comparison_metrics(df, kl_cols, track_names,
                               fdr_col="significant_fdr0.05",
                               top_n=500):
    """For one comparison, compute count/driver/contribution metrics
    under both FDR-significance and top-N thresholds.

    Convention: rows with direction == '+' are "gained in A"
                rows with direction == '-' are "gained in B"
    """
    out = {}

    # --- Threshold A: FDR-significant bins
    sig_mask = (df[fdr_col].astype(int) == 1) if fdr_col in df.columns else \
        pd.Series([False] * len(df))
    sig_df = df[sig_mask]

    # --- Threshold B: top-N by |bearing_score| in each direction
    # We take top-N per direction (gained-A and gained-B) so the two
    # subsets are matched in size.
    by_dir = {"+": df[df["direction"] == "+"].nlargest(top_n, "bearing_score_tested"),
              "-": df[df["direction"] == "-"].nlargest(top_n, "bearing_score_tested")}
    topn_df = pd.concat([by_dir["+"], by_dir["-"]], ignore_index=False)

    for label, subset in [("fdr", sig_df), ("topN", topn_df)]:
        for direction_label, dir_sign in [("A", "+"), ("B", "-")]:
            sub = subset[subset["direction"] == dir_sign]
            d = {"n_bins": len(sub)}
            if len(sub) > 0:
                kl_mat = sub[kl_cols].to_numpy()
                # Driver track per bin = argmax of |kl_i|
                abs_kl = np.abs(kl_mat)
                drivers = np.argmax(abs_kl, axis=1)
                driver_counts = np.bincount(drivers, minlength=len(kl_cols))
                d["drivers"] = dict(zip(track_names, driver_counts.tolist()))
                d["driver_frac"] = dict(zip(
                    track_names,
                    (driver_counts / max(len(sub), 1)).tolist()
                ))
                # Mean signed kl per track (can be near-zero if signs cancel)
                d["mean_kl"] = dict(zip(track_names,
                                         kl_mat.mean(axis=0).tolist()))
                # Median signed kl per track (more robust)
                d["median_kl"] = dict(zip(track_names,
                                           np.median(kl_mat, axis=0).tolist()))
                # Mean |kl_i| -- magnitude regardless of direction. Better
                # measure of "how much each track contributes" when sub
                # contains mixed-direction bins.
                d["mean_abs_kl"] = dict(zip(track_names,
                                              abs_kl.mean(axis=0).tolist()))
                # Fraction of total |BES| explained by each track.
                # = sum(|kl_i|) / sum(|BES|), where sum(|BES|) = sum(sum(|kl|))
                # over bins. This is the cleanest "which track drives this
                # comparison" metric.
                total_abs_bes = abs_kl.sum()
                if total_abs_bes > 0:
                    track_frac = abs_kl.sum(axis=0) / total_abs_bes
                else:
                    track_frac = np.zeros(len(kl_cols))
                d["frac_of_abs_bes"] = dict(zip(track_names,
                                                  track_frac.tolist()))
            else:
                d["drivers"] = dict(zip(track_names,
                                         [0] * len(track_names)))
                d["driver_frac"] = dict(zip(track_names,
                                              [0.0] * len(track_names)))
                d["mean_kl"] = dict(zip(track_names,
                                          [0.0] * len(track_names)))
                d["median_kl"] = dict(zip(track_names,
                                            [0.0] * len(track_names)))
                d["mean_abs_kl"] = dict(zip(track_names,
                                              [0.0] * len(track_names)))
                d["frac_of_abs_bes"] = dict(zip(track_names,
                                                  [0.0] * len(track_names)))
            out[(label, direction_label)] = d

    return out


# -----------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------

def plot_bin_count_heatmap(metrics, comparisons, ax, title):
    """Heatmap: rows = comparisons, cols = [A-fdr, B-fdr, A-topN, B-topN]"""
    rows = [f"{a}\nvs\n{b}" for a, b in comparisons]
    col_labels = ["A gained\n(FDR)", "B gained\n(FDR)",
                  "A gained\n(top-N)", "B gained\n(top-N)"]
    data = np.zeros((len(comparisons), 4), dtype=float)
    for i, (a, b) in enumerate(comparisons):
        m = metrics[(a, b)]
        data[i, 0] = m[("fdr", "A")]["n_bins"]
        data[i, 1] = m[("fdr", "B")]["n_bins"]
        data[i, 2] = m[("topN", "A")]["n_bins"]
        data[i, 3] = m[("topN", "B")]["n_bins"]

    # Log color for the wide dynamic range
    with np.errstate(divide="ignore"):
        log_data = np.log10(data + 1)
    im = ax.imshow(log_data, aspect="auto", cmap="viridis")
    ax.set_xticks(range(4))
    ax.set_xticklabels(col_labels, fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=8)
    ax.set_title(title)
    # Annotate cells with counts
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = int(data[i, j])
            c = "white" if log_data[i, j] < log_data.max() * 0.5 else "black"
            ax.text(j, i, f"{v:,}", ha="center", va="center", color=c,
                    fontsize=7)
    plt.colorbar(im, ax=ax, label="log10(n_bins+1)")


def plot_driver_stacked(metrics, comparisons, track_names, track_colors,
                         label, ax, title):
    """Stacked bars: each comparison has two bars (A-gained, B-gained).
    Each bar segment = driver-track count.
    """
    n = len(comparisons)
    x = np.arange(n)
    width = 0.4

    # Build matrices: rows = comparisons, cols = tracks
    A_mat = np.zeros((n, len(track_names)))
    B_mat = np.zeros((n, len(track_names)))
    for i, (a, b) in enumerate(comparisons):
        m_a = metrics[(a, b)][(label, "A")]["drivers"]
        m_b = metrics[(a, b)][(label, "B")]["drivers"]
        for j, t in enumerate(track_names):
            A_mat[i, j] = m_a.get(t, 0)
            B_mat[i, j] = m_b.get(t, 0)

    # Stacked bars
    bottoms_a = np.zeros(n)
    bottoms_b = np.zeros(n)
    handles = []
    for j, t in enumerate(track_names):
        ba = ax.bar(x - width / 2, A_mat[:, j], width, bottom=bottoms_a,
                    color=track_colors[j], edgecolor="white", linewidth=0.5,
                    label=t)
        ax.bar(x + width / 2, B_mat[:, j], width, bottom=bottoms_b,
               color=track_colors[j], edgecolor="white", linewidth=0.5)
        handles.append(ba)
        bottoms_a += A_mat[:, j]
        bottoms_b += B_mat[:, j]

    # X-axis labels
    labels = [f"{a}\nvs {b}" for a, b in comparisons]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Driver bin count")
    ax.set_title(title)
    ax.legend(handles=handles, labels=list(track_names), loc="upper right",
              fontsize=7, ncol=2, framealpha=0.9)

    # Annotate which bar is A vs B
    for i in range(n):
        ax.text(i - width / 2, bottoms_a[i] + 3, "A",
                ha="center", va="bottom", fontsize=7, color="black",
                fontweight="bold")
        ax.text(i + width / 2, bottoms_b[i] + 3, "B",
                ha="center", va="bottom", fontsize=7, color="black",
                fontweight="bold")


def plot_track_contrib_heatmap(metrics, comparisons, track_names,
                                label, ax, title):
    """Heatmap of mean signed kl_i per track per comparison direction.

    Rows = comparisons
    Cols = track × direction (12 cols for 6 tracks)
    Color = mean signed kl_i (positive = A gained more of that track,
            negative = B gained more)

    CAVEAT: this can be near-zero for many tracks even when those tracks
    contribute substantially, because signed contributions cancel across
    bins. For "which tracks drive this comparison?" use the frac_of_abs_bes
    metric instead (plot_track_frac_heatmap).
    """
    n = len(comparisons)
    K = len(track_names)
    data = np.zeros((n, K * 2))
    col_labels = []
    for j, t in enumerate(track_names):
        col_labels.extend([f"{t}\n(A)", f"{t}\n(B)"])

    for i, (a, b) in enumerate(comparisons):
        m_a = metrics[(a, b)][(label, "A")]["mean_kl"]
        m_b = metrics[(a, b)][(label, "B")]["mean_kl"]
        for j, t in enumerate(track_names):
            data[i, j * 2] = m_a.get(t, 0)
            data[i, j * 2 + 1] = m_b.get(t, 0)

    vmax = np.nanmax(np.abs(data)) if data.size > 0 else 1.0
    vmax = max(vmax, 0.01)
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(K * 2))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=7)
    rows = [f"{a} vs {b}" for a, b in comparisons]
    ax.set_yticks(range(n))
    ax.set_yticklabels(rows, fontsize=8)
    ax.set_title(title)
    # Annotate cells
    for i in range(n):
        for j in range(K * 2):
            v = data[i, j]
            if abs(v) > 0.05:
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=6, color="black" if abs(v) < vmax * 0.5
                        else "white")
    plt.colorbar(im, ax=ax, label="Mean signed kl_i\n(+ favors A, - favors B)")


def plot_track_frac_heatmap(metrics, comparisons, track_names,
                             label, ax, title):
    """Heatmap of fraction of total |BES| explained by each track in each
    comparison direction.

    Rows = comparisons
    Cols = track × direction (12 cols for 6 tracks)
    Color = fraction of sum(|BES|) attributable to track i (0-1)

    Each row sums to 2.0 (1.0 for A direction + 1.0 for B direction).
    This is the cleanest "which tracks drive this comparison?" metric
    because it doesn't suffer from sign-cancellation across bins.
    """
    n = len(comparisons)
    K = len(track_names)
    data = np.zeros((n, K * 2))
    col_labels = []
    for j, t in enumerate(track_names):
        col_labels.extend([f"{t}\n(A)", f"{t}\n(B)"])

    for i, (a, b) in enumerate(comparisons):
        m_a = metrics[(a, b)][(label, "A")]["frac_of_abs_bes"]
        m_b = metrics[(a, b)][(label, "B")]["frac_of_abs_bes"]
        for j, t in enumerate(track_names):
            data[i, j * 2] = m_a.get(t, 0)
            data[i, j * 2 + 1] = m_b.get(t, 0)

    # Color: 0 -> white, 1/K (equal share) -> mid, 1.0 -> deep
    vmax = max(0.5, np.nanmax(data)) if data.size > 0 else 0.5
    im = ax.imshow(data, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=vmax)
    ax.set_xticks(range(K * 2))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=7)
    rows = [f"{a} vs {b}" for a, b in comparisons]
    ax.set_yticks(range(n))
    ax.set_yticklabels(rows, fontsize=8)
    ax.set_title(title + "\n(each cell = fraction of sum-of-|BES| from this track)")
    # Annotate cells
    equal_share = 1.0 / K
    for i in range(n):
        for j in range(K * 2):
            v = data[i, j]
            if v > 0.02:
                color = "black" if v < vmax * 0.6 else "white"
                marker = "*" if v > equal_share * 1.5 else ""
                ax.text(j, i, f"{v:.2f}{marker}", ha="center", va="center",
                        fontsize=6, color=color)
    plt.colorbar(im, ax=ax,
                  label=f"Fraction of |BES| from track\n"
                        f"(equal share={equal_share:.2f}; * = >1.5× equal share)")


def plot_per_track_violin(diff_data, track_names, track_colors, ax):
    """One violin per track showing the distribution of |kl_i| across all
    significant bins across all comparisons. Identifies tracks with broad
    vs narrow contribution ranges.
    """
    parts = []
    for j, t in enumerate(track_names):
        all_vals = []
        for (a, b), df_kl_dict in diff_data.items():
            df, kl_cols = df_kl_dict
            sig_df = df[df.get("significant_fdr0.05", pd.Series(
                [0] * len(df))).astype(int) == 1]
            if len(sig_df):
                all_vals.append(np.abs(sig_df[kl_cols[j]].to_numpy()))
        if all_vals:
            parts.append(np.concatenate(all_vals))
        else:
            parts.append(np.array([0.0]))

    # Plot
    positions = np.arange(len(track_names))
    vp = ax.violinplot(parts, positions=positions, widths=0.7,
                        showmedians=True, showextrema=False)
    for body, color in zip(vp["bodies"], track_colors):
        body.set_facecolor(color)
        body.set_alpha(0.7)
        body.set_edgecolor("black")
    ax.set_xticks(positions)
    ax.set_xticklabels(track_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("|kl_i| across FDR-significant bins\n(all comparisons)")
    ax.set_title("Per-track contribution distribution\n"
                 "(width = how often that track drives sig bins)")


def plot_comparison_legend(track_names, track_colors, ax):
    """Standalone legend panel — color = track, label = name."""
    ax.axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, edgecolor="black")
               for c in track_colors]
    ax.legend(handles, track_names, loc="center", fontsize=10,
              ncol=2, title="Track colors", title_fontsize=10)


# -----------------------------------------------------------------------
# Top-level orchestration
# -----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--diff-dir", required=True,
                    help="Directory containing diff_*.stats.tsv files")
    ap.add_argument("--categories", required=True,
                    help="Path to a *_cats.json file (any sample's; "
                         "defines track names and colors)")
    ap.add_argument("--out", required=True,
                    help="Output prefix (creates <out>.pdf, "
                         "<out>.bin_counts.tsv, etc.)")
    ap.add_argument("--top-n", type=int, default=500,
                    help="Top-N bins per direction by |bearing_score| "
                         "(default 500)")
    ap.add_argument("--fdr-col", default="significant_fdr0.05",
                    help="Column name for FDR-significance indicator")
    ap.add_argument("--filename-pattern", default="diff_*.stats.tsv",
                    help="Glob pattern within --diff-dir")
    args = ap.parse_args()

    # Parse track definitions
    cats = parse_categories(args.categories)
    track_names = [c[0] for c in cats]
    track_colors = [c[1] for c in cats]
    n_tracks = len(track_names)
    print(f"Tracks ({n_tracks}): {track_names}", flush=True)

    # Enumerate comparisons
    diff_files = sorted(Path(args.diff_dir).glob(args.filename_pattern))
    if not diff_files:
        sys.exit(f"ERROR: no diff files matched in {args.diff_dir}")
    print(f"Found {len(diff_files)} diff files", flush=True)

    # Load all
    diff_data = {}     # (A, B) -> (df, kl_cols)
    metrics = {}       # (A, B) -> per-comparison dict
    comparisons = []   # ordered list of (A, B) pairs

    for f in diff_files:
        try:
            a, b = parse_comparison_name(str(f))
        except ValueError as ex:
            print(f"  WARN: skipping {f}: {ex}", file=sys.stderr)
            continue
        print(f"  Loading {a} vs {b} from {f}...", flush=True)
        df, kl_cols = load_diff_tsv(str(f), n_tracks)
        diff_data[(a, b)] = (df, kl_cols)
        m = compute_comparison_metrics(df, kl_cols, track_names,
                                        fdr_col=args.fdr_col,
                                        top_n=args.top_n)
        metrics[(a, b)] = m
        comparisons.append((a, b))

    # --- Write TSV summaries ---
    out_prefix = args.out

    # 1) Bin counts table
    rows = []
    for a, b in comparisons:
        for label in ("fdr", "topN"):
            for dir_label in ("A", "B"):
                d = metrics[(a, b)][(label, dir_label)]
                rows.append({
                    "comparison_A": a,
                    "comparison_B": b,
                    "threshold": label,
                    "gained_in": dir_label,
                    "n_bins": d["n_bins"],
                })
    pd.DataFrame(rows).to_csv(f"{out_prefix}.bin_counts.tsv",
                               sep="\t", index=False)
    print(f"Wrote {out_prefix}.bin_counts.tsv", flush=True)

    # 2) Driver counts
    rows = []
    for a, b in comparisons:
        for label in ("fdr", "topN"):
            for dir_label in ("A", "B"):
                d = metrics[(a, b)][(label, dir_label)]
                for t, n in d["drivers"].items():
                    rows.append({
                        "comparison_A": a,
                        "comparison_B": b,
                        "threshold": label,
                        "gained_in": dir_label,
                        "track": t,
                        "n_bins_driven_by_track": n,
                    })
    pd.DataFrame(rows).to_csv(f"{out_prefix}.driver_counts.tsv",
                               sep="\t", index=False)
    print(f"Wrote {out_prefix}.driver_counts.tsv", flush=True)

    # 3) Per-track contribution means
    rows = []
    for a, b in comparisons:
        for label in ("fdr", "topN"):
            for dir_label in ("A", "B"):
                d = metrics[(a, b)][(label, dir_label)]
                for t in track_names:
                    rows.append({
                        "comparison_A": a,
                        "comparison_B": b,
                        "threshold": label,
                        "gained_in": dir_label,
                        "track": t,
                        "mean_kl": d["mean_kl"][t],
                        "median_kl": d["median_kl"][t],
                        "mean_abs_kl": d["mean_abs_kl"][t],
                        "frac_of_abs_bes": d["frac_of_abs_bes"][t],
                        "driver_frac": d["driver_frac"][t],
                    })
    pd.DataFrame(rows).to_csv(f"{out_prefix}.track_contrib.tsv",
                               sep="\t", index=False)
    print(f"Wrote {out_prefix}.track_contrib.tsv", flush=True)

    # --- Build PDF report ---
    pdf_path = f"{out_prefix}.pdf"
    with PdfPages(pdf_path) as pdf:

        # Page 1: bin counts heatmap
        fig, axes = plt.subplots(1, 2, figsize=(16, max(6, 0.5 * len(comparisons))))
        plot_bin_count_heatmap(metrics, comparisons, axes[0],
                                "Bin counts: FDR-sig and top-N")
        plot_comparison_legend(track_names, track_colors, axes[1])
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 2: HEADLINE - fraction of |BES| explained by each track (FDR)
        fig, ax = plt.subplots(figsize=(14, max(6, 0.5 * len(comparisons))))
        plot_track_frac_heatmap(metrics, comparisons, track_names,
                                  "fdr", ax,
                                  "Track contribution to compositional shift\n"
                                  "(FDR-significant bins)")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 3: HEADLINE - fraction of |BES| explained by each track (top-N)
        fig, ax = plt.subplots(figsize=(14, max(6, 0.5 * len(comparisons))))
        plot_track_frac_heatmap(metrics, comparisons, track_names,
                                  "topN", ax,
                                  f"Track contribution to compositional shift\n"
                                  f"(top-{args.top_n} bins per direction)")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 4: driver-track stacked bars (FDR)
        fig, ax = plt.subplots(figsize=(14, 6))
        plot_driver_stacked(metrics, comparisons, track_names, track_colors,
                             "fdr", ax,
                             "Driver-track decomposition: FDR-significant bins\n"
                             "(A=A-gained-bar, B=B-gained-bar; color=dominant track)")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 5: driver-track stacked bars (top-N)
        fig, ax = plt.subplots(figsize=(14, 6))
        plot_driver_stacked(metrics, comparisons, track_names, track_colors,
                             "topN", ax,
                             f"Driver-track decomposition: top-{args.top_n} "
                             f"bins per direction by |BES|")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 6: per-track violin distribution
        fig, ax = plt.subplots(figsize=(10, 6))
        plot_per_track_violin(diff_data, track_names, track_colors, ax)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Page 7-8: mean signed kl (moved to end - less useful diagnostic)
        fig, ax = plt.subplots(figsize=(14, max(6, 0.5 * len(comparisons))))
        plot_track_contrib_heatmap(metrics, comparisons, track_names,
                                     "fdr", ax,
                                     "Mean signed kl_i per track (FDR bins)\n"
                                     "CAVEAT: near-zero values may reflect "
                                     "sign cancellation across bins,\n"
                                     "not lack of contribution. See pages 2-3 for "
                                     "track-fraction summary.")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(14, max(6, 0.5 * len(comparisons))))
        plot_track_contrib_heatmap(metrics, comparisons, track_names,
                                     "topN", ax,
                                     f"Mean signed kl_i per track "
                                     f"(top-{args.top_n} bins)\n"
                                     f"CAVEAT: near-zero values may reflect "
                                     f"sign cancellation across bins")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"Wrote {pdf_path}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
