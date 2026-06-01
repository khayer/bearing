#!/usr/bin/env python3
"""
bearing_hic_plot_pval_overlay.py
================================
Combined epilogos + p-value overlay drawing functions.

Drop-in additions for bearing_hic_plot.py. Each function renders:

  1. KL stacked bars with BINARY OPACITY MASK
       - Significant bins (p_adj <= alpha): full opacity (alpha=0.92)
       - Non-significant bins:              low opacity (alpha=0.18)
     Desaturated bins are still visible so the score magnitude is never
     hidden, but the reader immediately sees where significance falls.

  2. -log10(p) LINE + FILL on a TWIN AXIS
       Plotted in a muted purple (#7b4fa0) below the bars.
       The significance cutoff (-log10(alpha)) is shown as a dashed line.
       The fill area is semi-transparent so the KL bars remain the visual
       focus; the line track answers "where is this significant?" without
       competing with the colour-coded biology.

FUNCTIONS PROVIDED
------------------
  draw_epilogos_with_pval_horizontal(
      ax, positions, score_mat, num_states,
      region_start, region_end, categories,
      pval_positions, pval_values,
      pval_alpha=0.05,
      highlights=None, label="", y_max=None,
  )

  draw_epilogos_with_pval_vertical(
      ax, positions, score_mat, num_states,
      region_start, region_end, categories,
      pval_positions, pval_values,
      pval_alpha=0.05,
      highlights=None, label="", y_max=None,
  )

  align_pval_to_qcat(pval_positions, pval_values, qcat_positions)
      Helper to interpolate/snap a BigWig p-value array onto the qcat
      bin positions so they can be used together.

USAGE IN bearing_hic_plot.py
-----------------------------
Replace calls to:
    draw_epilogos_horizontal(ax_epi, pos, scores, ...)
    draw_pval_horizontal(ax_pval, pval_pos, pval_vals, ...)

With a single call to:
    draw_epilogos_with_pval_horizontal(
        ax_epi, pos, scores, num_states,
        region_start, region_end, categories,
        pval_positions=pval_pos,
        pval_values=pval_vals,
        pval_alpha=0.05,
    )

This collapses two separate rows into one row, saving figure height.

DEPENDENCIES
------------
Same as bearing_hic_plot.py (numpy, matplotlib).
"""

import math
import numpy as np


# ── Significance colours ─────────────────────────────────────────────────────
PVAL_LINE_COLOR  = "#5b21b6"   # jewel purple — matches PDF accent colour
PVAL_FILL_ALPHA  = 0.20
PVAL_LINE_WIDTH  = 0.9
SIG_BAR_ALPHA    = 0.92        # opacity for significant bins
INSIG_BAR_ALPHA  = 0.15        # opacity for non-significant bins


# ---------------------------------------------------------------------------
# Helper: align pval BigWig values onto qcat bin positions
# ---------------------------------------------------------------------------

def align_pval_to_qcat(pval_positions, pval_values, qcat_positions):
    """
    Snap/interpolate a -log10(p) BigWig value array onto qcat bin positions.

    pval_positions : 1-D array of BigWig bin starts (from load_bigwig_values)
    pval_values    : 1-D array of -log10(p) values at those positions
    qcat_positions : list of (bin_start, bin_end) tuples from load_qcat_scores

    Returns a 1-D numpy array of length len(qcat_positions) with the
    nearest-neighbour -log10(p) value for each qcat bin.
    """
    if pval_positions is None or len(pval_positions) == 0:
        return np.zeros(len(qcat_positions), dtype=np.float64)

    pval_pos_arr = np.asarray(pval_positions, dtype=np.float64)
    pval_val_arr = np.asarray(pval_values, dtype=np.float64)

    result = np.zeros(len(qcat_positions), dtype=np.float64)
    for i, (bs, be) in enumerate(qcat_positions):
        # Use the bin midpoint for matching
        mid = (bs + be) / 2.0
        idx = np.searchsorted(pval_pos_arr, mid, side="right") - 1
        idx = int(np.clip(idx, 0, len(pval_val_arr) - 1))
        result[i] = pval_val_arr[idx]

    return result


# ---------------------------------------------------------------------------
# Combined horizontal track
# ---------------------------------------------------------------------------

def draw_epilogos_with_pval_horizontal(
        ax, positions, score_mat, num_states,
        region_start, region_end, categories,
        pval_positions=None, pval_values=None,
        pval_alpha=0.05,
        highlights=None, label="", y_max=None):
    """
    Horizontal epilogos stacked-bar track with embedded -log10(p) overlay.

    KL bars: significant bins at full opacity, non-significant at low opacity.
    -log10(p) line + fill: drawn on a twin y-axis (right side, muted purple).

    Parameters
    ----------
    ax              : matplotlib Axes
    positions       : list of (bin_start, bin_end) from load_qcat_scores
    score_mat       : (n_bins, num_states) float array
    num_states      : int
    region_start/end: int genomic coordinates
    categories      : list of (name, color) tuples
    pval_positions  : 1-D array of BigWig bin starts (from load_bigwig_values)
                      or None to skip the overlay
    pval_values     : 1-D array of -log10(p) values
                      or None to skip the overlay
    pval_alpha      : significance threshold (default 0.05)
    highlights      : list of (start, end, hex_color) for shading
    label           : y-axis label
    y_max           : fixed y limit for KL axis
    """
    # Compute -log10(alpha) threshold
    cutoff_neglog10 = -math.log10(pval_alpha) if pval_alpha > 0 else 1.3

    # Align pval onto qcat bins
    has_pval = (pval_positions is not None and pval_values is not None
                and len(pval_positions) > 0)
    if has_pval:
        pval_per_bin = align_pval_to_qcat(
            pval_positions, pval_values, positions)
        sig_mask = pval_per_bin >= cutoff_neglog10
    else:
        pval_per_bin = np.zeros(len(positions), dtype=np.float64)
        sig_mask = np.ones(len(positions), dtype=bool)  # all significant

    # ── KL axis setup ─────────────────────────────────────────────────────
    ax.set_xlim(0, 1)
    _y_max = (y_max if y_max is not None
              else (float(score_mat.max()) * 1.15
                    if score_mat.size > 0 else 1.0))
    ax.set_ylim(0, _y_max)

    # Keep KL bars above the p-value twin axis so bars remain the primary layer.
    ax.set_zorder(2)
    ax.patch.set_alpha(0.0)

    # Highlights (behind bars)
    if highlights:
        for hs, he, hcol in highlights:
            x0 = _gax(max(hs, region_start), region_start, region_end)
            x1 = _gax(min(he, region_end),   region_start, region_end)
            ax.axvspan(x0, x1, color=hcol, alpha=0.18, zorder=0)

    # ── Twin axis: -log10(p) line + fill (background layer) ──────────────
    if has_pval:
        ax2 = ax.twinx()
        ax2.set_zorder(1)
        ax2.set_xlim(0, 1)
        pval_max = max(float(pval_per_bin.max()) * 1.15, cutoff_neglog10 * 1.5)
        ax2.set_ylim(0, pval_max)

        xs = np.array([
            _gax((bs + be) / 2, region_start, region_end)
            for bs, be in positions
        ])

        ax2.fill_between(xs, 0, pval_per_bin,
                         color=PVAL_LINE_COLOR, alpha=PVAL_FILL_ALPHA,
                         linewidth=0, zorder=1)
        ax2.plot(xs, pval_per_bin,
                 color=PVAL_LINE_COLOR, linewidth=PVAL_LINE_WIDTH,
                 zorder=1.1)
        ax2.axhline(cutoff_neglog10,
                    color=PVAL_LINE_COLOR, linestyle="--",
                    linewidth=0.7, alpha=0.6, zorder=1.2)

        ax2.set_ylabel(f"-log10(p)", fontsize=5.5, color=PVAL_LINE_COLOR,
                       labelpad=2)
        ax2.tick_params(axis="y", labelsize=5, colors=PVAL_LINE_COLOR,
                        right=True, left=False)
        ax2.yaxis.set_label_position("right")
        ax2.spines["top"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.spines["bottom"].set_visible(False)
        ax2.spines["right"].set_color(PVAL_LINE_COLOR)
        ax2.spines["right"].set_linewidth(0.5)
        ax2.spines["right"].set_alpha(0.5)

    # ── Draw KL stacked bars with opacity mask ────────────────────────────
    for i, (bs, be) in enumerate(positions):
        x0    = _gax(bs, region_start, region_end)
        x1    = _gax(be, region_start, region_end)
        width = max(x1 - x0, 0.0002)
        row   = score_mat[i]
        alpha = SIG_BAR_ALPHA if sig_mask[i] else INSIG_BAR_ALPHA
        bottom = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val <= 0:
                continue
            color = (categories[si][1] if si < len(categories)
                     else "#cccccc")
            ax.bar(x0, val, width=width, bottom=bottom,
                   color=color, alpha=alpha,
                   align="edge", linewidth=0, zorder=2)
            bottom += val

    # ── KL axis cosmetics ─────────────────────────────────────────────────
    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


# ---------------------------------------------------------------------------
# Combined vertical track
# ---------------------------------------------------------------------------

def draw_epilogos_with_pval_vertical(
        ax, positions, score_mat, num_states,
        region_start, region_end, categories,
        pval_positions=None, pval_values=None,
        pval_alpha=0.05,
        highlights=None, label="", y_max=None):
    """
    Vertical epilogos stacked-bar track with embedded -log10(p) overlay.

    Same design as the horizontal version, rotated 90 degrees to align
    with the Hi-C square (genomic position on Y axis, scores on X axis).

    Significant bins rendered at full opacity; non-significant at low opacity.
    -log10(p) line drawn on a twin x-axis (top side).
    """
    cutoff_neglog10 = -math.log10(pval_alpha) if pval_alpha > 0 else 1.3

    has_pval = (pval_positions is not None and pval_values is not None
                and len(pval_positions) > 0)
    if has_pval:
        pval_per_bin = align_pval_to_qcat(
            pval_positions, pval_values, positions)
        sig_mask = pval_per_bin >= cutoff_neglog10
    else:
        pval_per_bin = np.zeros(len(positions), dtype=np.float64)
        sig_mask = np.ones(len(positions), dtype=bool)

    # ── KL axis setup (Y = genomic, X = score) ───────────────────────────
    ax.set_ylim(1, 0)   # top = region_start (matches Hi-C)
    _x_max = (y_max if y_max is not None
              else (float(score_mat.max()) * 1.15
                    if score_mat.size > 0 else 1.0))
    ax.set_xlim(_x_max, 0)  # bars extend left toward Hi-C square

    # Keep KL bars above the p-value twin axis so bars remain the primary layer.
    ax.set_zorder(2)
    ax.patch.set_alpha(0.0)

    if highlights:
        for hs, he, hcol in highlights:
            y0 = _gax(max(hs, region_start), region_start, region_end)
            y1 = _gax(min(he, region_end),   region_start, region_end)
            ax.axhspan(y0, y1, color=hcol, alpha=0.18, zorder=0)

    # ── Twin axis: -log10(p) line (background layer) ─────────────────────
    if has_pval:
        ax2 = ax.twiny()
        ax2.set_zorder(1)
        ax2.set_ylim(1, 0)  # keep same Y orientation
        pval_max = max(float(pval_per_bin.max()) * 1.15, cutoff_neglog10 * 1.5)
        # Match the KL axis direction: larger values toward the Hi-C square
        # so the vertical overlay is not mirrored relative to stacked bars.
        ax2.set_xlim(pval_max, 0)

        ys = np.array([
            _gax((bs + be) / 2, region_start, region_end)
            for bs, be in positions
        ])

        ax2.fill_betweenx(ys, 0, pval_per_bin,
                          color=PVAL_LINE_COLOR, alpha=PVAL_FILL_ALPHA,
                          linewidth=0, zorder=1)
        ax2.plot(pval_per_bin, ys,
                 color=PVAL_LINE_COLOR, linewidth=PVAL_LINE_WIDTH,
                 zorder=1.1)
        ax2.axvline(cutoff_neglog10,
                    color=PVAL_LINE_COLOR, linestyle="--",
                    linewidth=0.7, alpha=0.6, zorder=1.2)

        ax2.set_xlabel(f"-log10(p)", fontsize=5.5, color=PVAL_LINE_COLOR,
                       labelpad=2)
        ax2.tick_params(axis="x", labelsize=5, colors=PVAL_LINE_COLOR,
                        top=True, bottom=False)
        ax2.xaxis.set_label_position("top")
        ax2.spines["bottom"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_color(PVAL_LINE_COLOR)
        ax2.spines["top"].set_linewidth(0.5)
        ax2.spines["top"].set_alpha(0.5)

    # ── Draw KL stacked bars (horizontal bars, vertical layout) ──────────
    for i, (bs, be) in enumerate(positions):
        y0     = _gax(bs, region_start, region_end)
        y1     = _gax(be, region_start, region_end)
        height = max(y1 - y0, 0.0002)
        row    = score_mat[i]
        alpha  = SIG_BAR_ALPHA if sig_mask[i] else INSIG_BAR_ALPHA
        left = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val <= 0:
                continue
            color = (categories[si][1] if si < len(categories)
                     else "#cccccc")
            ax.barh(y0, val, height=height, left=left,
                    color=color, alpha=alpha,
                    align="edge", linewidth=0, zorder=2)
            left += val

    # ── KL axis cosmetics ─────────────────────────────────────────────────
    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")

    # Genomic tick bars on Y axis (matching draw_epilogos_vertical)
    _add_vertical_tick_bars(ax, region_start, region_end)


# ---------------------------------------------------------------------------
# Private helpers (duplicated to keep this file self-contained)
# ---------------------------------------------------------------------------

def _gax(pos, region_start, region_end):
    """Convert genomic coordinate to [0, 1] axis fraction."""
    return (pos - region_start) / (region_end - region_start)


def _add_vertical_tick_bars(ax, region_start, region_end):
    """Draw unlabeled genomic tick bars on a vertical track's Y axis."""
    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        if 4 <= span / interval <= 12:
            break
    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)
    ticks_ax = [_gax(t, region_start, region_end) for t in ticks_genomic]
    ax.set_yticks(ticks_ax)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", left=True, right=False,
                   length=2.5, width=0.6, colors="#777777")


# ---------------------------------------------------------------------------
# Standalone demo (run directly to see what the overlay looks like)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(42)
    n = 200
    positions = [(i * 200, (i + 1) * 200) for i in range(n)]
    region_start, region_end = 0, n * 200

    # Simulated scores: 5 tracks
    score_mat = np.zeros((n, 5), dtype=np.float32)
    score_mat[40:60, 0] = rng.uniform(0.8, 2.0, 20)  # ATAC peak
    score_mat[80:100, 2] = rng.uniform(1.2, 2.5, 20) # H3K27ac peak
    score_mat[130:140, 1] = rng.uniform(1.5, 3.0, 10) # CTCF peak
    score_mat += rng.uniform(0, 0.1, score_mat.shape)  # background noise

    # Simulated pval: peaks match KL peaks
    pval_vals = np.zeros(n, dtype=np.float64)
    pval_vals[40:60] = rng.uniform(2, 8, 20)
    pval_vals[80:100] = rng.uniform(3, 10, 20)
    pval_vals[130:140] = rng.uniform(5, 12, 10)
    pval_vals += rng.uniform(0, 0.5, n)
    pval_positions = np.array([p[0] for p in positions], dtype=np.float64)

    categories = [
        ("ATAC",    "#00b050"),
        ("CTCF",    "#ff2200"),
        ("H3K27ac", "#00c864"),
        ("RNAseq+", "#6495ed"),
        ("RNAseq-", "#1a3a8f"),
    ]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 3), dpi=120)
    fig.patch.set_facecolor("#f4f3ef")

    # Top: with pval overlay
    draw_epilogos_with_pval_horizontal(
        ax1, positions, score_mat, 5,
        region_start, region_end, categories,
        pval_positions=pval_positions,
        pval_values=pval_vals,
        pval_alpha=0.05,
        label="With p-value overlay",
    )

    # Bottom: without (all bins full opacity)
    draw_epilogos_with_pval_horizontal(
        ax2, positions, score_mat, 5,
        region_start, region_end, categories,
        pval_positions=None,
        label="No p-value (all bins full opacity)",
    )

    fig.suptitle("draw_epilogos_with_pval_horizontal demo",
                 fontsize=9, y=1.01)
    fig.tight_layout()
    out = "/tmp/pval_overlay_demo.png"
    fig.savefig(out, dpi=120, bbox_inches="tight", facecolor="#f4f3ef")
    plt.close(fig)
    print(f"Demo saved: {out}")
