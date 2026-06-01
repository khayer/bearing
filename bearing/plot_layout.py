from __future__ import annotations

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt


def create_main_figure(use_pval_overlay: bool, has_loops: bool = True):
    hic_sz = 5.5
    epi_w = 1.0
    pval_w = 0.08 if use_pval_overlay else 0.55
    loop_w = 0.55
    diffp_w = 0.65
    axis_w = 0.45
    epi_h = 1.0
    pval_h = 0.08 if use_pval_overlay else 0.50
    loop_h = 0.55 if has_loops else 0.12
    diff_h = 0.8
    diffp_h = 0.60
    gene_h = 0.30
    axis_h = 0.28

    col_widths = [hic_sz, epi_w, pval_w, loop_w, epi_w, diffp_w, epi_w, axis_w]
    row_heights = [hic_sz, epi_h, pval_h, loop_h, diff_h, diffp_h, gene_h, axis_h]

    fig_w = sum(col_widths) + 0.6
    fig_h = sum(row_heights) + 0.6

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
    gs = gridspec.GridSpec(
        8, 8,
        figure=fig,
        width_ratios=col_widths,
        height_ratios=row_heights,
        hspace=0.04,
        wspace=0.01,
        left=0.07, right=0.97, top=0.97, bottom=0.04,
    )

    axes = {
        "hic": fig.add_subplot(gs[0, 0]),
        "epi_a": fig.add_subplot(gs[1, 0]),
        "pval_a": fig.add_subplot(gs[2, 0]),
        "loop_a": fig.add_subplot(gs[3, 0]),
        "diff_a": fig.add_subplot(gs[4, 0]),
        "diff_pval_a": fig.add_subplot(gs[5, 0]),
        "gene_h": fig.add_subplot(gs[6, 0]),
        "axis_h": fig.add_subplot(gs[7, 0]),
        "epi_b": fig.add_subplot(gs[0, 1]),
        "pval_b": fig.add_subplot(gs[0, 2]),
        "loop_b": fig.add_subplot(gs[0, 3]),
        "diff_b": fig.add_subplot(gs[0, 4]),
        "diff_pval_b": fig.add_subplot(gs[0, 5]),
        "gene_v": fig.add_subplot(gs[0, 6]),
        "axis_v": fig.add_subplot(gs[0, 7]),
        "legend": fig.add_subplot(gs[4:8, 5:8]),
    }
    return fig, axes
