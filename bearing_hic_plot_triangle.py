#!/usr/bin/env python3
"""
bearing_hic_plot_triangle.py
============================
Triangle-oriented Hi-C comparison figure for two conditions.

This module keeps the same input surface as bearing_hic_plot.py, but it uses
an alternate layout:
- an upward-pointing rotated Hi-C triangle at the top
- tracks in order: qcat A, qcat B, diff qcat B-A, diff p-value B-A, genes, coordinates
- a downward-pointing inverted triangle at the bottom

The comparison direction is always current condition minus the reference
condition supplied by batch_bearing_hic_plots.py.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from bearing_hic_plot import (
    ALL_CATEGORIES,
    HIC_CMAP,
    draw_diff_horizontal,
    draw_epilogos_horizontal,
    draw_loops_horizontal,
    draw_gene_track,
    draw_genomic_axis,
    draw_legend,
    draw_pval_diff_horizontal,
    load_categories_yaml,
    load_contact_matrix,
    load_genes,
    load_genes_gtf,
    load_highlights,
    load_loops,
    load_qcat_scores,
    load_regions_file,
    load_pval_track_values,
    load_stats_tsv_values,
    load_stats_tsv_with_categories,
    make_rgb_hic,
    parse_region,
    _resolve_contacts_for_region_resolution,
)
from bearing.plot_loaders import load_bed_for_region
from bearing.plot_tracks import draw_bed_track
from bearing_hic_plot_pval_overlay import draw_epilogos_with_pval_horizontal
from bearing.plot_tracks import draw_pval_manhattan_horizontal


TRIANGLE_TOP_VERTICES = np.array(
    [(0.50, 1.00), (0.00, 0.00), (1.00, 0.00)],
    dtype=float,
)
TRIANGLE_BOTTOM_VERTICES = np.array(
    [(0.00, 1.00), (1.00, 1.00), (0.50, 0.00)],
    dtype=float,
)


def _safe_filename(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in name)


def _draw_scalar_triangle(ax, matrix, inverted=False, hic_vmax=None):
    """Draw one matrix half as a true triangle without geometric warping."""
    arr = np.asarray(matrix, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[1] == 0:
        return

    n = min(arr.shape[0], arr.shape[1])
    arr = arr[:n, :n]

    # Upper-triangle (including diagonal) maps to an isosceles triangle:
    # x = (i + j) / (2n), y = (j - i) / n
    i_edges, j_edges = np.meshgrid(
        np.arange(n + 1, dtype=np.float64),
        np.arange(n + 1, dtype=np.float64),
        indexing="ij",
    )
    x = (i_edges + j_edges) / (2.0 * n)
    y = (j_edges - i_edges) / n
    if inverted:
        y = 1.0 - y

    c = np.ma.array(arr, mask=np.tril(np.ones((n, n), dtype=bool), k=-1))
    ax.pcolormesh(
        x,
        y,
        c,
        cmap=HIC_CMAP,
        vmin=0.0,
        vmax=hic_vmax,
        shading="flat",
        antialiased=False,
        zorder=1,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)


def _draw_rgb_triangle(ax, image, inverted=False):
    """Draw RGB Hi-C as an unwarped triangle via per-bin quads."""
    arr = np.asarray(image)
    if arr.ndim != 3 or arr.shape[2] < 3 or arr.shape[0] == 0 or arr.shape[1] == 0:
        return

    n = min(arr.shape[0], arr.shape[1])
    arr = arr[:n, :n, :3].astype(np.float64) / 255.0

    i_edges, j_edges = np.meshgrid(
        np.arange(n + 1, dtype=np.float64),
        np.arange(n + 1, dtype=np.float64),
        indexing="ij",
    )
    x = (i_edges + j_edges) / (2.0 * n)
    y = (j_edges - i_edges) / n
    if inverted:
        y = 1.0 - y

    verts = []
    colors = []
    for i in range(n):
        for j in range(i, n):
            verts.append([
                (x[i, j], y[i, j]),
                (x[i + 1, j], y[i + 1, j]),
                (x[i + 1, j + 1], y[i + 1, j + 1]),
                (x[i, j + 1], y[i, j + 1]),
            ])
            colors.append(arr[i, j])

    pc = PolyCollection(verts, facecolors=colors, edgecolors="none", antialiaseds=False, zorder=1)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)


def _palette_ab_colors(palette):
    """Return (color_a, color_b) hex matching the RGB Hi-C palette, so loop
    arcs/markers and condition labels use the same colors as the contact map."""
    return {
        "magenta-green": ("#cc00cc", "#00aa00"),
        "magenta-green-white": ("#cc00cc", "#00aa00"),
        "red-green": ("#cc0000", "#00aa00"),
        "blue-red": ("#1f5fcc", "#cc0000"),
        "green-blue": ("#00aa00", "#1f5fcc"),
    }.get(palette, ("#cc00cc", "#00aa00"))


def _overlay_loops_on_triangle(ax, loops, region_start, region_end, used_res,
                               matrix_n, inverted=False, color="#1f77b4",
                               marker="o", label=None):
    """Overlay loop apexes on a rotated Hi-C triangle.

    The triangle maps bin (i,j) to normalized axes x=(i+j)/(2n), y=(j-i)/n
    (see _draw_rgb_triangle). A loop connects genomic anchors p1<p2; convert to
    bin indices relative to region_start at the matrix resolution and place a
    marker at the apex of that bin pair. Markers are drawn filled with a white
    halo so they stay visible against the busy RGB contact background.
    """
    import matplotlib.patheffects as pe
    if not loops or matrix_n <= 0 or used_res <= 0:
        return
    halo = [pe.withStroke(linewidth=2.2, foreground="white")]
    drawn = 0
    for (s1, e1, s2, e2, _score) in loops:
        a = ((s1 + e1) // 2)
        b = ((s2 + e2) // 2)
        lo, hi = (a, b) if a <= b else (b, a)
        i = (lo - region_start) // used_res
        j = (hi - region_start) // used_res
        if i < 0 or j < 0 or i >= matrix_n or j >= matrix_n:
            continue
        # bin centers (the quad for (i,j) spans i..i+1, j..j+1)
        ic, jc = i + 0.5, j + 0.5
        x = (ic + jc) / (2.0 * matrix_n)
        y = (jc - ic) / float(matrix_n)
        if inverted:
            y = 1.0 - y
        ax.plot(x, y, marker=marker, markersize=7, markerfacecolor=color,
                markeredgecolor="white", markeredgewidth=1.0, alpha=0.95,
                zorder=6, path_effects=halo,
                label=label if drawn == 0 else None)
        drawn += 1
    if drawn:
        print(f"  overlaid {drawn} loop apex(es) ({label or 'loops'})")


def _draw_triangle_hic(
    ax,
    image,
    inverted: bool = False,
    rgb_mode: bool = False,
    hic_vmax=None,
):
    ax.set_axis_off()
    if rgb_mode:
        _draw_rgb_triangle(ax, image, inverted=inverted)
    else:
        _draw_scalar_triangle(ax, image, inverted=inverted, hic_vmax=hic_vmax)


def _triangle_figure_layout(num_beds: int = 0, bed_styles: list | None = None):
    fig = plt.figure(figsize=(12.0, 14.0), dpi=150)
    # Base rows before bed tracks: hic_top, loops_a, loops_b, qcat_a, qcat_b,
    # diff, pval_diff, pval_manhattan
    pre_bed = [4.3, 0.45, 0.45, 0.95, 0.95, 0.95, 0.75, 0.75]
    # Default gene/axis/bottom triangle heights
    gene_h = 0.40
    axis_h = 0.35
    bottom_h = 4.3

    bed_heights = []
    if num_beds and bed_styles:
        for i in range(num_beds):
            sty = bed_styles[i] if i < len(bed_styles) else "itemRgb"
            if sty == "cbe":
                bed_heights.append(0.12)
            else:
                bed_heights.append(0.28)
    else:
        bed_heights = [0.12] * num_beds

    height_ratios = pre_bed + bed_heights + [gene_h, axis_h, bottom_h]
    rows = len(height_ratios)

    gs = gridspec.GridSpec(
        rows,
        2,
        figure=fig,
        width_ratios=[7.3, 2.2],
        height_ratios=height_ratios,
        hspace=0.10,
        wspace=0.06,
        left=0.06,
        right=0.97,
        top=0.97,
        bottom=0.08,
    )

    axes = {}
    axes["hic_top"] = fig.add_subplot(gs[0, 0])
    axes["loops_a"] = fig.add_subplot(gs[1, 0])
    axes["loops_b"] = fig.add_subplot(gs[2, 0])
    axes["qcat_a"] = fig.add_subplot(gs[3, 0])
    axes["qcat_b"] = fig.add_subplot(gs[4, 0])
    axes["diff"] = fig.add_subplot(gs[5, 0])
    axes["pval_diff"] = fig.add_subplot(gs[6, 0])
    axes["pval_manhattan"] = fig.add_subplot(gs[7, 0])

    # bed rows start at index 8
    bed_start = 8
    for i in range(num_beds):
        axes[f"bed_{i}"] = fig.add_subplot(gs[bed_start + i, 0])

    gene_idx = bed_start + num_beds
    axes["gene"] = fig.add_subplot(gs[gene_idx, 0])
    axes["axis"] = fig.add_subplot(gs[gene_idx + 1, 0])
    axes["hic_bottom"] = fig.add_subplot(gs[gene_idx + 2, 0])

    # legend spans from qcat_a through gene (adjust for beds)
    legend_top = 1
    legend_bottom = gene_idx
    axes["legend"] = fig.add_subplot(gs[legend_top:legend_bottom, 1])
    return fig, axes


def make_figure_triangle(
    hic_a_path,
    hic_b_path,
    qcat_a_path,
    qcat_b_path,
    region_str,
    resolution,
    out_path,
    loops_path=None,
    genes_path=None,
    gtf_path=None,
    highlights_path=None,
    diff_qcat_path=None,
    pval_a_path=None,
    pval_b_path=None,
    pval_diff_path=None,
    loops_a_path=None,
    loops_b_path=None,
    pval_cutoff=0.05,
    label_a="Condition A",
    label_b="Condition B",
    categories=None,
    hic_fmt=None,
    pval_overlay=False,
    rgb_hic=False,
    rgb_palette="magenta-green",
    loops_on_triangle=True,
    beds=None,
    bed_style_overrides=None,
):
    chrom, region_start, region_end = parse_region(region_str)
    span = region_end - region_start
    print(f"\nRegion: {chrom}:{region_start:,}-{region_end:,}  ({span/1e6:.2f} Mb)")

    if categories is None:
        categories = ALL_CATEGORIES

    print("Loading Hi-C matrices...")
    mat_a, used_res_a = load_contact_matrix(
        hic_a_path, chrom, region_start, region_end,
        resolution, fmt=hic_fmt,
    )
    mat_b, used_res_b = load_contact_matrix(
        hic_b_path, chrom, region_start, region_end,
        resolution, fmt=hic_fmt,
    )
    matrix_n = min(mat_a.shape[0], mat_b.shape[0])
    print(f"  Matrix size: {matrix_n} x {matrix_n} bins")
    if used_res_a == used_res_b:
        print(f"  Hi-C resolution used: {used_res_a:,} bp")
        resolution_note = f"Resolution: {used_res_a:,} bp"
    else:
        print(
            "  WARNING: Hi-C resolutions differ between A and B: "
            f"A={used_res_a:,} bp, B={used_res_b:,} bp"
        )
        resolution_note = f"Resolution: A={used_res_a:,} bp, B={used_res_b:,} bp"

    print("Loading Bearing epilogos scores...")
    pos_a, scores_a, ns_a = load_qcat_scores(
        qcat_a_path, chrom, region_start, region_end)
    pos_b, scores_b, ns_b = load_qcat_scores(
        qcat_b_path, chrom, region_start, region_end)
    num_states = max(ns_a, ns_b)
    if scores_a.shape[1] < num_states:
        scores_a = np.pad(scores_a, ((0, 0), (0, num_states - scores_a.shape[1])))
    if scores_b.shape[1] < num_states:
        scores_b = np.pad(scores_b, ((0, 0), (0, num_states - scores_b.shape[1])))
    print(f"  Cond A: {len(pos_a)} bins, Cond B: {len(pos_b)} bins, {num_states} states")

    loops = None
    if loops_path:
        print("Loading loops...")
        loops = load_loops(loops_path, chrom, region_start, region_end)
        print(f"  {len(loops)} loops in region")

    loops_a = None
    if loops_a_path:
        print("Loading condition A loops...")
        loops_a = load_loops(loops_a_path, chrom, region_start, region_end)
        print(f"  {len(loops_a)} loops in region")

    loops_b = None
    if loops_b_path:
        print("Loading condition B loops...")
        loops_b = load_loops(loops_b_path, chrom, region_start, region_end)
        print(f"  {len(loops_b)} loops in region")

    genes = None
    if genes_path:
        print("Loading genes...")
        genes = load_genes(genes_path, chrom, region_start, region_end)
        print(f"  {len(genes)} gene records in region")
    elif gtf_path:
        print("Loading genes from GTF...")
        genes = load_genes_gtf(gtf_path, chrom, region_start, region_end)
        print(f"  {len(genes)} gene records in region")

    highlights = None
    if highlights_path:
        print("Loading highlights...")
        highlights = load_highlights(highlights_path, chrom, region_start, region_end)
        print(f"  {len(highlights)} highlight regions")

    pos_diff, scores_diff, ns_diff = [], np.zeros((0, 1), dtype=np.float32), 1
    has_diff = diff_qcat_path is not None
    if has_diff:
        print("Loading differential epilogos scores...")
        pos_diff, scores_diff, ns_diff = load_qcat_scores(
            diff_qcat_path, chrom, region_start, region_end)
        # Keep triangle mode consistent with batch comparison direction:
        # render diff as current-reference (label_b - label_a).
        scores_diff = -scores_diff
        ns_diff_full = max(ns_diff, num_states)
        if scores_diff.shape[1] < ns_diff_full:
            scores_diff = np.pad(
                scores_diff,
                ((0, 0), (0, ns_diff_full - scores_diff.shape[1])),
            )
        print(f"  Diff: {len(pos_diff)} bins")

    pos_pval_diff, vals_pval_diff = None, None
    pval_diff_kl_scores = None
    pval_diff_cat_names = None
    if pval_diff_path is not None:
        print("Loading signed differential p-value track...")
        # Check if it's a TSV to load category KL scores
        if str(pval_diff_path).endswith(".tsv"):
            pos_pval_diff, vals_pval_diff, pval_diff_kl_scores, pval_diff_cat_names = load_stats_tsv_with_categories(
                pval_diff_path, chrom, region_start, region_end)
        else:
            pos_pval_diff, vals_pval_diff = load_pval_track_values(
                pval_diff_path, chrom, region_start, region_end)
        # Keep p-value diff direction consistent with B-A display labels.
        vals_pval_diff = -vals_pval_diff
        if pval_diff_kl_scores is not None and pval_diff_kl_scores.size > 0:
            # Also negate the KL scores to match the direction
            pval_diff_kl_scores = -pval_diff_kl_scores
        print(f"  {len(pos_pval_diff)} bins")

    pos_pval_a, vals_pval_a = None, None
    pos_pval_b, vals_pval_b = None, None
    if pval_overlay:
        if pval_a_path is not None:
            print("Loading condition A p-value track for overlay...")
            pos_pval_a, vals_pval_a = load_pval_track_values(
                pval_a_path, chrom, region_start, region_end)
            print(f"  {len(pos_pval_a)} bins")
        if pval_b_path is not None:
            print("Loading condition B p-value track for overlay...")
            pos_pval_b, vals_pval_b = load_pval_track_values(
                pval_b_path, chrom, region_start, region_end)
            print(f"  {len(pos_pval_b)} bins")

    pval_cutoff_value = None
    if pval_cutoff is not None and pval_cutoff > 0:
        pval_cutoff_value = -math.log10(pval_cutoff)

    has_genes = genes is not None and len(genes) > 0
    has_pval_diff = pos_pval_diff is not None and len(pos_pval_diff) > 0
    # Manhattan row should display TSV-based differential p-values (diff stats TSV)
    has_pval_manhattan = pos_pval_diff is not None and len(pos_pval_diff) > 0
    # Load BED overlays (if any) and decide per-file styles
    beds = beds or []
    bed_style_overrides = bed_style_overrides or {}
    bed_features_list = []
    bed_styles = []
    for bed_path in beds:
        feats = load_bed_for_region(bed_path, chrom, region_start, region_end)
        # determine style: override by exact filename match or basename
        style = None
        key1 = str(bed_path)
        key2 = Path(bed_path).name
        if key1 in bed_style_overrides:
            style = bed_style_overrides[key1]
        elif key2 in bed_style_overrides:
            style = bed_style_overrides[key2]
        else:
            # auto-detect: if any feature has item_rgb, prefer itemRgb
            if any(f.get("item_rgb") is not None for f in feats) and len(feats) > 0:
                style = "itemRgb"
            else:
                style = "cbe"
        bed_features_list.append(feats)
        bed_styles.append(style)

    fig, axes = _triangle_figure_layout(num_beds=len(beds), bed_styles=bed_styles)
    ax_hic_top = axes["hic_top"]
    ax_loops_a = axes["loops_a"]
    ax_loops_b = axes["loops_b"]
    ax_qcat_a = axes["qcat_a"]
    ax_qcat_b = axes["qcat_b"]
    ax_diff = axes["diff"]
    ax_pval_diff = axes["pval_diff"]
    ax_pval_manhattan = axes["pval_manhattan"]
    ax_gene = axes["gene"]
    ax_axis = axes["axis"]
    ax_hic_bottom = axes["hic_bottom"]
    ax_legend = axes["legend"]

    # Colors matched to the RGB Hi-C palette so labels, arc tracks, and triangle
    # markers all agree (e.g. red-green -> A=red, B=green).
    color_a, color_b = _palette_ab_colors(rgb_palette) if rgb_hic else ("#cc00cc", "#00a83a")

    # draw bed tracks if present
    for i, bed_path in enumerate(beds or []):
        ax_bed = axes.get(f"bed_{i}")
        if ax_bed is None:
            continue
        style = bed_styles[i]
        label = Path(bed_path).stem
        draw_bed_track(
            ax_bed,
            bed_features_list[i],
            region_start,
            region_end,
            style=style,
            label=label,
        )

    print("Rendering top triangle panel...")
    if rgb_hic:
        top_image = make_rgb_hic(mat_a, mat_b, palette=rgb_palette)
        bottom_image = make_rgb_hic(mat_b, mat_a, palette=rgb_palette)
    else:
        top_image = mat_a
        bottom_image = mat_b

    hic_vmax = None
    if not rgb_hic:
        finite_vals = np.concatenate([
            np.asarray(mat_a, dtype=np.float64).ravel(),
            np.asarray(mat_b, dtype=np.float64).ravel(),
        ])
        finite_vals = finite_vals[np.isfinite(finite_vals)]
        if finite_vals.size > 0:
            hic_vmax = float(np.percentile(finite_vals, 97))
            if hic_vmax <= 0:
                hic_vmax = float(np.max(finite_vals)) if finite_vals.size else 1.0

    _draw_triangle_hic(
        ax_hic_top,
        top_image,
        inverted=False,
        rgb_mode=rgb_hic,
        hic_vmax=hic_vmax,
    )
    # Optional loop apex markers on the top triangle (toggle via loops_on_triangle).
    # Arc tracks below always show the loops; the triangle markers are additive.
    if loops_on_triangle:
        _overlay_loops_on_triangle(
            ax_hic_top, loops_a, region_start, region_end, used_res_a, matrix_n,
            inverted=False, color=color_a, marker="o", label="%s loops" % label_a)
        if rgb_hic:
            _overlay_loops_on_triangle(
                ax_hic_top, loops_b, region_start, region_end, used_res_b, matrix_n,
                inverted=False, color=color_b, marker="s", label="%s loops" % label_b)
        if loops is not None:
            _overlay_loops_on_triangle(
                ax_hic_top, loops, region_start, region_end, used_res_a, matrix_n,
                inverted=False, color="#222222", marker="o", label="loops")

    # Two arc tracks (condition A, condition B) below the triangle, colored to
    # match the Hi-C palette.
    draw_loops_horizontal(
        ax_loops_a, loops_a or [], region_start, region_end,
        label="%s loops" % label_a, color=color_a, anchor_color=color_a)
    draw_loops_horizontal(
        ax_loops_b, loops_b or [], region_start, region_end,
        label="%s loops" % label_b, color=color_b, anchor_color=color_b)
    if rgb_hic:
        ax_hic_top.text(
            0.06,
            0.94,
            label_a,
            transform=ax_hic_top.transAxes,
            fontsize=8,
            color=color_a,
            ha="left",
            va="top",
            fontweight="bold",
        )
        ax_hic_top.text(
            0.94,
            0.94,
            label_b,
            transform=ax_hic_top.transAxes,
            fontsize=8,
            color=color_b,
            ha="right",
            va="top",
            fontweight="bold",
        )
    else:
        ax_hic_top.text(
            0.50,
            0.94,
            label_a,
            transform=ax_hic_top.transAxes,
            fontsize=8,
            color="#222222",
            ha="center",
            va="top",
            fontweight="bold",
        )
    ax_hic_top.text(
        0.02,
        0.99,
        resolution_note,
        transform=ax_hic_top.transAxes,
        fontsize=6,
        ha="left",
        va="top",
        color="#222222",
        zorder=7,
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
    )

    y_max_shared = 0.0
    if scores_a.size:
        y_max_shared = max(y_max_shared, float(scores_a.max()))
    if scores_b.size:
        y_max_shared = max(y_max_shared, float(scores_b.max()))
    if has_diff and scores_diff.size:
        y_max_shared = max(y_max_shared, float(np.abs(scores_diff).max()))
    y_max_shared = 1.0 if y_max_shared == 0 else y_max_shared * 1.05

    y_max_shared_pval_diff = 0.0
    if has_pval_diff and vals_pval_diff.size:
        y_max_shared_pval_diff = float(np.abs(vals_pval_diff).max())
    y_max_shared_pval_diff = 1.0 if y_max_shared_pval_diff == 0 else y_max_shared_pval_diff * 1.05

    print("Rendering horizontal tracks...")
    if pval_overlay and pos_pval_a is not None and vals_pval_a is not None:
        draw_epilogos_with_pval_horizontal(
            ax_qcat_a, pos_a, scores_a, num_states,
            region_start, region_end,
            categories[:num_states],
            pval_positions=pos_pval_a,
            pval_values=vals_pval_a,
            pval_alpha=pval_cutoff,
            highlights=highlights,
            label=label_a,
            y_max=y_max_shared,
        )
    else:
        draw_epilogos_horizontal(
            ax_qcat_a, pos_a, scores_a, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label=label_a,
            y_max=y_max_shared,
        )

    if pval_overlay and pos_pval_b is not None and vals_pval_b is not None:
        draw_epilogos_with_pval_horizontal(
            ax_qcat_b, pos_b, scores_b, num_states,
            region_start, region_end,
            categories[:num_states],
            pval_positions=pos_pval_b,
            pval_values=vals_pval_b,
            pval_alpha=pval_cutoff,
            highlights=highlights,
            label=label_b,
            y_max=y_max_shared,
        )
    else:
        draw_epilogos_horizontal(
            ax_qcat_b, pos_b, scores_b, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label=label_b,
            y_max=y_max_shared,
        )

    if has_diff:
        draw_diff_horizontal(
            ax_diff, pos_diff, scores_diff, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label="diff qcat " + label_b + " - " + label_a,
            diff_max=y_max_shared,
        )
    else:
        ax_diff.set_axis_off()

    if has_pval_diff:
        draw_pval_diff_horizontal(
            ax_pval_diff, pos_pval_diff, vals_pval_diff,
            region_start, region_end,
            highlights=highlights,
            label="diff p-value " + label_b + " - " + label_a,
            y_max=y_max_shared_pval_diff,
            cutoff_value=pval_cutoff_value,
            diff_score_positions=pos_diff if has_diff else None,
            diff_score_matrix=scores_diff if has_diff else None,
            categories=categories[:num_states] if has_diff else None,
        )
    else:
        ax_pval_diff.set_axis_off()

    if has_pval_manhattan:
        # Use the TSV-based differential p-values (pos_pval_diff / vals_pval_diff)
        # Color Manhattan dots by dominant category from the TSV KL scores if available,
        # otherwise fall back to matching against qcat A scores.
        y_max_pval_manhattan = 0.0
        if vals_pval_diff is not None and vals_pval_diff.size > 0:
            y_max_pval_manhattan = float(np.abs(vals_pval_diff).max())
        y_max_pval_manhattan = 1.0 if y_max_pval_manhattan == 0 else y_max_pval_manhattan * 1.05

        # Color Manhattan by the dominant category in condition A (reference)
        draw_pval_manhattan_horizontal(
            ax_pval_manhattan,
            pos_pval_diff if pos_pval_diff is not None else [],
            vals_pval_diff if vals_pval_diff is not None else [],
            region_start, region_end,
            highlights=highlights,
            label="diff p-value (TSV) " + label_b + " - " + label_a,
            score_positions=pos_a,
            score_matrix=scores_a,
            categories=categories[:num_states],
            y_max=y_max_pval_manhattan,
            cutoff_value=pval_cutoff_value,
            kl_scores=pval_diff_kl_scores,
        )
    else:
        ax_pval_manhattan.set_axis_off()

    if has_genes:
        draw_gene_track(
            ax_gene, genes,
            region_start, region_end,
            highlights=highlights,
            label="Genes",
        )
    else:
        ax_gene.set_axis_off()

    draw_genomic_axis(ax_axis, region_start, region_end, chrom)

    print("Rendering bottom triangle panel...")
    # When RGB Hi-C mode is on, we don't render the lower Hi-C triangle
    # to save rendering time and because the RGB image is already shown above.
    if rgb_hic:
        ax_hic_bottom.set_axis_off()
    else:
        _draw_triangle_hic(
            ax_hic_bottom,
            bottom_image,
            inverted=True,
            rgb_mode=rgb_hic,
            hic_vmax=hic_vmax,
        )
        if loops_on_triangle:
            _overlay_loops_on_triangle(
                ax_hic_bottom, loops_b, region_start, region_end, used_res_b,
                matrix_n, inverted=True, color=color_b, marker="s",
                label="%s loops" % label_b)
        ax_hic_bottom.text(
            0.50,
            0.06,
            label_b,
            transform=ax_hic_bottom.transAxes,
            fontsize=8,
            color="#222222",
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    print("Rendering legend...")
    draw_legend(
        ax_legend,
        categories[:num_states],
        num_states,
        hic_vmax=hic_vmax,
        rgb_mode=rgb_hic,
        label_a=label_a,
        label_b=label_b,
        rgb_palette=rgb_palette,
    )

    print(f"Saving to {out_path}...")
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Done.  Figure: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Composite Hi-C + Bearing epilogos triangle figure "
            "for two conditions at a genomic locus."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--contact-a", required=True, metavar="FILE")
    parser.add_argument("--contact-b", required=True, metavar="FILE")
    parser.add_argument("--format", dest="hic_fmt", choices=["cool", "hic"], default=None, metavar="FMT")
    parser.add_argument("--qcat-a", required=True, metavar="FILE")
    parser.add_argument("--qcat-b", required=True, metavar="FILE")
    parser.add_argument("--diff-qcat", metavar="FILE", default=None)
    parser.add_argument("--pval-a", metavar="FILE", default=None)
    parser.add_argument("--pval-b", metavar="FILE", default=None)
    parser.add_argument("--pval-diff", metavar="FILE", default=None)
    parser.add_argument("--pval-cutoff", type=float, default=0.05, metavar="P")
    parser.add_argument("--pval-overlay", action="store_true")
    parser.add_argument("--rgb-hic", action="store_true")
    parser.add_argument(
        "--rgb-palette",
        choices=["magenta-green", "red-green", "blue-red", "green-blue", "magenta-green-white"],
        default="magenta-green",
    )
    parser.add_argument(
        "--no-loops-on-triangle", dest="loops_on_triangle", action="store_false",
        help="Show loops only as arc tracks, not as apex markers on the triangle.")
    parser.set_defaults(loops_on_triangle=True)
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument("--region", metavar="CHR:START-END")
    region_group.add_argument("--regions-file", metavar="TSV")
    parser.add_argument("--resolution", type=int, default=10000, metavar="BP")
    parser.add_argument("--out", metavar="FILE")
    parser.add_argument("--outdir", metavar="DIR", default=".")
    parser.add_argument("--label-a", default="Condition A", metavar="STR")
    parser.add_argument("--label-b", default="Condition B", metavar="STR")
    parser.add_argument("--loops", metavar="BEDPE")
    parser.add_argument("--loops-a", metavar="BEDPE", default=None)
    parser.add_argument("--loops-b", metavar="BEDPE", default=None)
    parser.add_argument("--genes", metavar="BED")
    parser.add_argument("--gtf", metavar="GTF")
    parser.add_argument("--highlights", metavar="BED")
    parser.add_argument(
        "--categories", metavar="YAML",
        help="YAML file defining category names and colors.",
    )
    parser.add_argument(
        "--bed",
        action="append",
        default=[],
        metavar="FILE",
        help=(
            "BED file to overlay as a track row. Repeat for multiple files. "
            "Files are stacked above the gene track in the order given."
        ),
    )
    parser.add_argument(
        "--bed-style",
        action="append",
        default=[],
        metavar="FILE=STYLE",
        help=(
            "Per-file rendering style override. STYLE is 'cbe' or 'itemRgb'. "
            "Default: itemRgb if BED9 with color column, otherwise cbe. Repeat for each file you want to override."
        ),
    )

    args = parser.parse_args()

    if not (0 < args.pval_cutoff <= 1):
        parser.error("--pval-cutoff must be > 0 and <= 1")
    if args.genes and args.gtf:
        parser.error("Use only one of --genes or --gtf")
    if args.region and not args.out:
        parser.error("--out is required when using --region")

    cli_categories = None
    if args.categories:
        cli_categories, _ = load_categories_yaml(args.categories)
        print(f"Categories loaded from: {args.categories}  ({len(cli_categories)} states)")

    print(f"\n{'=' * 60}")
    print("  bearing_hic_plot_triangle.py")
    print(f"  Format:     {args.hic_fmt or 'auto-detect'}")
    print(f"  Cond A:     {args.label_a}  ({args.contact_a})")
    print(f"  Cond B:     {args.label_b}  ({args.contact_b})")
    if args.loops:
        print(f"  Loops:      {args.loops}")
    if args.loops_a:
        print(f"  Loops A:    {args.loops_a}")
    if args.loops_b:
        print(f"  Loops B:    {args.loops_b}")
    if args.genes:
        print(f"  Genes:      {args.genes}")
    if args.gtf:
        print(f"  GTF:        {args.gtf}")
    if args.highlights:
        print(f"  Highlights: {args.highlights}")
    if args.diff_qcat:
        print(f"  Diff qcat:  {args.diff_qcat}")
    if args.pval_a:
        print(f"  P-val A:    {args.pval_a}")
    if args.pval_b:
        print(f"  P-val B:    {args.pval_b}")
    if args.pval_diff:
        print(f"  P-val diff: {args.pval_diff}")
    if args.pval_a or args.pval_b or args.pval_diff:
        print(f"  P cutoff:   {args.pval_cutoff} (-log10={-math.log10(args.pval_cutoff):.3f})")
    if args.pval_overlay:
        print("  P overlay:  enabled")
    if args.rgb_hic:
        print("  Hi-C mode:  RGB")
        print(f"  RGB palette: {args.rgb_palette}")
    print(f"{'=' * 60}\n")

    if args.regions_file:
        regions = load_regions_file(args.regions_file)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"Batch mode: processing {len(regions)} region(s) from {args.regions_file}")
        print(f"Output directory: {outdir}\n")

        # parse bed-style overrides into mapping
        bed_style_overrides = {}
        for s in args.bed_style:
            if "=" not in s:
                continue
            left, right = s.split("=", 1)
            style = ("itemRgb" if right.lower() == "itemrgb" else "cbe")
            bed_style_overrides[left] = style

        for i, reg in enumerate(regions, start=1):
            name = reg["name"]
            region_str = reg["region"]
            resolution = reg["resolution"] if reg["resolution"] else args.resolution
            label = reg["label"] if reg["label"] else name

            hic_a_path, hic_b_path, resolve_note = _resolve_contacts_for_region_resolution(
                args.contact_a,
                args.contact_b,
                resolution,
                hic_fmt=args.hic_fmt,
            )

            if reg["out"]:
                out_path = outdir / reg["out"]
            else:
                out_path = outdir / f"{_safe_filename(name)}_triangle.pdf"

            print(f"[{i}/{len(regions)}] {name}: {region_str} @ {resolution:,} bp -> {out_path}")
            print(f"  Hi-C files: {hic_a_path} | {hic_b_path}")
            print(f"  Resolution match: {resolve_note}")

            try:
                make_figure_triangle(
                    hic_a_path=hic_a_path,
                    hic_b_path=hic_b_path,
                    qcat_a_path=args.qcat_a,
                    qcat_b_path=args.qcat_b,
                    region_str=region_str,
                    resolution=resolution,
                    out_path=out_path,
                    loops_path=args.loops,
                    loops_a_path=args.loops_a,
                    loops_b_path=args.loops_b,
                    genes_path=args.genes,
                    gtf_path=args.gtf,
                    highlights_path=args.highlights,
                    diff_qcat_path=args.diff_qcat,
                    pval_a_path=args.pval_a,
                    pval_b_path=args.pval_b,
                    pval_diff_path=args.pval_diff,
                    pval_cutoff=args.pval_cutoff,
                    label_a=args.label_a,
                    label_b=args.label_b,
                    hic_fmt=args.hic_fmt,
                    categories=cli_categories,
                    rgb_hic=args.rgb_hic,
                    rgb_palette=args.rgb_palette,
                    loops_on_triangle=args.loops_on_triangle,
                    beds=args.bed,
                    bed_style_overrides=bed_style_overrides,
                    pval_overlay=args.pval_overlay,
                )
            except Exception as e:
                print(f"  ERROR: Failed to process region {name}: {e}")
                continue

        print(f"\nBatch complete. {len(regions)} figure(s) written to {outdir}/")
    else:
        print("Single region mode:")
        print(f"  Region:     {args.region}")
        print(f"  Resolution: {args.resolution:,} bp")
        print(f"  Output:     {args.out}\n")

        make_figure_triangle(
            hic_a_path=args.contact_a,
            hic_b_path=args.contact_b,
            qcat_a_path=args.qcat_a,
            qcat_b_path=args.qcat_b,
            region_str=args.region,
            resolution=args.resolution,
            out_path=args.out,
            loops_path=args.loops,
            loops_a_path=args.loops_a,
            loops_b_path=args.loops_b,
            genes_path=args.genes,
            gtf_path=args.gtf,
            highlights_path=args.highlights,
            diff_qcat_path=args.diff_qcat,
            pval_a_path=args.pval_a,
            pval_b_path=args.pval_b,
            pval_diff_path=args.pval_diff,
            pval_cutoff=args.pval_cutoff,
            pval_overlay=args.pval_overlay,
            label_a=args.label_a,
            label_b=args.label_b,
            hic_fmt=args.hic_fmt,
            categories=cli_categories,
            rgb_hic=args.rgb_hic,
            rgb_palette=args.rgb_palette,
            loops_on_triangle=args.loops_on_triangle,
            beds=args.bed,
            bed_style_overrides=bed_style_overrides,
        )


if __name__ == "__main__":
    main()
