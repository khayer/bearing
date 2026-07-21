#!/usr/bin/env python3
"""
bearing_hic_combined_plot.py
============================
Combined two-condition BEARING + Hi-C browser figure.

Merges the two figure styles you have been using:

  From bearing_hic_plot_triangle.py (rendered with the SAME shared helpers,
  so these look identical to what you already produce):
    - per-condition qcat epilogos tracks (qcat A, qcat B)
    - diff qcat track (label_a - label_b)
    - Manhattan "lollipop" differential p-value track
    - gene track, BED overlay rows, category legend
    - optional combined RGB Hi-C triangle (--rgb-hic)

  From bearing_hic_kl_track_plot_v10.py (ported here as self-contained code):
    - single + diff Hi-C panels as genomic-coordinate rotated triangles
      (--show-hic A,B,diff), each with an inset colorbar
    - insulation A vs B and |delta insulation| with percentile significance
      markers, matched genome-wide (or per-chrom) via --threshold-scope
    - per-track KL decomposition: one bar panel per chromatin category, with
      two-tier percentile markers (default top 1% light / top 0.1% dark)

Layout, top to bottom:
    [combined RGB Hi-C triangle]   (optional, --rgb-hic)
    Hi-C A / Hi-C B / Hi-C diff    (genomic triangles, per --show-hic)
    [insulation A vs B]            (if --insul-A/--insul-B given)
    [|delta insulation| + sig]
    qcat A
    qcat B
    diff qcat (label_a - label_b)
    diff p-value (Manhattan lollipops)
    [diff p-value fill]            (optional, --pval-fill)
    BED rows ...
    genes
    per-track KL decomposition ...
    genomic coordinate axis

Drop-in for batch_bearing_hic_plots.py: accepts the same CLI surface as
bearing_hic_plot_triangle.py, plus the v10-style options (--insul-A/B,
--show-hic, --max-distance, --diff-stats, --track-pctile-*). The genomic
Hi-C triangles and insulation are read directly with cooler from the
--contact-a/--contact-b .cool files (.hic is not supported in this mode).

Direction convention: the qcat diff, Manhattan p-values, KL decomposition,
and Hi-C diff are all shown as label_a - label_b (reference minus condition,
e.g. DN - S3T3) so the whole figure reads in one direction (positive =
enriched in A, the reference). This matches the native sign of the diff-qcat
and stats files and the per-track decomposition.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm

try:
    import cooler
    HAVE_COOLER = True
except ImportError:
    HAVE_COOLER = False

from bearing_hic_plot import (
    ALL_CATEGORIES,
    HIC_CMAP,
    draw_diff_horizontal,
    draw_epilogos_horizontal,
    draw_gene_track,
    draw_genomic_axis,
    draw_legend,
    draw_pval_diff_horizontal,
    draw_loops_horizontal,
    load_categories_yaml,
    load_genes,
    load_genes_gtf,
    load_highlights,
    load_loops,
    load_qcat_scores,
    load_regions_file,
    load_pval_track_values,
    load_stats_tsv_with_categories,
    make_rgb_hic,
    parse_region,
)
from bearing.plot_loaders import load_bed_for_region
from bearing.plot_tracks import draw_bed_track
from bearing_hic_plot_pval_overlay import draw_epilogos_with_pval_horizontal
from bearing.plot_tracks import draw_pval_manhattan_horizontal


# Autosomes for genome-wide percentile thresholds (mouse chr1-19 and human
# chr1-22; extra names are harmless).
AUTOSOMES = {"chr{}".format(i) for i in range(1, 23)}

# Approximate width of column 0 (the data column) in inches, used to size the
# rotated Hi-C triangles at a correct genomic:distance aspect ratio.
COL0_WIDTH_IN = 8.3


def _safe_filename(name):
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in name)


# ---------------------------------------------------------------------
# Genomic-coordinate Hi-C triangles (ported from v10)
# ---------------------------------------------------------------------

def fetch_cool_region(path, chrom, start, end, balance=True):
    if not HAVE_COOLER:
        sys.exit("cooler not installed: pip install cooler")
    c = cooler.Cooler(str(path))
    binsize = c.binsize
    region = "{}:{}-{}".format(chrom, start, end)
    M = c.matrix(balance=balance, sparse=True).fetch(region).toarray()
    bins = c.bins().fetch(region).reset_index(drop=True)
    return M, bins["start"].to_numpy(), binsize


def plot_hic_triangle(ax, mat, bin_starts, binsize, region_start, region_end,
                      max_distance, cmap, norm, label=None):
    """Plot mat as a 45-degree rotated triangle (genomic axis horizontal,
    contact distance vertical, upper triangle only)."""
    n = mat.shape[0]
    if n == 0:
        return None
    max_d_bins = max(1, int(max_distance / binsize))
    hw = binsize / 2.0
    polys = []
    values = []
    for i in range(n):
        bi = bin_starts[i]
        jmax = min(n, i + max_d_bins + 1)
        for j in range(i, jmax):
            v = mat[i, j]
            if not np.isfinite(v) or v == 0:
                continue
            bj = bin_starts[j]
            cx = (bi + bj) / 2.0 + hw
            cy = (bj - bi) / 2.0
            polys.append(((cx - hw, cy), (cx, cy + hw),
                          (cx + hw, cy), (cx, cy - hw)))
            values.append(v)
    if not polys:
        return None
    pc = PolyCollection(polys, array=np.asarray(values),
                        cmap=cmap, norm=norm,
                        edgecolor="none", linewidth=0, antialiased=False)
    ax.add_collection(pc)
    ax.set_xlim(region_start, region_end)
    # Cap the vertical extent at the triangle's actual apex: contacts can't span
    # farther than the region itself, so for small regions use region_w/2 rather
    # than max_distance/2 (which would leave whitespace above the triangle).
    eff_max_d = min(max_distance, region_end - region_start)
    ax.set_ylim(0, eff_max_d / 2.0 + binsize)
    ax.set_yticks([])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if label:
        ax.set_ylabel(label, fontsize=8, rotation=0, ha="right", va="center")
    return pc


def _add_inset_colorbar(ax, pc, cbar_label):
    if pc is None:
        return
    # Use the native Axes.inset_axes (axes-fraction bounds) rather than
    # axes_grid1.inset_locator, whose anchored locator can drift at draw time
    # (frame box correct but the color gradient offset). [x0, y0, w, h] in
    # axes fraction places the colorbar reliably at the upper-right.
    cax = ax.inset_axes([0.80, 0.90, 0.18, 0.045])
    cb = ax.figure.colorbar(pc, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=5, length=2, pad=1)
    cb.set_label(cbar_label, fontsize=6, labelpad=1)
    vmin, vmax = pc.get_clim()
    cb.set_ticks([vmin, vmax])
    cb.ax.set_xticklabels(["{:.2g}".format(vmin), "{:.2g}".format(vmax)])


# ---------------------------------------------------------------------
# Combined RGB / scalar triangle (normalized coords, from triangle.py)
# ---------------------------------------------------------------------

def _draw_rgb_triangle(ax, image, inverted=False):
    arr = np.asarray(image)
    if arr.ndim != 3 or arr.shape[2] < 3 or arr.shape[0] == 0 or arr.shape[1] == 0:
        return
    n = min(arr.shape[0], arr.shape[1])
    arr = arr[:n, :n, :3].astype(np.float64) / 255.0
    i_edges, j_edges = np.meshgrid(
        np.arange(n + 1, dtype=np.float64),
        np.arange(n + 1, dtype=np.float64),
        indexing="ij")
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
    pc = PolyCollection(verts, facecolors=colors, edgecolors="none",
                        antialiaseds=False, zorder=1)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)


# ---------------------------------------------------------------------
# Insulation (ported from v10)
# ---------------------------------------------------------------------

def load_insul_bm(path, chrom, start, end):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit(str(path) + ": .bm needs >=4 columns")
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    out["score"] = df.iloc[:, 3:].mean(axis=1).values
    out["center"] = (out["start"] + out["end"]) / 2
    out = out[(out["chrom"] == chrom) & (out["end"] > start) & (out["start"] < end)]
    return out.reset_index(drop=True)


def load_insul_bm_chr(path, chrom):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit(str(path) + ": .bm needs >=4 columns")
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    out["score"] = df.iloc[:, 3:].mean(axis=1).values
    out["center"] = (out["start"] + out["end"]) / 2
    out = out[out["chrom"] == chrom].reset_index(drop=True)
    return out


def load_insul_bm_genome(path, autosomes_only=True):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit(str(path) + ": .bm needs >=4 columns")
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    out["score"] = df.iloc[:, 3:].mean(axis=1).values
    out["center"] = (out["start"] + out["end"]) / 2
    if autosomes_only:
        out = out[out["chrom"].isin(AUTOSOMES)]
    return out.reset_index(drop=True)


def compute_insul_pctile_threshold(ins_A, ins_B, pctile=95.0):
    """Per-bin |delta insul| significance via percentile threshold over the
    bins passed in (chrom-restricted or genome-wide)."""
    a = ins_A[["chrom", "center", "score"]].rename(columns={"score": "score_A"})
    b = ins_B[["chrom", "center", "score"]].rename(columns={"score": "score_B"})
    merged = pd.merge(a, b, on=["chrom", "center"]).sort_values(
        ["chrom", "center"]).reset_index(drop=True)
    obs = (merged["score_A"] - merged["score_B"]).abs().values
    finite = obs[np.isfinite(obs)]
    threshold = float(np.nanpercentile(finite, pctile)) if len(finite) else np.nan
    merged["obs"] = obs
    merged["threshold"] = threshold
    merged["pctile"] = pctile
    merged["sig"] = np.isfinite(obs) & (obs >= threshold)
    return merged


def write_insul_sig_table(df, out_path, resolution, threshold_scope, chrom,
                          insul_a_path, insul_b_path, pctile):
    """Export the per-bin |delta insulation| table that
    compute_insul_pctile_threshold() builds.

    Export-only: this reads the already-computed table and writes it to disk. It
    does not change the plot, the threshold, or any computation. `df` must be the
    FULL pre-restriction (genome-wide or chrom-wide) table, i.e. the same scope
    the threshold was computed over -- captured BEFORE the plotted-chromosome
    filter -- so the exported `sig` column and threshold are self-consistent.

    If out_path ends in .bedgraph/.bg, writes a 4-column bedgraph of the |delta|
    values (chrom, start, end, delta_abs); otherwise a full TSV. `start`/`end`
    are derived from `center` using `resolution`.
    """
    res = int(resolution)
    thr = float(df["threshold"].iloc[0]) if len(df) else float("nan")
    scope_desc = ("genome-wide (autosomes)" if threshold_scope == "genome"
                  else "chrom-restricted (%s)" % chrom)
    header = ("# bearing_hic_combined_plot.py |delta insulation| export | "
              "insul_A=%s insul_B=%s resolution=%d percentile=%g scope=%s "
              "threshold=%.6g" % (insul_a_path, insul_b_path, res, pctile,
                                  scope_desc, thr))
    low = str(out_path).lower()
    is_bg = low.endswith(".bedgraph") or low.endswith(".bg")
    with open(out_path, "w") as fh:
        fh.write(header + "\n")
        if is_bg:
            for _, r in df.iterrows():
                c = float(r["center"])
                start = int(round(c - res / 2.0))
                end = int(round(c + res / 2.0))
                fh.write("%s\t%d\t%d\t%.6g\n" % (r["chrom"], start, end, r["obs"]))
        else:
            fh.write("chrom\tstart\tend\tcenter\tscore_A\tscore_B\t"
                     "delta_abs\tthreshold\tpctile\tsig\n")
            for _, r in df.iterrows():
                c = float(r["center"])
                start = int(round(c - res / 2.0))
                end = int(round(c + res / 2.0))
                cen = int(c) if c.is_integer() else c
                fh.write("%s\t%d\t%d\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%g\t%d\n" % (
                    r["chrom"], start, end, cen,
                    float(r["score_A"]), float(r["score_B"]), float(r["obs"]),
                    float(r["threshold"]), float(r["pctile"]), int(bool(r["sig"]))))


# ---------------------------------------------------------------------
# Per-track KL decomposition (ported from v10)
# ---------------------------------------------------------------------

def parse_categories_json(cats_json_path):
    """Ordered [(name, color), ...] from a BEARING categories JSON; keys are
    integer-string indices matching the kl_0, kl_1, ... TSV columns."""
    with open(cats_json_path) as f:
        d = json.load(f)
    cats = d["categories"]
    return [(cats[k][0], cats[k][1]) for k in sorted(cats.keys(), key=int)]


def _detect_kl_columns(diff_tsv_path):
    with open(diff_tsv_path) as f:
        header = f.readline().rstrip("\n").split("\t")
    kl_cols = sorted([c for c in header if c.startswith("kl_")],
                     key=lambda c: int(c.split("_")[1]))
    return header, kl_cols


def load_diff_kl_table(diff_tsv_path):
    header, kl_cols = _detect_kl_columns(diff_tsv_path)
    if not kl_cols:
        return None, []
    want = [c for c in ("chrom", "start", "end") if c in header] + kl_cols
    df = pd.read_csv(diff_tsv_path, sep="\t", usecols=want, low_memory=False)
    return df, kl_cols


def compute_track_thresholds(df_full, kl_cols, names, pctile_low, pctile_high,
                             scope="genome", chrom=None):
    if scope == "genome":
        df_ref = df_full[df_full["chrom"].isin(AUTOSOMES)] if "chrom" in df_full else df_full
        scope_label, scope_short = "genome (autosomes)", "genome"
    else:
        df_ref = df_full[df_full["chrom"] == chrom] if "chrom" in df_full else df_full
        scope_label, scope_short = chrom, chrom
    thresholds = {}
    print("Per-track {} |kl| thresholds (top {:g}% / top {:g}%):".format(
        scope_label, 100 - pctile_low, 100 - pctile_high))
    for col, name in zip(kl_cols, names):
        abs_vals = df_ref[col].abs().to_numpy()
        finite = abs_vals[np.isfinite(abs_vals)]
        if len(finite) == 0:
            thresholds[col] = (np.inf, np.inf)
            print("  {:<10} (no data)".format(name))
        else:
            tl = float(np.nanpercentile(finite, pctile_low))
            th = float(np.nanpercentile(finite, pctile_high))
            thresholds[col] = (tl, th)
            print("  {:<10} >= {:.3f} (light) / >= {:.3f} (dark)".format(name, tl, th))
    return thresholds, scope_short


def draw_kl_decomposition_track(ax, centers, vals, color, track_name,
                                region_start, region_end, bin_width,
                                thr_low, thr_high, pctile_low, pctile_high,
                                scope_short, show_legend):
    vals = np.asarray(vals, dtype=float)
    centers = np.asarray(centers, dtype=float)
    nonzero = vals != 0
    pos = nonzero & (vals > 0)
    neg = nonzero & (vals < 0)
    ax.bar(centers[pos], vals[pos], width=bin_width,
           color=color, edgecolor="none", alpha=0.9)
    ax.bar(centers[neg], vals[neg], width=bin_width,
           color=color, edgecolor="black", linewidth=0.3, alpha=0.5)
    ax.axhline(0, color="black", lw=0.5)
    abs_vals = np.abs(vals)
    light = (abs_vals >= thr_low) & (abs_vals < thr_high)
    dark = abs_vals >= thr_high
    drew_legend = False
    if light.any():
        ax.scatter(centers[light], np.sign(vals[light]) * (abs_vals[light] + 0.2),
                   marker="v", color="#ff9999", edgecolor="#cc0000",
                   linewidth=0.3, s=18, zorder=10,
                   label=(None if not show_legend
                          else "top {:g}% {}".format(100 - pctile_low, scope_short)))
    if dark.any():
        ax.scatter(centers[dark], np.sign(vals[dark]) * (abs_vals[dark] + 0.2),
                   marker="v", color="#990000", edgecolor="black",
                   linewidth=0.4, s=24, zorder=11,
                   label=(None if not show_legend
                          else "top {:g}% {}".format(100 - pctile_high, scope_short)))
    if show_legend and (light.any() or dark.any()):
        ax.legend(loc="upper right", fontsize=6, ncol=2,
                  handletextpad=0.2, columnspacing=0.5, framealpha=0.8)
        drew_legend = True
    ax.set_xlim(region_start, region_end)
    ax.set_ylabel("{}\nkl".format(track_name), fontsize=8, rotation=0,
                  ha="right", va="center", color=color, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return drew_legend


# ---------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------

def _filter_genes_by_name(genes, whitelist):
    """Curate the gene track down to a whitelist (case-insensitive exact OR
    substring match on the gene name). Returns genes unchanged if whitelist is
    empty, so the default figure is untouched. `genes` is a list of
    (start, end, name, strand) tuples."""
    if not whitelist:
        return genes
    wanted = {str(w).lower() for w in whitelist}
    kept = []
    for g in genes:
        nm = str(g[2] or "").lower()
        if nm in wanted or any(w in nm for w in wanted):
            kept.append(g)
    return kept


def _solve_hic_height(layout_kwargs, data_aspect, hic_height_in,
                      n_hic_rows, iters=4):
    """Find the Hi-C row height (inches) that makes each Hi-C triangle box
    PHYSICALLY match its data aspect (region_w : eff_max_d/2), so the rotated
    diamonds render square instead of stretched.

    The stacked figure shares one genomic x-axis across every panel, so we
    cannot shrink the Hi-C box width to fix the aspect (that would misalign the
    tracks below); the box must instead be given the correct HEIGHT. The exact
    height depends on hspace, the figure margins and the other rows, which are
    awkward to model by hand -- so we build throwaway probe figures with the
    same layout, let matplotlib do the margin math, measure the real rendered
    box, and correct the height by fixed-point iteration. The probes use empty
    axes, so this needs no Hi-C data and never touches the .cool files."""
    h = float(hic_height_in)
    if n_hic_rows == 0 or data_aspect <= 0:
        return h
    for _ in range(max(1, iters)):
        kw = dict(layout_kwargs)
        kw["hic_height_in"] = h
        fig, axes = _combined_figure_layout(**kw)
        try:
            fig.canvas.draw()
            hkey = None
            for cand in ("hic_A", "hic_B", "hic_diff"):
                if cand in axes:
                    hkey = cand
                    break
            if hkey is None:
                return h
            bb = axes[hkey].get_window_extent()
            w_in = bb.width / float(fig.dpi)
            h_in = bb.height / float(fig.dpi)
        finally:
            plt.close(fig)
        if h_in <= 0 or w_in <= 0:
            return h
        box_aspect = w_in / h_in            # >data_aspect => box too wide/short
        h *= box_aspect / data_aspect       # so grow the Hi-C row height
        h = max(1.4, min(7.0, h))
    return h


def _combined_figure_layout(show_rgb, show_hic_list, hic_height_in,
                            has_insul, num_beds, bed_styles,
                            num_decomp, has_pval_fill, fig_width_in=12.0,
                            has_loops_a=False, has_loops_b=False,
                            rgb_slim=False):
    """Build the combined GridSpec. All data panels live in column 0 so
    genomic coordinates stay vertically aligned across the Hi-C, insulation,
    qcat and KL-decomposition panels. The legend spans column 1 over the
    qcat..genes section. Heights are in approximate inches and the figure
    height is their sum, so panels stay consistent regardless of count.
    fig_width_in is auto-scaled by the caller so the Hi-C triangle keeps a
    sensible aspect (wide regions widen the figure rather than squishing)."""
    bed_styles = bed_styles or []
    rgb_h = 1.9 if rgb_slim else 4.3
    loop_h = 0.42
    qcat_h = 0.95
    manhattan_h = 0.85
    pval_fill_h = 0.75
    insul_h = 1.0
    delta_insul_h = 0.85
    gene_h = 0.50
    decomp_h = 0.55
    axis_h = 0.35

    bed_heights = []
    for i in range(num_beds):
        sty = bed_styles[i] if i < len(bed_styles) else "itemRgb"
        bed_heights.append(0.12 if sty == "cbe" else 0.28)

    rows = []
    if show_rgb:
        rows.append(("rgb", rgb_h))
    for k in show_hic_list:
        rows.append(("hic_{}".format(k), hic_height_in))
    if has_loops_a:
        rows.append(("loops_a", loop_h))
    if has_loops_b:
        rows.append(("loops_b", loop_h))
    if has_insul:
        rows.append(("insul", insul_h))
        rows.append(("delta_insul", delta_insul_h))
    rows.append(("qcat_a", qcat_h))
    rows.append(("qcat_b", qcat_h))
    rows.append(("diff", qcat_h))
    rows.append(("manhattan", manhattan_h))
    if has_pval_fill:
        rows.append(("pval_fill", pval_fill_h))
    for i in range(num_beds):
        rows.append(("bed_{}".format(i), bed_heights[i]))
    rows.append(("gene", gene_h))
    for i in range(num_decomp):
        rows.append(("decomp_{}".format(i), decomp_h))
    rows.append(("axis", axis_h))

    height_ratios = [h for _, h in rows]
    total_h = sum(height_ratios) + 1.2
    fig = plt.figure(figsize=(fig_width_in, total_h), dpi=150)
    gs = gridspec.GridSpec(
        len(rows), 2, figure=fig,
        width_ratios=[7.3, 2.2], height_ratios=height_ratios,
        hspace=0.18, wspace=0.06,
        left=0.07, right=0.97, top=0.985, bottom=0.03)
    key_to_index = {key: i for i, (key, _) in enumerate(rows)}
    axes = {key: fig.add_subplot(gs[idx, 0]) for key, idx in key_to_index.items()}
    legend_top = key_to_index["qcat_a"]
    legend_bottom = key_to_index["gene"] + 1
    axes["legend"] = fig.add_subplot(gs[legend_top:legend_bottom, 1])
    return fig, axes


# ---------------------------------------------------------------------
# Main figure builder
# ---------------------------------------------------------------------

def make_combined_figure(
    hic_a_path, hic_b_path, qcat_a_path, qcat_b_path,
    region_str, out_path,
    diff_stats_path=None,
    genes_path=None, gtf_path=None, highlights_path=None,
    diff_qcat_path=None,
    pval_a_path=None, pval_b_path=None, pval_diff_path=None,
    pval_cutoff=0.05,
    label_a="Condition A", label_b="Condition B",
    categories=None, categories_json=None,
    pval_overlay=False,
    rgb_hic=False, rgb_palette="magenta-green",
    rgb_overview=False,
    show_hic=("A", "B", "diff"),
    loops_a_path=None, loops_b_path=None,
    label_genes=None, tidy_labels=False,
    insul_a_path=None, insul_b_path=None, insul_pctile=95.0,
    insul_sig_out=None,
    max_distance=500000, balance=True,
    hic_vmax_arg=None, diff_hic_vmax_arg=None,
    beds=None, bed_style_overrides=None,
    track_pctile_low=99.0, track_pctile_high=99.9,
    threshold_scope="genome", bin_width=200,
    decomp_tracks=None, diff_sign="flip",
    show_pval_fill=False,
):
    chrom, region_start, region_end = parse_region(region_str)
    span = region_end - region_start
    print("\nRegion: {}:{:,}-{:,}  ({:.2f} Mb)".format(
        chrom, region_start, region_end, span / 1e6))
    if categories is None:
        categories = ALL_CATEGORIES

    show_hic_list = [s for s in ("A", "B", "diff") if s in set(show_hic)]
    # RGB triangle: rgb_hic replaces the panels with a big overview; rgb_overview
    # adds a SLIM strip on top of the per-condition panels. Both use the same
    # draw; only the row height differs. Default (both False) => unchanged.
    show_rgb = bool(rgb_hic or rgb_overview)
    rgb_slim = bool(rgb_overview and not rgb_hic)

    # ---- Hi-C via cooler (genomic-coordinate triangles + RGB overview) ----
    M_A = M_B = M_diff = None
    bin_starts = None
    binsize = None
    hic_norm = diff_norm = None
    need_hic = bool(show_hic_list) or show_rgb
    if need_hic:
        print("Loading Hi-C from cool files...")
        M_A, bin_starts, bs_A = fetch_cool_region(
            hic_a_path, chrom, region_start, region_end, balance)
        M_B, bin_starts_B, bs_B = fetch_cool_region(
            hic_b_path, chrom, region_start, region_end, balance)
        if bs_A != bs_B:
            sys.exit("cool binsizes differ: A={} B={}".format(bs_A, bs_B))
        if M_A.shape != M_B.shape:
            sys.exit("cool matrix shapes differ: {} vs {}".format(M_A.shape, M_B.shape))
        binsize = bs_A
        # Diff in the figure's A - B direction (reference minus condition;
        # matches qcat diff / lollipops / decomposition).
        M_diff = M_A - M_B
        print("  binsize={}, {} bins".format(binsize, M_A.shape[0]))

        if hic_vmax_arg is None:
            posv = M_A[np.isfinite(M_A) & (M_A > 0)]
            hic_vmax = float(np.nanpercentile(posv, 99)) if len(posv) else 1.0
        else:
            hic_vmax = hic_vmax_arg
        hic_norm = Normalize(vmin=0, vmax=hic_vmax)

        if diff_hic_vmax_arg is None:
            d_abs = np.abs(M_diff[np.isfinite(M_diff)])
            diff_vmax = float(np.nanpercentile(d_abs, 99)) if len(d_abs) else 1.0
        else:
            diff_vmax = diff_hic_vmax_arg
        if diff_vmax <= 0:
            diff_vmax = max(hic_vmax * 0.5, 1e-9)
        diff_norm = TwoSlopeNorm(vmin=-diff_vmax, vcenter=0, vmax=diff_vmax)

    # ---- qcat epilogos ----
    print("Loading BEARING epilogos scores...")
    pos_a, scores_a, ns_a = load_qcat_scores(qcat_a_path, chrom, region_start, region_end)
    pos_b, scores_b, ns_b = load_qcat_scores(qcat_b_path, chrom, region_start, region_end)
    num_states = max(ns_a, ns_b)
    if scores_a.shape[1] < num_states:
        scores_a = np.pad(scores_a, ((0, 0), (0, num_states - scores_a.shape[1])))
    if scores_b.shape[1] < num_states:
        scores_b = np.pad(scores_b, ((0, 0), (0, num_states - scores_b.shape[1])))
    print("  Cond A: {} bins, Cond B: {} bins, {} states".format(
        len(pos_a), len(pos_b), num_states))

    # ---- genes / highlights ----
    genes = None
    if genes_path:
        genes = load_genes(genes_path, chrom, region_start, region_end)
        print("Loaded {} gene records".format(len(genes)))
    elif gtf_path:
        genes = load_genes_gtf(gtf_path, chrom, region_start, region_end)
        print("Loaded {} gene records from GTF".format(len(genes)))
    if genes is not None and label_genes:
        _n0 = len(genes)
        genes = _filter_genes_by_name(genes, label_genes)
        print("  curated gene track: {} -> {} genes ({})".format(
            _n0, len(genes), ", ".join(label_genes)))
    highlights = None
    if highlights_path:
        highlights = load_highlights(highlights_path, chrom, region_start, region_end)

    # ---- loop anchors (optional; drawn as arc tracks under the Hi-C) ----
    loops_a = load_loops(loops_a_path, chrom, region_start, region_end) if loops_a_path else None
    loops_b = load_loops(loops_b_path, chrom, region_start, region_end) if loops_b_path else None
    has_loops_a = loops_a_path is not None
    has_loops_b = loops_b_path is not None
    def _loops_msg(loops, path, label):
        n = len(loops or [])
        if n > 0:
            print("Loaded {} {} loop(s) in region".format(n, label))
            return
        # 0 loops: report the file's chromosome naming so a chr6-vs-6 or
        # wrong-region mismatch is obvious instead of a silently empty track.
        seen = []
        try:
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    p = line.split("\t")
                    if len(p) < 6:
                        p = line.split()
                    if len(p) >= 6:
                        seen.append(p[0])
                    if len(seen) >= 200:
                        break
        except OSError:
            print("  WARNING: could not open {} loops file {}".format(label, path))
            return
        chroms = sorted(set(seen))[:8]
        print("  WARNING: 0 {} loops in {} for {} -- file chroms seen: {} "
              "(want {})".format(label, path, chrom, chroms, chrom))

    if has_loops_a:
        _loops_msg(loops_a, loops_a_path, label_a)
    if has_loops_b:
        _loops_msg(loops_b, loops_b_path, label_b)

    # ---- diff qcat ----
    pos_diff, scores_diff, ns_diff = [], np.zeros((0, 1), dtype=np.float32), 1
    has_diff = diff_qcat_path is not None
    if has_diff:
        pos_diff, scores_diff, ns_diff = load_qcat_scores(
            diff_qcat_path, chrom, region_start, region_end)
        # File is natively label_a - label_b (reference - condition); keep it.
        ns_full = max(ns_diff, num_states)
        if scores_diff.shape[1] < ns_full:
            scores_diff = np.pad(scores_diff, ((0, 0), (0, ns_full - scores_diff.shape[1])))
        print("Loaded diff qcat: {} bins".format(len(pos_diff)))

    # ---- diff p-value (Manhattan + optional fill) ----
    pos_pval_diff, vals_pval_diff = None, None
    pval_diff_kl_scores = None
    if pval_diff_path is not None:
        if str(pval_diff_path).endswith(".tsv"):
            (pos_pval_diff, vals_pval_diff,
             pval_diff_kl_scores, _) = load_stats_tsv_with_categories(
                pval_diff_path, chrom, region_start, region_end)
        else:
            pos_pval_diff, vals_pval_diff = load_pval_track_values(
                pval_diff_path, chrom, region_start, region_end)
        # Keep native label_a - label_b (reference - condition) sign.
        print("Loaded diff p-value: {} bins".format(len(pos_pval_diff)))

    pos_pval_a = vals_pval_a = pos_pval_b = vals_pval_b = None
    if pval_overlay:
        if pval_a_path is not None:
            pos_pval_a, vals_pval_a = load_pval_track_values(
                pval_a_path, chrom, region_start, region_end)
        if pval_b_path is not None:
            pos_pval_b, vals_pval_b = load_pval_track_values(
                pval_b_path, chrom, region_start, region_end)
    pval_cutoff_value = -math.log10(pval_cutoff) if (pval_cutoff and pval_cutoff > 0) else None

    # ---- insulation (+ stats-matched significance) ----
    ins_A = ins_B = None
    insul_sig_df = None
    has_insul = insul_a_path is not None and insul_b_path is not None
    do_insul_sig = 0 < insul_pctile < 100
    if has_insul:
        print("Loading insulation...")
        ins_A = load_insul_bm(insul_a_path, chrom, region_start, region_end)
        ins_B = load_insul_bm(insul_b_path, chrom, region_start, region_end)
        if do_insul_sig:
            scope_label = ("genome (autosomes)" if threshold_scope == "genome" else chrom)
            print("Computing {} |delta insul| {}th percentile...".format(
                scope_label, insul_pctile))
            if threshold_scope == "genome":
                ins_A_ref = load_insul_bm_genome(insul_a_path)
                ins_B_ref = load_insul_bm_genome(insul_b_path)
            else:
                ins_A_ref = load_insul_bm_chr(insul_a_path, chrom)
                ins_B_ref = load_insul_bm_chr(insul_b_path, chrom)
            insul_sig_df = compute_insul_pctile_threshold(
                ins_A_ref, ins_B_ref, pctile=insul_pctile)
            # Export the FULL-scope table (before the plotted-chromosome filter
            # below), so the exported sig/threshold reflect the scope the
            # threshold was computed over. Export-only; does not touch plotting.
            if insul_sig_out:
                _res_insul = int(round(
                    (ins_A_ref["end"] - ins_A_ref["start"]).median()))
                write_insul_sig_table(
                    insul_sig_df, insul_sig_out, _res_insul,
                    threshold_scope, chrom, insul_a_path, insul_b_path,
                    insul_pctile)
                print("  Wrote |delta insul| table ({} bins) -> {}".format(
                    len(insul_sig_df), insul_sig_out))
            insul_sig_df = insul_sig_df[insul_sig_df["chrom"] == chrom].reset_index(drop=True)
            thr_val = insul_sig_df["threshold"].iloc[0] if len(insul_sig_df) else float("nan")
            print("  threshold |delta insul| = {:.3f}".format(thr_val))

    # ---- per-track KL decomposition source ----
    decomp_source = None
    if diff_stats_path is not None:
        decomp_source = diff_stats_path
    elif pval_diff_path is not None and str(pval_diff_path).endswith(".tsv"):
        decomp_source = pval_diff_path

    decomp_names, decomp_colors, decomp_region_vals = [], [], []
    decomp_centers = None
    track_thresholds = {}
    scope_short = "genome" if threshold_scope == "genome" else chrom
    if decomp_source is not None and categories_json is not None:
        print("Loading per-track KL decomposition from {}...".format(decomp_source))
        cat_pairs = parse_categories_json(categories_json)
        all_names = [c[0] for c in cat_pairs]
        all_colors = [c[1] for c in cat_pairs]
        df_full, kl_cols = load_diff_kl_table(decomp_source)
        if df_full is None or not kl_cols:
            print("  WARNING: no kl_* columns; skipping decomposition.")
        else:
            n_use = min(len(kl_cols), len(all_names))
            if len(kl_cols) != len(all_names):
                print("  WARNING: {} kl_* cols vs {} categories; using first {}.".format(
                    len(kl_cols), len(all_names), n_use))
            kl_cols = kl_cols[:n_use]
            sel_names = all_names[:n_use]
            sel_colors = all_colors[:n_use]
            if decomp_tracks:
                keep = [i for i, nm in enumerate(sel_names) if nm in set(decomp_tracks)]
                if not keep:
                    print("  WARNING: --decomp-tracks matched nothing; using all.")
                    keep = list(range(len(sel_names)))
                kl_cols = [kl_cols[i] for i in keep]
                sel_names = [sel_names[i] for i in keep]
                sel_colors = [sel_colors[i] for i in keep]
            track_thresholds, scope_short = compute_track_thresholds(
                df_full, kl_cols, sel_names, track_pctile_low, track_pctile_high,
                scope=threshold_scope, chrom=chrom)
            mask = ((df_full["chrom"] == chrom)
                    & (df_full["end"] > region_start)
                    & (df_full["start"] < region_end))
            df_r = df_full[mask].copy().reset_index(drop=True)
            decomp_centers = (df_r["start"].to_numpy() + df_r["end"].to_numpy()) / 2.0
            # Decomposition is fixed to the A - B convention (raw diff_i =
            # KL_A_i - KL_B_i; positive = more active in condition A), so signal
            # stronger in A (the reference, e.g. DN) points ABOVE the 0-line.
            # This is intentionally independent of diff_sign, which only governs
            # the qcat-diff / lollipop panels.
            for col in kl_cols:
                decomp_region_vals.append(df_r[col].to_numpy(dtype=float))
            decomp_names, decomp_colors = sel_names, sel_colors
            print("  {} decomposition tracks, {} bins in region".format(
                len(decomp_names), len(df_r)))
    elif decomp_source is not None and categories_json is None:
        print("  NOTE: decomposition needs --categories (JSON cats file); skipping.")
    num_decomp = len(decomp_names)

    # ---- BED overlays ----
    beds = beds or []
    bed_style_overrides = bed_style_overrides or {}
    bed_features_list, bed_styles, bed_paths_kept = [], [], []
    for bed_path in beds:
        feats = load_bed_for_region(bed_path, chrom, region_start, region_end)
        if not feats:
            # Locus-specific annotation (e.g. CBE at Tcrb, AgR genes at the
            # antigen-receptor loci) -- skip the track entirely where it has no
            # features rather than drawing an empty row.
            print("  bed {}: no features in region, skipping track".format(
                Path(bed_path).name))
            continue
        k1, k2 = str(bed_path), Path(bed_path).name
        if k1 in bed_style_overrides:
            style = bed_style_overrides[k1]
        elif k2 in bed_style_overrides:
            style = bed_style_overrides[k2]
        elif any(f.get("item_rgb") is not None for f in feats) and len(feats) > 0:
            style = "itemRgb"
        else:
            style = "cbe"
        bed_features_list.append(feats)
        bed_styles.append(style)
        bed_paths_kept.append(bed_path)

    # ---- layout ----
    # Auto-scale so the Hi-C triangle fills its panel without squishing or
    # leaving whitespace. The triangle only rises to (eff_max_d/2), where
    # eff_max_d caps the contact distance at the region's own span -- so for a
    # small region the panel must be SHORTER, and for a wide region WIDER.
    region_w = region_end - region_start
    eff_max_d = min(max_distance, region_w) if region_w > 0 else max_distance
    # True triangle data aspect width:height = region_w : (eff_max_d/2). This is
    # the SAME extent plot_hic_triangle renders (it caps ylim at eff_max_d/2), so
    # matching the physical box to it makes the diamonds square.
    data_aspect = (region_w / (eff_max_d / 2.0)) if eff_max_d > 0 else 2.0
    # Fix the figure WIDTH from a sensible nominal box, then solve the Hi-C row
    # HEIGHT so the box aspect equals data_aspect (see _solve_hic_height). We
    # never squeeze the width, so the genomic axis stays aligned with the tracks
    # below; only the Hi-C rows grow/shrink vertically.
    hic_height_nominal = 3.0
    col0_w = hic_height_nominal * min(max(data_aspect, 1.2), 6.0)
    col0_w = max(7.0, min(16.0, col0_w))             # sane Hi-C column width
    fig_width_in = col0_w * (7.3 + 2.2) / 7.3
    fig_width_in = max(11.0, min(22.0, fig_width_in))
    layout_kwargs = dict(
        show_rgb=show_rgb, show_hic_list=show_hic_list,
        hic_height_in=hic_height_nominal,
        has_insul=has_insul, num_beds=len(bed_features_list), bed_styles=bed_styles,
        num_decomp=num_decomp,
        has_pval_fill=show_pval_fill and (pos_pval_diff is not None),
        fig_width_in=fig_width_in,
        has_loops_a=has_loops_a, has_loops_b=has_loops_b, rgb_slim=rgb_slim)
    hic_height_in = _solve_hic_height(
        layout_kwargs, data_aspect, hic_height_nominal, len(show_hic_list))
    layout_kwargs["hic_height_in"] = hic_height_in
    fig, axes = _combined_figure_layout(**layout_kwargs)

    # ---- combined RGB triangle (optional overview) ----
    if show_rgb:
        print("Rendering combined RGB Hi-C triangle...")
        rgb_image = make_rgb_hic(np.nan_to_num(M_A, nan=0.0),
                                 np.nan_to_num(M_B, nan=0.0), palette=rgb_palette,
                                 joint_norm=True)
        ax_rgb = axes["rgb"]
        ax_rgb.set_axis_off()
        _draw_rgb_triangle(ax_rgb, rgb_image, inverted=False)
        # Label colors follow the palette's A/B channels so they match the image.
        _rgb_label_colors = {
            "red-green": ("#cc0000", "#00a83a"),
            "magenta-green": ("#cc00cc", "#00a83a"),
            "magenta-green-white": ("#cc00cc", "#00a83a"),
            "blue-red": ("#1f4fd8", "#cc0000"),
            "green-blue": ("#00a83a", "#1f4fd8"),
        }
        _ca, _cb = _rgb_label_colors.get(rgb_palette, ("#cc0000", "#00a83a"))
        ax_rgb.text(0.06, 0.94, label_a, transform=ax_rgb.transAxes, fontsize=8,
                    color=_ca, ha="left", va="top", fontweight="bold")
        ax_rgb.text(0.94, 0.94, label_b, transform=ax_rgb.transAxes, fontsize=8,
                    color=_cb, ha="right", va="top", fontweight="bold")

    # ---- genomic Hi-C triangles (single A, single B, diff) ----
    for k in show_hic_list:
        ax_h = axes["hic_{}".format(k)]
        if k == "A":
            pc = plot_hic_triangle(ax_h, M_A, bin_starts, binsize,
                                   region_start, region_end, max_distance,
                                   HIC_CMAP, hic_norm,
                                   label="{}\nHi-C".format(label_a))
            _add_inset_colorbar(ax_h, pc, "contact")
        elif k == "B":
            pc = plot_hic_triangle(ax_h, M_B, bin_starts, binsize,
                                   region_start, region_end, max_distance,
                                   HIC_CMAP, hic_norm,
                                   label="{}\nHi-C".format(label_b))
            _add_inset_colorbar(ax_h, pc, "contact")
        else:  # diff (B - A)
            pc = plot_hic_triangle(ax_h, M_diff, bin_starts, binsize,
                                   region_start, region_end, max_distance,
                                   "RdBu_r", diff_norm,
                                   label="{} - {}\nHi-C diff".format(label_a, label_b))
            _add_inset_colorbar(ax_h, pc, "{} - {}".format(label_a, label_b))
        ax_h.tick_params(axis="x", labelbottom=False)

    # ---- loop-anchor arc tracks (optional) ----
    # NOTE: draw_loops_horizontal works in axis-fraction coords (it sets
    # xlim(0,1) and maps genomic->[0,1] internally). Do NOT override the xlim
    # to genomic coordinates here -- that squeezes the arcs into an invisible
    # sliver at x~0 (this is exactly why the loop track came up empty).
    if "loops_a" in axes:
        draw_loops_horizontal(
            axes["loops_a"], loops_a or [], region_start, region_end,
            highlights=highlights, label="{} loops".format(label_a),
            color="#d62728", anchor_color="#d62728")
        axes["loops_a"].tick_params(axis="x", labelbottom=False)
    if "loops_b" in axes:
        draw_loops_horizontal(
            axes["loops_b"], loops_b or [], region_start, region_end,
            highlights=highlights, label="{} loops".format(label_b),
            color="#2ca02c", anchor_color="#2ca02c")
        axes["loops_b"].tick_params(axis="x", labelbottom=False)

    # ---- insulation panels ----
    if has_insul:
        ax_ins = axes["insul"]
        ax_ins.plot(ins_A["center"], ins_A["score"], color="#1f77b4",
                    label=label_a, lw=1.2)
        ax_ins.plot(ins_B["center"], ins_B["score"], color="#d62728",
                    label=label_b, lw=1.2)
        ax_ins.axhline(0, color="gray", lw=0.5, ls="--")
        ax_ins.legend(loc="upper right", fontsize=7)
        ax_ins.set_xlim(region_start, region_end)
        ax_ins.set_ylabel("Insulation", fontsize=8, rotation=0, ha="right", va="center")
        ax_ins.tick_params(axis="x", labelbottom=False)
        for sp in ("top", "right"):
            ax_ins.spines[sp].set_visible(False)

        ax_di = axes["delta_insul"]
        merged = pd.merge(ins_A[["center", "score"]], ins_B[["center", "score"]],
                          on="center", suffixes=("_A", "_B"))
        delta = (merged["score_A"] - merged["score_B"]).abs()
        ax_di.fill_between(merged["center"], 0, delta, color="steelblue", alpha=0.6)
        ax_di.plot(merged["center"], delta, color="steelblue", lw=0.9)
        ax_di.set_xlim(region_start, region_end)
        ax_di.set_ylabel("|delta\ninsul|", fontsize=8, rotation=0, ha="right", va="center")
        ax_di.tick_params(axis="x", labelbottom=False)
        for sp in ("top", "right"):
            ax_di.spines[sp].set_visible(False)
        if insul_sig_df is not None and len(insul_sig_df) > 0:
            thr = float(insul_sig_df["threshold"].iloc[0])
            ax_di.axhline(thr, color="gold", lw=0.8, ls="--", alpha=0.8, zorder=5)
            sub = insul_sig_df[(insul_sig_df["center"] >= region_start)
                               & (insul_sig_df["center"] <= region_end)
                               & insul_sig_df["sig"]]
            if not sub.empty:
                y_off = max(0.05, float(np.nanmax(delta)) * 0.08) if len(delta) else 0.05
                ax_di.scatter(sub["center"], sub["obs"] + y_off, marker="*",
                              color="gold", s=40, edgecolor="black", linewidth=0.5,
                              zorder=10, label="top {:g}% {}".format(
                                  100 - insul_pctile, scope_short))
                ax_di.legend(loc="upper right", fontsize=7)

    # ---- shared y-scale for qcat A/B/diff ----
    y_max_shared = 0.0
    if scores_a.size:
        y_max_shared = max(y_max_shared, float(scores_a.max()))
    if scores_b.size:
        y_max_shared = max(y_max_shared, float(scores_b.max()))
    if has_diff and scores_diff.size:
        y_max_shared = max(y_max_shared, float(np.abs(scores_diff).max()))
    y_max_shared = 1.0 if y_max_shared == 0 else y_max_shared * 1.05

    # ---- qcat A / B ----
    print("Rendering qcat tracks...")
    if pval_overlay and pos_pval_a is not None and vals_pval_a is not None:
        draw_epilogos_with_pval_horizontal(
            axes["qcat_a"], pos_a, scores_a, num_states, region_start, region_end,
            categories[:num_states], pval_positions=pos_pval_a, pval_values=vals_pval_a,
            pval_alpha=pval_cutoff, highlights=highlights, label=label_a, y_max=y_max_shared)
    else:
        draw_epilogos_horizontal(
            axes["qcat_a"], pos_a, scores_a, num_states, region_start, region_end,
            categories[:num_states], highlights=highlights, label=label_a, y_max=y_max_shared)
    if pval_overlay and pos_pval_b is not None and vals_pval_b is not None:
        draw_epilogos_with_pval_horizontal(
            axes["qcat_b"], pos_b, scores_b, num_states, region_start, region_end,
            categories[:num_states], pval_positions=pos_pval_b, pval_values=vals_pval_b,
            pval_alpha=pval_cutoff, highlights=highlights, label=label_b, y_max=y_max_shared)
    else:
        draw_epilogos_horizontal(
            axes["qcat_b"], pos_b, scores_b, num_states, region_start, region_end,
            categories[:num_states], highlights=highlights, label=label_b, y_max=y_max_shared)

    # ---- diff qcat ----
    if has_diff:
        _diff_lbl = ("diff qcat\n{}-{}".format(label_a, label_b) if tidy_labels
                     else "diff qcat " + label_a + " - " + label_b)
        draw_diff_horizontal(
            axes["diff"], pos_diff, scores_diff, num_states, region_start, region_end,
            categories[:num_states], highlights=highlights,
            label=_diff_lbl, diff_max=y_max_shared)
    else:
        axes["diff"].set_axis_off()

    # ---- Manhattan lollipops ----
    has_pm = pos_pval_diff is not None and len(pos_pval_diff) > 0
    if has_pm:
        y_max_pm = float(np.abs(vals_pval_diff).max()) if vals_pval_diff.size else 0.0
        y_max_pm = 1.0 if y_max_pm == 0 else y_max_pm * 1.05
        _pm_lbl = ("diff p\n{}-{}".format(label_a, label_b) if tidy_labels
                   else "diff p-value " + label_a + " - " + label_b)
        draw_pval_manhattan_horizontal(
            axes["manhattan"], pos_pval_diff, vals_pval_diff, region_start, region_end,
            highlights=highlights, label=_pm_lbl,
            score_positions=pos_a, score_matrix=scores_a, categories=categories[:num_states],
            y_max=y_max_pm, cutoff_value=pval_cutoff_value, kl_scores=pval_diff_kl_scores)
    else:
        axes["manhattan"].set_axis_off()

    # ---- optional filled signed diff p-value ----
    if "pval_fill" in axes:
        if has_pm:
            y_max_pf = float(np.abs(vals_pval_diff).max()) if vals_pval_diff.size else 0.0
            y_max_pf = 1.0 if y_max_pf == 0 else y_max_pf * 1.05
            draw_pval_diff_horizontal(
                axes["pval_fill"], pos_pval_diff, vals_pval_diff, region_start, region_end,
                highlights=highlights, label="diff p-value (fill) " + label_a + " - " + label_b,
                y_max=y_max_pf, cutoff_value=pval_cutoff_value,
                diff_score_positions=pos_diff if has_diff else None,
                diff_score_matrix=scores_diff if has_diff else None,
                categories=categories[:num_states] if has_diff else None)
        else:
            axes["pval_fill"].set_axis_off()

    # ---- BED rows ----
    for i, bed_path in enumerate(bed_paths_kept):
        ax_bed = axes.get("bed_{}".format(i))
        if ax_bed is None:
            continue
        draw_bed_track(ax_bed, bed_features_list[i], region_start, region_end,
                       style=bed_styles[i], label=Path(bed_path).stem)

    # ---- genes ----
    if genes is not None and len(genes) > 0:
        draw_gene_track(axes["gene"], genes, region_start, region_end,
                        highlights=highlights, label="Genes")
    else:
        axes["gene"].set_axis_off()

    # ---- per-track KL decomposition ----
    if num_decomp > 0:
        print("Rendering per-track KL decomposition ({} tracks)...".format(num_decomp))
        kl_col_order = list(track_thresholds.keys())
        legend_drawn = False
        for i in range(num_decomp):
            ax_d = axes["decomp_{}".format(i)]
            col = kl_col_order[i] if i < len(kl_col_order) else None
            thr_low, thr_high = track_thresholds.get(col, (np.inf, np.inf))
            drew = draw_kl_decomposition_track(
                ax_d, decomp_centers, decomp_region_vals[i], decomp_colors[i],
                decomp_names[i], region_start, region_end, bin_width,
                thr_low, thr_high, track_pctile_low, track_pctile_high,
                scope_short, show_legend=(not legend_drawn))
            if drew:
                legend_drawn = True
            ax_d.tick_params(axis="x", labelbottom=False)
            if i == 0:
                # Make the A - B direction explicit: up = more in A (e.g. DN).
                ax_d.annotate(
                    "+ {}".format(label_a), xy=(0.004, 0.97),
                    xycoords="axes fraction", fontsize=6, ha="left", va="top",
                    color="#444444")
                ax_d.annotate(
                    "- {}".format(label_b), xy=(0.004, 0.03),
                    xycoords="axes fraction", fontsize=6, ha="left", va="bottom",
                    color="#444444")

    # ---- genomic axis (bottom) ----
    draw_genomic_axis(axes["axis"], region_start, region_end, chrom)

    # ---- legend ----
    print("Rendering legend...")
    draw_legend(axes["legend"], categories[:num_states], num_states,
                hic_vmax=None, rgb_mode=rgb_hic, label_a=label_a, label_b=label_b)

    print("Saving to {}...".format(out_path))
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Done.  Figure: {}".format(out_path))


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=("Combined Hi-C + BEARING figure: single/diff Hi-C triangles, "
                     "insulation, qcat epilogos tracks, Manhattan lollipops, and "
                     "per-track KL decomposition with top-1%/0.1% markers."),
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    # triangle.py-compatible surface
    parser.add_argument("--contact-a", required=True, metavar="COOL")
    parser.add_argument("--contact-b", required=True, metavar="COOL")
    parser.add_argument("--format", dest="hic_fmt", choices=["cool", "hic"],
                        default=None, metavar="FMT",
                        help="Accepted for batch compatibility; this script reads "
                             "cool directly via cooler.")
    parser.add_argument("--qcat-a", required=True, metavar="FILE")
    parser.add_argument("--qcat-b", required=True, metavar="FILE")
    parser.add_argument("--diff-qcat", metavar="FILE", default=None)
    parser.add_argument("--pval-a", metavar="FILE", default=None)
    parser.add_argument("--pval-b", metavar="FILE", default=None)
    parser.add_argument("--pval-diff", metavar="FILE", default=None)
    parser.add_argument("--pval-cutoff", type=float, default=0.05, metavar="P")
    parser.add_argument("--pval-overlay", action="store_true")
    parser.add_argument("--rgb-hic", action="store_true",
                        help="Also draw the combined magenta/green RGB overview "
                             "triangle above the single/diff Hi-C panels.")
    parser.add_argument("--rgb-overview", action="store_true",
                        help="Add a SLIM RGB overview strip on top of the "
                             "quantitative single/diff Hi-C panels (does not "
                             "replace them; good for a compact overview).")
    parser.add_argument("--rgb-palette",
                        choices=["magenta-green", "red-green", "blue-red",
                                 "green-blue", "magenta-green-white"],
                        default="magenta-green")
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument("--region", metavar="CHR:START-END")
    region_group.add_argument("--regions-file", metavar="TSV")
    parser.add_argument("--resolution", type=int, default=10000, metavar="BP",
                        help="Informational only; binsize comes from the cool file.")
    parser.add_argument("--out", metavar="FILE")
    parser.add_argument("--outdir", metavar="DIR", default=".")
    parser.add_argument("--label-a", default="Condition A", metavar="STR")
    parser.add_argument("--label-b", default="Condition B", metavar="STR")
    parser.add_argument("--genes", metavar="BED")
    parser.add_argument("--gtf", metavar="GTF")
    parser.add_argument("--label-genes", nargs="+", default=None, metavar="NAME",
                        help="Curate the gene track to only these gene names "
                             "(case-insensitive, substring match), to de-clutter "
                             "dense loci. Default: show all genes.")
    parser.add_argument("--tidy-labels", action="store_true",
                        help="Use compact two-line panel labels for the "
                             "differential tracks (reduces left-margin crowding).")
    parser.add_argument("--highlights", metavar="BED")
    parser.add_argument("--categories", metavar="JSON/YAML",
                        help="Category names/colors. A JSON cats file is required "
                             "for the KL decomposition (maps kl_* columns).")
    parser.add_argument("--bed", action="append", default=[], metavar="FILE")
    parser.add_argument("--bed-style", action="append", default=[], metavar="FILE=STYLE")
    parser.add_argument("--loops-a", metavar="BEDPE", default=None,
                        help="BEDPE loop calls for condition A; drawn as an arc "
                             "track under the Hi-C panels.")
    parser.add_argument("--loops-b", metavar="BEDPE", default=None,
                        help="BEDPE loop calls for condition B; drawn as an arc "
                             "track under the Hi-C panels.")
    # Accepted for batch drop-in compatibility (single combined loop file); the
    # per-condition --loops-a/--loops-b above are what get drawn.
    parser.add_argument("--loops", metavar="BEDPE", default=None, help=argparse.SUPPRESS)

    # Hi-C panels (genomic triangles)
    parser.add_argument("--show-hic", default="A,B,diff",
                        help="Subset of A,B,diff genomic Hi-C triangles to draw "
                             "(default A,B,diff). Use '' for none.")
    parser.add_argument("--max-distance", type=int, default=500000,
                        help="Max Hi-C contact distance to plot (bp, default 500 kb).")
    parser.add_argument("--no-balance", action="store_true",
                        help="Skip cooler KR balancing (default: balanced).")
    parser.add_argument("--hic-vmax", type=float, default=None,
                        help="Color-scale max for single Hi-C (default 99th pctile).")
    parser.add_argument("--diff-hic-vmax", type=float, default=None,
                        help="Color-scale max for diff Hi-C (default 99th pctile of |B-A|).")

    # Insulation
    parser.add_argument("--insul-A", dest="insul_a", metavar="BM", default=None,
                        help="HiCExplorer .bm TAD/insulation score for condition A.")
    parser.add_argument("--insul-B", dest="insul_b", metavar="BM", default=None,
                        help="HiCExplorer .bm TAD/insulation score for condition B.")
    parser.add_argument("--insul-pctile", type=float, default=95.0,
                        help="Percentile threshold for |delta insul| significance "
                             "stars (default 95.0; 0 or 100 disables).")
    parser.add_argument("--insul-sig-out", dest="insul_sig_out", metavar="PATH",
                        default=None,
                        help="Export the per-bin |delta insulation| table (the "
                             "data behind the significance stars) to PATH as TSV. "
                             "If PATH ends in .bedgraph or .bg, write a 4-column "
                             "bedgraph of the |delta| values instead. Export-only; "
                             "does not change the figure. Table covers the full "
                             "--threshold-scope (genome or chrom), before the "
                             "plotted-region restriction.")

    # KL decomposition
    parser.add_argument("--diff-stats", metavar="TSV", default=None,
                        help="Differential stats TSV (kl_* cols) for the decomposition. "
                             "Defaults to --pval-diff when that is a .tsv.")
    parser.add_argument("--track-pctile-low", type=float, default=99.0,
                        help="Light-marker percentile (default 99.0 = top 1%%).")
    parser.add_argument("--track-pctile-high", type=float, default=99.9,
                        help="Dark-marker percentile (default 99.9 = top 0.1%%).")
    parser.add_argument("--threshold-scope", choices=["genome", "chrom"],
                        default="genome",
                        help="Reference distribution for the per-track AND insulation "
                             "percentile thresholds (default genome autosomes).")
    parser.add_argument("--bin-width", type=int, default=200,
                        help="BEARING bin width for decomposition bars (bp).")
    parser.add_argument("--decomp-tracks", nargs="+", default=None,
                        help="Subset of category names for the decomposition.")
    parser.add_argument("--diff-sign", choices=["flip", "asis"], default="flip",
                        help="Deprecated no-op. The per-track decomposition and the "
                             "diff qcat now share the native label_a - label_b "
                             "(reference - condition) sign; this flag is ignored.")
    parser.add_argument("--pval-fill", action="store_true",
                        help="Also show the filled signed diff p-value track.")

    args = parser.parse_args()

    if not (0 < args.pval_cutoff <= 1):
        parser.error("--pval-cutoff must be > 0 and <= 1")
    if args.genes and args.gtf:
        parser.error("Use only one of --genes or --gtf")
    if args.region and not args.out:
        parser.error("--out is required when using --region")

    show_hic = [s.strip() for s in args.show_hic.split(",") if s.strip()]
    for s in show_hic:
        if s not in ("A", "B", "diff"):
            parser.error("--show-hic values must be in {A, B, diff}: " + s)

    cli_categories = None
    categories_json = None
    if args.categories:
        cli_categories, _ = load_categories_yaml(args.categories)
        print("Categories loaded from: {}  ({} states)".format(
            args.categories, len(cli_categories)))
        if str(args.categories).lower().endswith(".json"):
            categories_json = args.categories
        else:
            print("  NOTE: --categories is not .json; KL decomposition will be skipped.")

    print("\n" + "=" * 60)
    print("  bearing_hic_combined_plot.py")
    print("  Cond A:     {}  ({})".format(args.label_a, args.contact_a))
    print("  Cond B:     {}  ({})".format(args.label_b, args.contact_b))
    print("  Hi-C panels: {}{}".format(
        ",".join(show_hic) if show_hic else "(none)",
        " + RGB overview" if args.rgb_hic else ""))
    if args.insul_a and args.insul_b:
        print("  Insulation: A={}  B={}  (sig top {:g}%)".format(
            args.insul_a, args.insul_b, 100 - args.insul_pctile))
    if args.gtf:
        print("  GTF:        {}".format(args.gtf))
    if args.diff_qcat:
        print("  Diff qcat:  {}".format(args.diff_qcat))
    if args.pval_diff:
        print("  P-val diff: {}".format(args.pval_diff))
    if args.diff_stats:
        print("  Diff stats: {}".format(args.diff_stats))
    print("  Decomp:     top {:g}% / top {:g}% ({} scope)".format(
        100 - args.track_pctile_low, 100 - args.track_pctile_high, args.threshold_scope))
    print("  Max dist:   {:,} bp".format(args.max_distance))
    print("=" * 60 + "\n")

    bed_style_overrides = {}
    for s in args.bed_style:
        if "=" not in s:
            continue
        left, right = s.split("=", 1)
        bed_style_overrides[left] = "itemRgb" if right.lower() == "itemrgb" else "cbe"

    common = dict(
        qcat_a_path=args.qcat_a, qcat_b_path=args.qcat_b,
        diff_stats_path=args.diff_stats,
        genes_path=args.genes, gtf_path=args.gtf, highlights_path=args.highlights,
        diff_qcat_path=args.diff_qcat,
        pval_a_path=args.pval_a, pval_b_path=args.pval_b, pval_diff_path=args.pval_diff,
        pval_cutoff=args.pval_cutoff,
        label_a=args.label_a, label_b=args.label_b,
        categories=cli_categories, categories_json=categories_json,
        pval_overlay=args.pval_overlay,
        rgb_hic=args.rgb_hic, rgb_palette=args.rgb_palette,
        rgb_overview=args.rgb_overview,
        show_hic=show_hic,
        loops_a_path=args.loops_a, loops_b_path=args.loops_b,
        label_genes=args.label_genes, tidy_labels=args.tidy_labels,
        insul_a_path=args.insul_a, insul_b_path=args.insul_b, insul_pctile=args.insul_pctile,
        insul_sig_out=args.insul_sig_out,
        max_distance=args.max_distance, balance=not args.no_balance,
        hic_vmax_arg=args.hic_vmax, diff_hic_vmax_arg=args.diff_hic_vmax,
        beds=args.bed, bed_style_overrides=bed_style_overrides,
        track_pctile_low=args.track_pctile_low, track_pctile_high=args.track_pctile_high,
        threshold_scope=args.threshold_scope, bin_width=args.bin_width,
        decomp_tracks=args.decomp_tracks, diff_sign=args.diff_sign,
        show_pval_fill=args.pval_fill,
    )

    if args.regions_file:
        regions = load_regions_file(args.regions_file)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        print("Batch mode: {} region(s) -> {}\n".format(len(regions), outdir))
        for i, reg in enumerate(regions, start=1):
            name = reg["name"]
            region_str = reg["region"]
            out_path = (outdir / reg["out"]) if reg["out"] else (
                outdir / "{}_combined.pdf".format(_safe_filename(name)))
            print("[{}/{}] {}: {} -> {}".format(i, len(regions), name, region_str, out_path))
            try:
                make_combined_figure(
                    hic_a_path=args.contact_a, hic_b_path=args.contact_b,
                    region_str=region_str, out_path=out_path, **common)
            except Exception as e:
                print("  ERROR: failed region {}: {}".format(name, e))
                continue
        print("\nBatch complete. {} figure(s) -> {}/".format(len(regions), outdir))
    else:
        print("Single region mode: {} -> {}\n".format(args.region, args.out))
        make_combined_figure(
            hic_a_path=args.contact_a, hic_b_path=args.contact_b,
            region_str=args.region, out_path=args.out, **common)


if __name__ == "__main__":
    main()
