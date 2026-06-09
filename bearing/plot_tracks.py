from __future__ import annotations

import numpy as np
import matplotlib.patches as mpatches
from matplotlib.patches import Arc, FancyArrowPatch

from .track_primitives import add_horizontal_highlights, add_vertical_highlights


def genomic_to_ax(pos, region_start, region_end):
    return (pos - region_start) / (region_end - region_start)


def add_vertical_tick_bars(ax, region_start, region_end):
    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        if 4 <= span / interval <= 12:
            break

    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)

    ticks_ax = [genomic_to_ax(t, region_start, region_end) for t in ticks_genomic]
    ax.set_yticks(ticks_ax)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", left=True, right=False, length=2.5, width=0.6, colors="#777777")


def draw_epilogos_horizontal(ax, positions, score_mat, num_states,
                             region_start, region_end,
                             categories, highlights=None, label="",
                             y_max=None):
    ax.set_xlim(0, 1)
    _y_max = y_max if y_max is not None else (score_mat.max() * 1.15 if score_mat.size > 0 else 1)
    ax.set_ylim(0, _y_max)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        x0 = genomic_to_ax(bs, region_start, region_end)
        x1 = genomic_to_ax(be, region_start, region_end)
        width = max(x1 - x0, 0.0002)
        row = score_mat[i]
        bottom = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val <= 0:
                continue
            color = categories[si][1] if si < len(categories) else "#cccccc"
            ax.bar(x0, val, width=width, bottom=bottom,
                   color=color, align="edge", linewidth=0, zorder=2)
            bottom += val

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_epilogos_vertical(ax, positions, score_mat, num_states,
                           region_start, region_end,
                           categories, highlights=None, label="",
                           y_max=None):
    ax.set_ylim(1, 0)
    _y_max = y_max if y_max is not None else (score_mat.max() * 1.15 if score_mat.size > 0 else 1)
    ax.set_xlim(_y_max, 0)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        y0 = genomic_to_ax(bs, region_start, region_end)
        y1 = genomic_to_ax(be, region_start, region_end)
        height = max(y1 - y0, 0.0002)
        row = score_mat[i]
        left = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val <= 0:
                continue
            color = categories[si][1] if si < len(categories) else "#cccccc"
            ax.barh(y0, val, height=height, left=left,
                    color=color, align="edge", linewidth=0, zorder=2)
            left += val

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_diff_horizontal(ax, positions, score_mat, num_states,
                         region_start, region_end,
                         categories, highlights=None, label="",
                         diff_max=None):
    if score_mat.size == 0:
        ax.set_axis_off()
        return

    if diff_max is not None:
        abs_max = diff_max
    else:
        abs_max = np.abs(score_mat).max()
        if abs_max == 0:
            abs_max = 1.0
    ax.set_xlim(0, 1)
    ax.set_ylim(abs_max * 1.15, -abs_max * 1.15)
    ax.axhline(0, color="#888888", linewidth=0.6, zorder=1)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        x0 = genomic_to_ax(bs, region_start, region_end)
        x1 = genomic_to_ax(be, region_start, region_end)
        width = max(x1 - x0, 0.0002)
        row = score_mat[i]
        pos_bottom = 0.0
        neg_bottom = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val == 0:
                continue
            color = categories[si][1] if si < len(categories) else "#cccccc"
            if val > 0:
                ax.bar(x0, val, width=width, bottom=pos_bottom,
                       color=color, align="edge", linewidth=0, zorder=2)
                pos_bottom += val
            else:
                ax.bar(x0, val, width=width, bottom=neg_bottom,
                       color=color, align="edge", linewidth=0, zorder=2)
                neg_bottom += val

    # Note: subtle baseline fill removed — the zero line is sufficient

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_diff_vertical(ax, positions, score_mat, num_states,
                       region_start, region_end,
                       categories, highlights=None, label="",
                       diff_max=None):
    if score_mat.size == 0:
        ax.set_axis_off()
        return

    if diff_max is not None:
        abs_max = diff_max
    else:
        abs_max = np.abs(score_mat).max()
        if abs_max == 0:
            abs_max = 1.0

    ax.set_ylim(1, 0)
    ax.set_xlim(abs_max * 1.15, -abs_max * 1.15)
    ax.axvline(0, color="#888888", linewidth=0.6, zorder=1)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        y0 = genomic_to_ax(bs, region_start, region_end)
        y1 = genomic_to_ax(be, region_start, region_end)
        height = max(y1 - y0, 0.0002)
        row = score_mat[i]
        pos_left = 0.0
        neg_left = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val == 0:
                continue
            color = categories[si][1] if si < len(categories) else "#cccccc"
            if val > 0:
                ax.barh(y0, val, height=height, left=pos_left,
                        color=color, align="edge", linewidth=0, zorder=2)
                pos_left += val
            else:
                ax.barh(y0, val, height=height, left=neg_left,
                        color=color, align="edge", linewidth=0, zorder=2)
                neg_left += val

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    ax.spines["bottom"].set_visible(False)
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_gene_track(ax, genes, region_start, region_end,
                    highlights=None, label="Genes"):
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 1.5)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    used = []
    for gs, ge, name, strand in genes:
        x0 = genomic_to_ax(max(gs, region_start), region_start, region_end)
        x1 = genomic_to_ax(min(ge, region_end), region_start, region_end)
        if x1 <= x0:
            continue

        yl = 0.2
        for prev_x0, prev_x1, prev_y in used:
            if not (x1 < prev_x0 or x0 > prev_x1) and prev_y == 0.2:
                yl = 0.8
                break
        used.append((x0, x1, yl))

        ax.plot([x0, x1], [yl, yl], color="#555555", lw=2.5, solid_capstyle="butt", zorder=2)
        arrow_x = x1 if strand == "+" else x0
        dx = 0.01 if strand == "+" else -0.01
        ax.annotate("", xy=(arrow_x + dx, yl), xytext=(arrow_x, yl),
                    arrowprops=dict(arrowstyle="->", color="#555555", lw=1.0, mutation_scale=6), zorder=3)

        if (x1 - x0) >= 0.004:
            mid_x = 0.5 * (x0 + x1)
            ax.text(mid_x, yl + 0.15, name, fontsize=6, ha="center", va="bottom",
                    color="#333333", style="italic", rotation=0, zorder=3)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


def draw_gene_track_vertical(ax, genes, region_start, region_end,
                             highlights=None, label="Genes"):
    ax.set_ylim(1, 0)
    x_levels = [0.64, 0.79, 0.92]
    ax.set_xlim(0.56, 1.0)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    used = []
    for gs, ge, name, strand in genes:
        y0 = genomic_to_ax(max(gs, region_start), region_start, region_end)
        y1 = genomic_to_ax(min(ge, region_end), region_start, region_end)
        if y1 <= y0:
            continue

        xl = x_levels[0]
        for lvl in x_levels:
            overlap = any(not (y1 < py0 or y0 > py1) and px == lvl for py0, py1, px in used)
            if not overlap:
                xl = lvl
                break
        used.append((y0, y1, xl))

        ax.plot([xl, xl], [y0, y1], color="#555555", lw=2.0, solid_capstyle="butt", zorder=2)
        arrow_y = y1 if strand == "+" else y0
        dy = 0.01 if strand == "+" else -0.01
        ax.annotate("", xy=(xl, arrow_y + dy), xytext=(xl, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#555555", lw=1.0, mutation_scale=6), zorder=3)

        mid_y = 0.5 * (y0 + y1)
        if (y1 - y0) >= 0.004:
            ax.text(max(0.57, xl - 0.025), mid_y, name,
                    fontsize=4, ha="right", va="center",
                    color="#333333", style="italic", rotation=90, zorder=3)

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_genomic_axis(ax, region_start, region_end, chrom):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        n_ticks = span / interval
        if 4 <= n_ticks <= 12:
            break

    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)

    ticks_ax = [genomic_to_ax(t, region_start, region_end) for t in ticks_genomic]
    ax.set_xticks(ticks_ax)
    ax.set_xticklabels([f"{t/1e6:.2f}" for t in ticks_genomic], fontsize=6, rotation=45, ha="right")
    ax.set_yticks([])
    ax.set_xlabel(f"{chrom} (Mb)", fontsize=7, labelpad=2)
    ax.tick_params(axis="x", length=3, width=0.6, labelsize=6)
    ax.spines["bottom"].set_linewidth(0.6)


def draw_genomic_axis_vertical(ax, region_start, region_end, chrom):
    ax.set_ylim(1, 0)
    ax.set_xlim(0, 1)
    for sp in ["top", "right", "bottom"]:
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)

    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        if 4 <= span / interval <= 12:
            break

    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)

    ticks_ax = [genomic_to_ax(t, region_start, region_end) for t in ticks_genomic]
    ax.set_yticks(ticks_ax)
    ax.set_yticklabels([f"{t/1e6:.2f}" for t in ticks_genomic], fontsize=5, rotation=45, ha="left")
    ax.set_xticks([])
    ax.set_ylabel("Mb", fontsize=6, labelpad=2)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", length=3, width=0.6, labelsize=5, pad=1)


def draw_pval_horizontal(ax, positions, values, region_start, region_end,
                         highlights=None, label="", color="#a43cca",
                         y_max=None, cutoff_value=None):
    ax.set_xlim(0, 1)
    if y_max is None:
        vmax = float(values.max()) * 1.1 if values.size > 0 and values.max() > 0 else 1.0
    else:
        vmax = y_max
    ax.set_ylim(vmax, 0)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    xs = np.array([genomic_to_ax(p, region_start, region_end) for p in positions])
    ax.fill_between(xs, 0, values, color=color, alpha=0.65, linewidth=0, zorder=2)
    ax.plot(xs, values, color=color, linewidth=0.6, zorder=3)

    if cutoff_value is not None and cutoff_value >= 0:
        ax.axhline(cutoff_value, color="#666666", linestyle="--", linewidth=0.8, zorder=4)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_pval_vertical(ax, positions, values, region_start, region_end,
                       highlights=None, label="", color="#a43cca",
                       x_max=None, cutoff_value=None):
    ax.set_ylim(1, 0)
    if x_max is None:
        vmax = float(values.max()) * 1.1 if values.size > 0 and values.max() > 0 else 1.0
    else:
        vmax = x_max
    ax.set_xlim(0, vmax)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    ys = np.array([genomic_to_ax(p, region_start, region_end) for p in positions])
    ax.fill_betweenx(ys, 0, values, color=color, alpha=0.65, linewidth=0, zorder=2)
    ax.plot(values, ys, color=color, linewidth=0.6, zorder=3)

    if cutoff_value is not None and cutoff_value >= 0:
        ax.axvline(cutoff_value, color="#666666", linestyle="--", linewidth=0.8, zorder=4)

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_bed_track(
    ax,
    features,
    region_start,
    region_end,
    style="itemRgb",
    label="",
    label_threshold_kb=200,
    label_density_threshold=30,
):
    """
    Draw a BED track row. styles: 'cbe' (vertical narrow bars) or 'itemRgb' (horizontal boxes).
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    add_horizontal_highlights(ax, None, region_start, region_end, genomic_to_ax)

    span = region_end - region_start
    width_kb = span / 1000.0
    suppress_labels = (width_kb > label_threshold_kb and len(features) > label_density_threshold)

    if style == "cbe":
        # vertical thin bars at feature start with directional strand indicators
        for feat in features:
            x = genomic_to_ax(max(feat["start"], region_start), region_start, region_end)
            strand = feat.get("strand")
            if strand is None or strand == ".":
                color = "#8b0000"
                is_plus = True
            else:
                is_plus = strand == "+"
                color = "#8b0000" if is_plus else "#1f3d99"
            ax.vlines(x, 0.1, 0.9, color=color, linewidth=0.6, zorder=2)
            
            # Draw directional triangle: right (>) for + strand, left (<) for - strand
            triangle_width = 0.02  # triangle width in axes coordinates
            if is_plus:
                # Right-pointing triangle (>): tip points right
                triangle = mpatches.Polygon(
                    [(x, 0.5), (x + triangle_width, 0.65), (x + triangle_width, 0.35)],
                    facecolor=color, edgecolor="none", zorder=2
                )
            else:
                # Left-pointing triangle (<): tip points left
                triangle = mpatches.Polygon(
                    [(x, 0.5), (x - triangle_width, 0.65), (x - triangle_width, 0.35)],
                    facecolor=color, edgecolor="none", zorder=2
                )
            ax.add_patch(triangle)
        # labels
        if not suppress_labels:
            # greedy placement left-to-right
            last_x = -1.0
            for feat in sorted(features, key=lambda f: f["start"]):
                x = genomic_to_ax(max(feat["start"], region_start), region_start, region_end)
                if last_x >= 0 and x - last_x < 0.02:
                    continue
                ax.text(x, 0.95, feat.get("name", ""), fontsize=6, ha="center", va="bottom", color="#222222", zorder=3)
                last_x = x
        ax.set_ylabel(label, fontsize=7, labelpad=3)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    else:
        # itemRgb: horizontal bars spanning feature extent
        for feat in features:
            x0 = genomic_to_ax(max(feat["start"], region_start), region_start, region_end)
            x1 = genomic_to_ax(min(feat["end"], region_end), region_start, region_end)
            if x1 <= x0:
                continue
            rgb = feat.get("item_rgb")
            if rgb is None:
                face = "#888888"
            else:
                face = f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
            height = 0.6
            rect_x = x0
            rect_w = max(x1 - x0, 0.0005)
            ax.add_patch(mpatches.Rectangle((rect_x, 0.2), rect_w, 0.6, facecolor=face, edgecolor="none", transform=ax.transAxes, zorder=2))

        if not suppress_labels:
            last_x = -1.0
            for feat in sorted(features, key=lambda f: f["start"]):
                x0 = genomic_to_ax(max(feat["start"], region_start), region_start, region_end)
                x1 = genomic_to_ax(min(feat["end"], region_end), region_start, region_end)
                mid = 0.5 * (x0 + x1)
                if last_x >= 0 and mid - last_x < 0.02:
                    continue
                ax.text(mid, 0.5, feat.get("name", ""), fontsize=6, ha="center", va="center", color="#111111", zorder=3)
                last_x = mid

        ax.set_ylabel(label, fontsize=7, labelpad=3)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_visible(False)


def _build_diff_pval_lollipop_specs(positions, values,
                                    diff_score_positions, diff_score_matrix,
                                    categories, cutoff_value,
                                    tolerance_bp=100,
                                    neutral_color="#9aa0a6"):
    """Return per-bin lollipop specs for significant differential p-value bins."""
    if (
        positions is None or values is None or diff_score_positions is None
        or diff_score_matrix is None or categories is None or cutoff_value is None
    ):
        return []

    positions = np.asarray(positions)
    values = np.asarray(values)
    diff_score_positions = np.asarray(diff_score_positions)
    diff_score_matrix = np.asarray(diff_score_matrix)

    if positions.size == 0 or values.size == 0:
        return []
    if diff_score_positions.size == 0 or diff_score_matrix.size == 0:
        return []

    n_ref = min(diff_score_positions.shape[0], diff_score_matrix.shape[0])
    diff_score_positions = diff_score_positions[:n_ref]
    diff_score_matrix = diff_score_matrix[:n_ref]

    if diff_score_positions.ndim != 1 or diff_score_matrix.ndim != 2:
        return []
    if diff_score_matrix.shape[0] != diff_score_positions.shape[0]:
        return []

    if diff_score_positions.shape[0] > 1 and np.any(np.diff(diff_score_positions) < 0):
        order = np.argsort(diff_score_positions)
        diff_score_positions = diff_score_positions[order]
        diff_score_matrix = diff_score_matrix[order, :]

    specs = []
    for pos, value in zip(positions, values):
        if abs(float(value)) < float(cutoff_value):
            continue

        insert_at = int(np.searchsorted(diff_score_positions, pos))
        candidates = []
        if insert_at < diff_score_positions.shape[0]:
            candidates.append(insert_at)
        if insert_at > 0:
            candidates.append(insert_at - 1)

        best_idx = None
        best_dist = None
        for cand in candidates:
            dist = abs(float(diff_score_positions[cand]) - float(pos))
            if best_dist is None or dist < best_dist:
                best_idx = cand
                best_dist = dist

        color = neutral_color
        if best_idx is not None and best_dist is not None and best_dist <= float(tolerance_bp):
            sgn = 1.0 if float(value) > 0 else -1.0
            per_track_signed = diff_score_matrix[best_idx, :] * sgn
            dominant_state = int(np.argmax(per_track_signed))
            if 0 <= dominant_state < len(categories):
                color = categories[dominant_state][1]

        specs.append((float(pos), float(value), color))

    return specs


def draw_pval_diff_horizontal(ax, positions, values, region_start, region_end,
                              highlights=None, label="",
                              pos_color="#0c4a6e", neg_color="#d97706",
                              y_max=None, cutoff_value=None,
                              diff_score_positions=None, diff_score_matrix=None,
                              categories=None):
    if values is None or len(values) == 0:
        ax.set_axis_off()
        return

    ax.set_xlim(0, 1)
    abs_max = y_max if y_max is not None else float(np.abs(values).max())
    if abs_max <= 0:
        abs_max = 1.0
    ax.set_ylim(abs_max * 1.15, -abs_max * 1.15)
    ax.axhline(0, color="#888888", linewidth=0.6, zorder=1)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    xs = np.array([genomic_to_ax(p, region_start, region_end) for p in positions])
    pos = np.clip(values, 0, None)
    neg = np.clip(values, None, 0)
    ax.fill_between(xs, 0, pos, color=pos_color, alpha=0.20, linewidth=0, zorder=2)
    ax.fill_between(xs, 0, neg, color=neg_color, alpha=0.20, linewidth=0, zorder=2)
    ax.plot(xs, values, color="#333333", linewidth=0.55, zorder=3)

    if cutoff_value is not None and cutoff_value >= 0:
        ax.axhline(cutoff_value, color="#444444", linewidth=0.8, zorder=4)
        ax.axhline(-cutoff_value, color="#444444", linewidth=0.8, zorder=4)
        cutoff_p = 10 ** (-cutoff_value)
        ax.text(
            0.995, cutoff_value + 0.05, f"p<{cutoff_p:.2g}",
            ha="right", va="bottom", fontsize=6, color="#444",
            transform=ax.get_yaxis_transform(),
        )
        ax.text(
            0.995, -cutoff_value - 0.05, f"p<{cutoff_p:.2g}",
            ha="right", va="top", fontsize=6, color="#444",
            transform=ax.get_yaxis_transform(),
        )

    for x, value, color in _build_diff_pval_lollipop_specs(
        positions, values,
        diff_score_positions, diff_score_matrix,
        categories, cutoff_value,
    ):
        ax.vlines(x, 0, value, color=color, linewidth=0.95, alpha=1.0, zorder=5)
        ax.scatter(
            [x], [value], s=9, c=[color], edgecolors="#222", linewidths=0.3,
            zorder=6,
        )

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_pval_manhattan_horizontal(ax, positions, values, region_start, region_end,
                                  highlights=None, label="",
                                  score_positions=None, score_matrix=None,
                                  categories=None, y_max=None,
                                  cutoff_value=None, kl_scores=None):
    """
    Draw a Manhattan-style plot for p-values with dots colored by dominant category.
    
    Can color dots in two ways:
    1. If kl_scores is provided: directly use the KL values at each position to determine dominant category
    2. Otherwise: use score_matrix and score_positions to look up dominant category
    
    Args:
        kl_scores: numpy array of shape (n_bins, n_categories) with per-category scores.
                   If provided, category colors are determined from the dominant category
                   at each position without needing position matching.
    """
    if values is None or len(values) == 0:
        ax.set_axis_off()
        return

    ax.set_xlim(0, 1)
    abs_max = y_max if y_max is not None else float(np.abs(values).max())
    if abs_max <= 0:
        abs_max = 1.0
    ax.set_ylim(abs_max * 1.15, -abs_max * 1.15)
    ax.axhline(0, color="#888888", linewidth=0.6, zorder=1)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    xs = np.array([genomic_to_ax(p, region_start, region_end) for p in positions])

    # Determine color for each position based on dominant category
    # in the direction of the per-bin differential value.
    # This matches the semantics in _build_diff_pval_lollipop_specs:
    # we color by "which track drove the differential, in the sign of the differential."
    values_arr = np.asarray(values, dtype=float)
    colors = []
    if kl_scores is not None and kl_scores.size > 0 and len(kl_scores) == len(positions):
        # Use KL scores directly; argmax in the direction of the differential value.
        for i, row in enumerate(kl_scores):
            color = "#9aa0a6"  # default grey
            if len(row) > 0 and categories is not None:
                sgn = 1.0 if values_arr[i] > 0 else -1.0
                per_track_signed = np.asarray(row, dtype=float) * sgn
                dom_cat_idx = int(np.argmax(per_track_signed))
                if 0 <= dom_cat_idx < len(categories):
                    cat = categories[dom_cat_idx]
                    color = cat[1] if isinstance(cat, (tuple, list)) else cat.color
            colors.append(color)
    elif categories is not None and score_matrix is not None and score_positions is not None:
        # Fallback path: nearest-neighbor lookup in score_matrix.
        # Same direction-aware argmax for consistency.
        score_pos_scalar = []
        for sp in score_positions:
            if isinstance(sp, (tuple, list)):
                score_pos_scalar.append(sp[0])
            else:
                score_pos_scalar.append(sp)

        score_pos_array = np.array(score_pos_scalar, dtype=float)

        for i, pos in enumerate(positions):
            color = "#9aa0a6"
            distances = np.abs(score_pos_array - float(pos))
            if len(distances) > 0:
                idx = np.argmin(distances)
                if 0 <= idx < len(score_matrix):
                    row = score_matrix[idx]
                    if len(row) > 0:
                        sgn = 1.0 if values_arr[i] > 0 else -1.0
                        per_track_signed = np.asarray(row, dtype=float) * sgn
                        dom_cat_idx = int(np.argmax(per_track_signed))
                        if 0 <= dom_cat_idx < len(categories):
                            cat = categories[dom_cat_idx]
                            color = cat[1] if isinstance(cat, (tuple, list)) else cat.color
            colors.append(color)
    else:
        colors = ["#9aa0a6"] * len(positions)

    # Draw stems and colored dots; emphasize values passing the cutoff
    sig_mask = None
    if cutoff_value is not None:
        sig_mask = np.abs(values) >= cutoff_value
    else:
        sig_mask = np.array([False] * len(values))

    # Calculate max absolute value for dynamic sizing
    abs_max = np.max(np.abs(values)) if len(values) > 0 else 1.0

    for i, (x, value, color) in enumerate(zip(xs, values, colors)):
        lw = 1.2 if sig_mask[i] else 0.9
        alpha = 0.95 if sig_mask[i] else 0.8
        ax.vlines(x, 0, value, color=color, linewidth=lw, alpha=alpha, zorder=4)
        s = max(4.0, 14.0 * (abs(value) / abs_max))
        ax.scatter([x], [value], s=s, c=[color], edgecolors="#222", linewidths=0.3, zorder=5)

    # Draw cutoff lines and labels for significance if requested
    if cutoff_value is not None and cutoff_value >= 0:
        ax.axhline(cutoff_value, color="#444444", linewidth=1.0, linestyle="--", zorder=3)
        ax.axhline(-cutoff_value, color="#444444", linewidth=1.0, linestyle="--", zorder=3)
        cutoff_p = 10 ** (-cutoff_value)
        ax.text(0.995, cutoff_value + 0.05, f"p<{cutoff_p:.2g}", ha="right", va="bottom", fontsize=6, color="#444",
                transform=ax.get_yaxis_transform())
        ax.text(0.995, -cutoff_value - 0.05, f"p<{cutoff_p:.2g}", ha="right", va="top", fontsize=6, color="#444",
                transform=ax.get_yaxis_transform())

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_pval_diff_vertical(ax, positions, values, region_start, region_end,
                            highlights=None, label="",
                            pos_color="#0c4a6e", neg_color="#d97706",
                            x_max=None, cutoff_value=None,
                            diff_score_positions=None, diff_score_matrix=None,
                            categories=None):
    if values is None or len(values) == 0:
        ax.set_axis_off()
        return

    ax.set_ylim(1, 0)
    abs_max = x_max if x_max is not None else float(np.abs(values).max())
    if abs_max <= 0:
        abs_max = 1.0
    ax.set_xlim(abs_max * 1.15, -abs_max * 1.15)
    ax.axvline(0, color="#888888", linewidth=0.6, zorder=1)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    ys = np.array([genomic_to_ax(p, region_start, region_end) for p in positions])
    pos = np.clip(values, 0, None)
    neg = np.clip(values, None, 0)
    ax.fill_betweenx(ys, 0, pos, color=pos_color, alpha=0.20, linewidth=0, zorder=2)
    ax.fill_betweenx(ys, 0, neg, color=neg_color, alpha=0.20, linewidth=0, zorder=2)
    ax.plot(values, ys, color="#333333", linewidth=0.55, zorder=3)

    if cutoff_value is not None and cutoff_value >= 0:
        ax.axvline(cutoff_value, color="#444444", linewidth=0.8, zorder=3)
        ax.axvline(-cutoff_value, color="#444444", linewidth=0.8, zorder=3)
        cutoff_p = 10 ** (-cutoff_value)
        ax.text(
            cutoff_value + 0.05, 0.995, f"p<{cutoff_p:.2g}",
            ha="left", va="top", fontsize=6, color="#444",
            transform=ax.get_xaxis_transform(),
        )
        ax.text(
            -cutoff_value - 0.05, 0.995, f"p<{cutoff_p:.2g}",
            ha="right", va="top", fontsize=6, color="#444",
            transform=ax.get_xaxis_transform(),
        )

    for x, value, color in _build_diff_pval_lollipop_specs(
        positions, values,
        diff_score_positions, diff_score_matrix,
        categories, cutoff_value,
    ):
        y = genomic_to_ax(x, region_start, region_end)
        ax.hlines(y, 0, value, color=color, linewidth=0.95, alpha=1.0, zorder=5)
        ax.scatter(
            [value], [y], s=9, c=[color], edgecolors="#222", linewidths=0.3,
            zorder=6,
        )

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_loops_horizontal(ax, loops, region_start, region_end,
                          highlights=None, label="", color="#4c78a8",
                          anchor_color="#4c78a8", anchor_height=0.10,
                          arc_alpha=0.45):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    if loops:
        max_dist = max(abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1))) for s1, e1, s2, e2, _ in loops)
        max_dist = max(max_dist, 1.0)
        max_score = max(lp[4] for lp in loops)
        max_score = max(max_score, 1e-9)

        for s1, e1, s2, e2, score in loops:
            a1s = genomic_to_ax(max(s1, region_start), region_start, region_end)
            a1e = genomic_to_ax(min(e1, region_end), region_start, region_end)
            a2s = genomic_to_ax(max(s2, region_start), region_start, region_end)
            a2e = genomic_to_ax(min(e2, region_end), region_start, region_end)
            if a1e <= 0 or a2e <= 0 or a1s >= 1 or a2s >= 1:
                continue

            # Short anchor boxes pinned to the baseline (default 10% height) so
            # they don't obscure the arcs; semi-transparent so overlaps read.
            ax.axvspan(a1s, a1e, ymin=0.0, ymax=anchor_height,
                       color=anchor_color, alpha=0.6, zorder=1)
            ax.axvspan(a2s, a2e, ymin=0.0, ymax=anchor_height,
                       color=anchor_color, alpha=0.6, zorder=1)

            m1 = genomic_to_ax(0.5 * (s1 + e1), region_start, region_end)
            m2 = genomic_to_ax(0.5 * (s2 + e2), region_start, region_end)
            left, right = sorted([m1, m2])
            width = right - left
            if width <= 0:
                continue
            dist_bp = abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1)))
            # Arc is the upper half of an ellipse centered at (xc, 0); its peak
            # sits at height/2. Use sqrt(distance) scaling so short-range loops
            # lift off the baseline and arcs spread vertically (peak ~0.30..0.95
            # of the track) instead of clustering flat near the anchors.
            norm = min(1.0, (dist_bp / max_dist) ** 0.5)
            height = 0.60 + 1.30 * norm

            # Constant low alpha on the arc so overlapping loops accumulate
            # visually (density) rather than the strongest one painting over.
            arc = Arc(((left + right) / 2.0, 0.0), width=width, height=height,
                      angle=0, theta1=0, theta2=180, lw=1.0,
                      color=color, alpha=arc_alpha, zorder=3)
            ax.add_patch(arc)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


def draw_loops_vertical(ax, loops, region_start, region_end,
                        highlights=None, label="", color="#e05c3a",
                        anchor_color="#e05c3a"):
    ax.set_ylim(1, 0)
    ax.set_xlim(0, 1)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    if loops:
        max_dist = max(abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1))) for s1, e1, s2, e2, _ in loops)
        max_dist = max(max_dist, 1.0)
        max_score = max(lp[4] for lp in loops)
        max_score = max(max_score, 1e-9)

        for s1, e1, s2, e2, score in loops:
            a1s = genomic_to_ax(max(s1, region_start), region_start, region_end)
            a1e = genomic_to_ax(min(e1, region_end), region_start, region_end)
            a2s = genomic_to_ax(max(s2, region_start), region_start, region_end)
            a2e = genomic_to_ax(min(e2, region_end), region_start, region_end)
            if a1e <= 0 or a2e <= 0 or a1s >= 1 or a2s >= 1:
                continue

            alpha = 0.35 + 0.55 * min(1.0, score / max_score)
            ax.axhspan(a1s, a1e, color=anchor_color, alpha=alpha, zorder=1)
            ax.axhspan(a2s, a2e, color=anchor_color, alpha=alpha, zorder=1)

            m1 = genomic_to_ax(0.5 * (s1 + e1), region_start, region_end)
            m2 = genomic_to_ax(0.5 * (s2 + e2), region_start, region_end)
            y1, y2 = sorted([m1, m2])
            if abs(y2 - y1) <= 0:
                continue
            dist_bp = abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1)))
            rad = 0.12 + 0.30 * (dist_bp / max_dist)

            arc = FancyArrowPatch((0.02, y1), (0.02, y2),
                                  connectionstyle=f"arc3,rad={rad}",
                                  arrowstyle="-", linewidth=0.9,
                                  color=color, alpha=alpha, zorder=3)
            ax.add_patch(arc)

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)
    ax.set_xticks([])
