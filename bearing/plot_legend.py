from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def draw_legend_panel(
    ax,
    categories,
    num_states,
    hic_cmap,
    hic_vmax=None,
    rgb_mode=False,
    label_a="Condition A",
    label_b="Condition B",
    rgb_palette="magenta-green",
):
    """Draw the legend panel for bearing_hic_plot outputs."""
    ax.set_axis_off()

    ax.add_patch(mpatches.FancyBboxPatch(
        (0.02, 0.02), 0.96, 0.96,
        boxstyle="round,pad=0.015",
        facecolor=(1, 1, 1, 0.90),
        edgecolor="#c7c7c7",
        linewidth=0.8,
        transform=ax.transAxes,
        zorder=0,
    ))

    colorbar_height = 0.0
    if rgb_mode:
        ax.text(0.5, 0.97, "RGB Hi-C legend",
                transform=ax.transAxes, fontsize=8,
                ha="center", va="top", color="#222222", fontweight="bold")

        cax_rgb = inset_axes(ax, width="96%", height="52%",
                             loc="upper center",
                 bbox_to_anchor=(0, -0.06, 1, 1),
                             bbox_transform=ax.transAxes,
                             borderpad=0)
        palette_stops = {
            "magenta-green": (
                np.array([1.0, 0.0, 1.0], dtype=np.float64),
                np.array([1.0, 1.0, 1.0], dtype=np.float64),
                np.array([0.0, 1.0, 0.0], dtype=np.float64),
            ),
            "red-green": (
                np.array([1.0, 0.0, 0.0], dtype=np.float64),
                np.array([1.0, 1.0, 0.0], dtype=np.float64),
                np.array([0.0, 1.0, 0.0], dtype=np.float64),
            ),
            "blue-red": (
                np.array([0.0, 0.0, 1.0], dtype=np.float64),
                np.array([1.0, 0.0, 1.0], dtype=np.float64),
                np.array([1.0, 0.0, 0.0], dtype=np.float64),
            ),
            "green-blue": (
                np.array([0.0, 1.0, 0.0], dtype=np.float64),
                np.array([0.0, 1.0, 1.0], dtype=np.float64),
                np.array([0.0, 0.0, 1.0], dtype=np.float64),
            ),
            "magenta-green-white": (
                np.array([1.0, 0.0, 1.0], dtype=np.float64),  # magenta (A only)
                np.array([0.2, 0.0, 0.3], dtype=np.float64),  # dark purple (both)
                np.array([0.0, 1.0, 0.0], dtype=np.float64),  # green (B only)
            ),
        }
        if rgb_palette not in palette_stops:
            raise ValueError(
                "Unknown RGB palette: "
                f"{rgb_palette}. Use one of: magenta-green, red-green, blue-red, green-blue, magenta-green-white"
            )
        c_left, c_both, c_right = palette_stops[rgb_palette]

        h, w = 256, 256
        x = np.linspace(0.0, 1.0, w, dtype=np.float64)
        y = np.linspace(0.0, 1.0, h, dtype=np.float64)
        xx, yy = np.meshgrid(x, y)

        # Triangle vertices in axes coordinates:
        # top-left (A high), top-right (B high), bottom-center (both high).
        # No-contact sits above the triangle.
        v_left = np.array([0.08, 0.82], dtype=np.float64)
        v_right = np.array([0.92, 0.82], dtype=np.float64)
        v_bottom = np.array([0.50, 0.12], dtype=np.float64)

        x1, y1 = v_left
        x2, y2 = v_right
        x3, y3 = v_bottom
        den = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)

        w_left = ((y2 - y3) * (xx - x3) + (x3 - x2) * (yy - y3)) / den
        w_right = ((y3 - y1) * (xx - x3) + (x1 - x3) * (yy - y3)) / den
        w_bottom = 1.0 - w_left - w_right
        inside = (w_left >= 0.0) & (w_right >= 0.0) & (w_bottom >= 0.0)

        if rgb_palette == "magenta-green-white":
            # Rasterize triangle interior using parameterization (u, v):
            #   - u in [0,1] moves along the top edge from left (0) -> right (1)
            #   - v in [0,1] moves from the top edge (v=0) down to the bottom apex (v=1)
            # Map intensities:
            #   intensity_A = v + (1-v) * left_scale
            #   intensity_B = v + (1-v) * right_scale
            # where left_scale = clip(1 - 2*u, 0,1) and right_scale = clip(2*u - 1,0,1)
            # This makes top-edge midpoint (u=0.5,v=0) -> (0,0), left->(1,0), right->(0,1), bottom->(1,1).
            c00 = np.array([1.0, 1.0, 1.0], dtype=np.float64)  # white (no contact)
            c10 = c_left  # magenta (A only)
            c01 = c_right  # green (B only)
            c11 = c_both  # dark purple (both)

            # Avoid division by zero when projecting to top-edge
            lr_sum = (w_left + w_right)
            u = np.where(lr_sum > 1e-12, w_right / lr_sum, 0.5)
            v = np.clip(w_bottom, 0.0, 1.0)

            left_scale = np.clip(1.0 - 2.0 * u, 0.0, 1.0)
            right_scale = np.clip(2.0 * u - 1.0, 0.0, 1.0)

            a = v + (1.0 - v) * left_scale
            b = v + (1.0 - v) * right_scale

            w00 = (1.0 - a) * (1.0 - b)
            w10 = a * (1.0 - b)
            w01 = (1.0 - a) * b
            w11 = a * b

            rgb_img = np.zeros((h, w, 4), dtype=np.float64)
            rgb_img[:, :, :3] = (
                w00[:, :, None] * c00 +
                w10[:, :, None] * c10 +
                w01[:, :, None] * c01 +
                w11[:, :, None] * c11
            )
            rgb_img[:, :, 3] = inside.astype(np.float64)
        else:
            # Standard palettes: piecewise hue with brightness modulation
            denom_lr = np.maximum(w_left + w_right, 1e-12)
            hue_frac = np.clip(w_right / denom_lr, 0.0, 1.0)
            # Flip vertically: brightest near the bottom apex (both), dimmest near the top edge.
            brightness = np.clip(0.22 + 0.78 * w_bottom, 0.0, 1.0)

            # Piecewise hue path: left -> both -> right, then modulate by vertical brightness.
            t_left = np.clip(hue_frac / 0.5, 0.0, 1.0)
            t_right = np.clip((hue_frac - 0.5) / 0.5, 0.0, 1.0)
            left_half = hue_frac <= 0.5

            base = np.zeros((h, w, 3), dtype=np.float64)
            base[left_half] = (
                (1.0 - t_left[left_half])[:, None] * c_left
                + t_left[left_half][:, None] * c_both
            )
            base[~left_half] = (
                (1.0 - t_right[~left_half])[:, None] * c_both
                + t_right[~left_half][:, None] * c_right
            )

            rgb_img = np.zeros((h, w, 4), dtype=np.float64)
            rgb_img[:, :, :3] = brightness[:, :, None] * base
            rgb_img[:, :, 3] = inside.astype(np.float64)

        cax_rgb.imshow(
            rgb_img,
            aspect="auto",
            origin="lower",
            interpolation="bicubic",
            extent=[0, 1, 0, 1],
        )

        # Triangle border and labels.
        border = np.array([v_left, v_right, v_bottom, v_left])
        cax_rgb.plot(border[:, 0], border[:, 1], color="#6f6f6f", linewidth=1.0)
        cax_rgb.set_xlim(0, 1)
        cax_rgb.set_ylim(0, 1)
        cax_rgb.set_xticks([])
        cax_rgb.set_yticks([])
        for spine in cax_rgb.spines.values():
            spine.set_visible(False)

        cax_rgb.text(v_bottom[0], v_bottom[1] - 0.02, "Both",
                     fontsize=6.6, color="#222222", ha="center", va="top", fontweight="bold")
        cax_rgb.text(v_left[0] - 0.01, v_left[1], label_a,
                     fontsize=6.5, color="#222222", ha="right", va="center", fontweight="bold")
        cax_rgb.text(v_right[0] + 0.01, v_right[1], label_b,
                     fontsize=6.5, color="#222222", ha="left", va="center", fontweight="bold")
        
        if rgb_palette == "magenta-green-white":
            # For magenta-green-white, "No contact" label sits above the triangle
            # where white region naturally appears
            cax_rgb.text(0.5, 0.94, "No contact",
                         fontsize=6.2, color="#222222", ha="center", va="bottom", fontweight="bold")
        else:
            # For standard palettes, place label above
            cax_rgb.text(0.5, 0.93, "No contact",
                         fontsize=6.2, color="#444444", ha="center", va="bottom")

        arrow_color = "#444444" if rgb_palette == "magenta-green-white" else "#8a8a8a"
        contact_text_color = "#444444" if rgb_palette == "magenta-green-white" else "#7a7a7a"

        cax_rgb.annotate(
            "",
            xy=(0.5, 0.18),
            xytext=(0.5, 0.72),
            arrowprops=dict(arrowstyle="-|>", color=arrow_color, lw=0.9),
        )
        cax_rgb.text(0.535, 0.43, "contact\nintensity",
                     fontsize=5.8, color=contact_text_color, ha="left", va="center")

        colorbar_height = 0.76
    elif hic_vmax is not None:
        cax = inset_axes(ax, width="80%", height="8%",
                         loc="upper center",
                         bbox_to_anchor=(0, 0.08, 1, 1),
                         bbox_transform=ax.transAxes,
                         borderpad=0)
        norm = mpl.colors.Normalize(vmin=0, vmax=hic_vmax)
        cb = mpl.colorbar.ColorbarBase(cax, cmap=hic_cmap, norm=norm,
                                       orientation="horizontal")
        cb.ax.tick_params(labelsize=6, length=2, width=0.5)
        cb.set_ticks([0, hic_vmax])
        cb.set_ticklabels(["0", f"{hic_vmax:.1f}"])
        cb.ax.xaxis.set_label_position("top")
        cb.ax.xaxis.tick_top()
        ax.text(0.5, 0.97, "Hi-C contacts (log1p)",
                transform=ax.transAxes, fontsize=8,
                ha="center", va="top", color="#222222", fontweight="bold")
        colorbar_height = 0.30

    n = num_states
    cols = 2
    rows = (n + cols - 1) // cols
    swatch_top = 1.0 - colorbar_height - 0.08
    row_step = swatch_top / (rows + 1)
    for i in range(n):
        cat_name, color = categories[i]
        if color.lower() == "#ffffff":
            color = "#d0d0d0"
        row = i % rows
        col = i // rows
        x = 0.05 + col * 0.5
        y = swatch_top - (row + 1) * row_step
        ax.add_patch(mpatches.Rectangle(
            (x, y - 0.024), 0.075, 0.048,
            facecolor=color, edgecolor="#555555", linewidth=0.3,
            transform=ax.transAxes,
            clip_on=False, zorder=2,
        ))
        ax.text(x + 0.095, y, cat_name,
                transform=ax.transAxes,
                fontsize=7, va="center", color="#222222")
    signal_title_y = swatch_top + (0.01 if rgb_mode else 0.04)
    ax.text(0.5, signal_title_y, "Signal Legend",
            transform=ax.transAxes,
            fontsize=8, ha="center", va="bottom",
            fontweight="bold", color="#222222")
