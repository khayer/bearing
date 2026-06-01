from __future__ import annotations


def add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax, alpha=0.18, zorder=0):
    """Draw horizontal-track highlight spans along x-axis."""
    if not highlights:
        return
    for hs, he, hcol in highlights:
        x0 = genomic_to_ax(max(hs, region_start), region_start, region_end)
        x1 = genomic_to_ax(min(he, region_end), region_start, region_end)
        ax.axvspan(x0, x1, color=hcol, alpha=alpha, zorder=zorder)


def add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax, alpha=0.18, zorder=0):
    """Draw vertical-track highlight spans along y-axis."""
    if not highlights:
        return
    for hs, he, hcol in highlights:
        y0 = genomic_to_ax(max(hs, region_start), region_start, region_end)
        y1 = genomic_to_ax(min(he, region_end), region_start, region_end)
        ax.axhspan(y0, y1, color=hcol, alpha=alpha, zorder=zorder)
