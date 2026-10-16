from __future__ import annotations


def _clipped_span(hs, he, region_start, region_end, genomic_to_ax):
    """Map a genomic interval, clipped to the region, to axis coordinates."""
    lo = genomic_to_ax(max(hs, region_start), region_start, region_end)
    hi = genomic_to_ax(min(he, region_end), region_start, region_end)
    return lo, hi


def _draw_highlights(span_fn, highlights, region_start, region_end, genomic_to_ax, alpha, zorder):
    for hs, he, hcol in highlights or ():
        lo, hi = _clipped_span(hs, he, region_start, region_end, genomic_to_ax)
        span_fn(lo, hi, color=hcol, alpha=alpha, zorder=zorder)


def add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax, alpha=0.18, zorder=0):
    """Draw horizontal-track highlight spans along x-axis."""
    _draw_highlights(ax.axvspan, highlights, region_start, region_end, genomic_to_ax, alpha, zorder)


def add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax, alpha=0.18, zorder=0):
    """Draw vertical-track highlight spans along y-axis."""
    _draw_highlights(ax.axhspan, highlights, region_start, region_end, genomic_to_ax, alpha, zorder)
