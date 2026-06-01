#!/usr/bin/env python3
"""
bearing_tad_compare.py

Compare TAD boundaries across multiple conditions and visualize boundary
conservation, shifts, and gain/loss patterns. Designed to test whether
specific perturbations (e.g., enhancer deletion) alter TAD structure at
defined loci.

Input: one TAD BED file per condition.
Output:
  - A per-locus stacked-TAD figure (one row per condition)
  - Optional insulation overlay
  - A boundary-shift table identifying which boundaries change

Boundary classification:
  - CONSERVED: present in all conditions within +/- tolerance bp
  - SHIFTED:   present in all conditions but with boundary movement
  - CONDITION-SPECIFIC: absent in one or more conditions

Usage:
  python bearing_tad_compare.py \\
      --tads DN:tads/DN_tads.bed DP:tads/DP_tads.bed EbKO:tads/EbKO_tads.bed \\
             ProB:tads/ProB_tads.bed S3T3:tads/S3T3_tads.bed \\
      --insulation DN:tads/DN_insulation.bm DP:tads/DP_insulation.bm \\
                   EbKO:tads/EbKO_insulation.bm \\
      --region chr6:40000000-47400000 \\
      --annotations tcrb_wide_annotations.bed \\
      --tolerance 50000 \\
      --out tcrb_wide_tad_compare

Outputs:
  <out>.pdf            Multi-panel stacked figure
  <out>.boundaries.tsv Boundary calls per condition + classification
  <out>.shifts.tsv     Boundary-shift events between condition pairs
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


def parse_region_str(s):
    m = re.match(r"^([^:]+):([0-9,]+)-([0-9,]+)$", s)
    if not m:
        raise argparse.ArgumentTypeError(f"Cannot parse region: {s}")
    return m.group(1), int(m.group(2).replace(",", "")), \
        int(m.group(3).replace(",", ""))


def parse_labeled_path(s):
    """Parse 'LABEL:PATH' format."""
    if ":" not in s:
        raise argparse.ArgumentTypeError(
            f"Need LABEL:PATH format, got: {s}"
        )
    label, path = s.split(":", 1)
    return (label.strip(), path.strip())


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
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
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
            regions.append((parts[0], int(parts[1]), int(parts[2]),
                             parts[3] if len(parts) >= 4
                             else f"r_{len(regions)+1}"))
    return regions


def extract_boundaries(tads):
    """A boundary = the END of one TAD / START of the next.
    Returns sorted list of unique boundary positions.
    """
    bounds = set()
    for _, s, e, _ in tads:
        bounds.add(s)
        bounds.add(e)
    return sorted(bounds)


def match_boundaries(cond_boundaries, tolerance):
    """Match boundary positions across conditions within tolerance.

    Input: dict {cond -> sorted list of boundary positions}
    Returns: list of "boundary events"; each event is a dict with:
      - representative_pos: median across matched conditions
      - by_cond: dict cond -> matched boundary position (or None)
      - status: "conserved" / "shifted" / "condition_specific"
      - max_shift_bp: max distance from representative
    """
    conds = sorted(cond_boundaries.keys())
    # Collect all boundaries with origin
    all_pts = []
    for c in conds:
        for p in cond_boundaries[c]:
            all_pts.append((p, c))
    all_pts.sort()

    # Group nearby boundaries
    events = []
    cur_group = []
    for p, c in all_pts:
        if not cur_group or p - cur_group[0][0] <= tolerance:
            cur_group.append((p, c))
        else:
            events.append(cur_group)
            cur_group = [(p, c)]
    if cur_group:
        events.append(cur_group)

    # Classify
    results = []
    n_conds = len(conds)
    for group in events:
        positions = [p for p, _ in group]
        cond_positions = {c: None for c in conds}
        for p, c in group:
            if cond_positions[c] is None:
                cond_positions[c] = p
            # If a cond has multiple in tolerance, use the first
        n_present = sum(1 for v in cond_positions.values() if v is not None)
        rep_pos = int(np.median(positions))
        present_positions = [p for p in cond_positions.values()
                              if p is not None]
        if len(present_positions) > 1:
            max_shift = max(present_positions) - min(present_positions)
        else:
            max_shift = 0

        if n_present == n_conds:
            status = "conserved" if max_shift <= tolerance else "shifted"
        else:
            status = "condition_specific"

        results.append({
            "representative_pos": rep_pos,
            "by_cond": cond_positions,
            "status": status,
            "max_shift_bp": max_shift,
            "n_present": n_present,
        })
    return results


def plot_tad_compare(tads_by_cond, insulation_by_cond, annotations,
                     boundary_events, region_start, region_end, chrom,
                     tolerance):
    """Stacked figure: features + (TAD bars + insulation) per condition +
    boundary classification track.
    """
    conds = sorted(tads_by_cond.keys())
    n_conds = len(conds)

    # Each condition gets a TAD-bar row and (optionally) an insulation row
    has_ins_by_cond = {c: c in insulation_by_cond and chrom in
                       insulation_by_cond[c] for c in conds}
    rows_per_cond = [(2 if has_ins_by_cond[c] else 1) for c in conds]

    n_panels = 1 + sum(rows_per_cond) + 1  # features + cond rows + bounds
    height_ratios = ([0.5] +
                      [v for c in conds for v in
                       ([0.5, 0.9] if has_ins_by_cond[c] else [0.5])] +
                      [0.6])

    fig = plt.figure(figsize=(14, max(7, sum(height_ratios) * 0.8)))
    gs = fig.add_gridspec(n_panels, 1, height_ratios=height_ratios,
                           hspace=0.15)

    # Panel 1: annotations
    ax_anno = fig.add_subplot(gs[0])
    ax_anno.set_ylim(0, 1.5)
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
                              xytext=((s2 + e2) / 2, 1.1),
                              ha="center", va="bottom", fontsize=7,
                              fontweight="bold",
                              arrowprops=dict(arrowstyle="-", color="navy",
                                               linewidth=0.5))
        else:
            ax_anno.text((s2 + e2) / 2, 0.5, n2,
                          ha="center", va="center",
                          fontsize=8, fontweight="bold")
    ax_anno.set_yticks([])
    ax_anno.set_ylabel("Features", fontsize=8, rotation=0,
                        ha="right", va="center")
    for sp in ["top", "right", "left"]:
        ax_anno.spines[sp].set_visible(False)
    ax_anno.tick_params(axis="x", labelbottom=False)

    # Per-condition panels
    panel_idx = 1
    cond_palette = ["#7799cc", "#88cc88", "#cc7777", "#cc8855", "#aa77cc"]
    for ci, cond in enumerate(conds):
        ax_tad = fig.add_subplot(gs[panel_idx], sharex=ax_anno)
        panel_idx += 1
        in_window = [t for t in tads_by_cond[cond]
                      if t[0] == chrom and t[2] > region_start
                      and t[1] < region_end]
        color_base = cond_palette[ci % len(cond_palette)]
        ax_tad.set_ylim(0, 1)
        for i, (_, t_start, t_end, _) in enumerate(in_window):
            face = color_base if i % 2 == 0 else "#cccccc"
            ax_tad.add_patch(plt.Rectangle(
                (t_start, 0.3), t_end - t_start, 0.4,
                facecolor=face, edgecolor="black", linewidth=0.5,
                alpha=0.75
            ))
        ax_tad.set_yticks([])
        ax_tad.set_ylabel(f"{cond}\nTADs (n={len(in_window)})",
                           fontsize=8, rotation=0, ha="right", va="center")
        ax_tad.tick_params(axis="x", labelbottom=False)

        if has_ins_by_cond[cond]:
            ax_ins = fig.add_subplot(gs[panel_idx], sharex=ax_anno)
            panel_idx += 1
            ins_data = insulation_by_cond[cond][chrom]
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
                                      color=color_base, alpha=0.3)
                ax_ins.axhline(0, color="gray", linestyle="--",
                                linewidth=0.5, alpha=0.5)
            ax_ins.set_ylabel(f"{cond}\nins.", fontsize=7, rotation=0,
                                ha="right", va="center")
            ax_ins.tick_params(axis="x", labelbottom=False)

    # Last panel: boundary classifications
    ax_bnd = fig.add_subplot(gs[panel_idx], sharex=ax_anno)
    ax_bnd.set_ylim(0, 1)
    color_map = {
        "conserved": "#888888",
        "shifted": "#dd9933",
        "condition_specific": "#cc3333",
    }
    in_window_events = [
        ev for ev in boundary_events
        if region_start <= ev["representative_pos"] <= region_end
    ]
    for ev in in_window_events:
        pos = ev["representative_pos"]
        col = color_map[ev["status"]]
        marker_size = 80 if ev["status"] == "conserved" else 140
        ax_bnd.scatter([pos], [0.5], c=col, s=marker_size,
                        zorder=3, edgecolor="black", linewidth=0.5)
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w",
                markerfacecolor=color_map["conserved"],
                markeredgecolor="black", markersize=8, label="Conserved"),
        Line2D([0], [0], marker="o", color="w",
                markerfacecolor=color_map["shifted"],
                markeredgecolor="black", markersize=10,
                label=f"Shifted (>{tolerance:,} bp)"),
        Line2D([0], [0], marker="o", color="w",
                markerfacecolor=color_map["condition_specific"],
                markeredgecolor="black", markersize=10,
                label="Condition-specific"),
    ]
    ax_bnd.legend(handles=legend_elements, loc="upper right", fontsize=7,
                    ncol=3)
    ax_bnd.set_yticks([])
    ax_bnd.set_ylabel("Boundary\nclassif.", fontsize=8, rotation=0,
                       ha="right", va="center")
    ax_bnd.set_xlim(region_start, region_end)
    ax_bnd.set_xlabel(f"{chrom} position (bp)", fontsize=10)
    ax_bnd.ticklabel_format(useOffset=False, style="plain", axis="x")

    fig.suptitle(
        f"Cross-condition TAD comparison: {chrom}:{region_start:,}-{region_end:,}\n"
        f"Conditions: {', '.join(conds)} | tolerance: {tolerance:,} bp",
        fontsize=11, y=0.995
    )
    return fig


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tads", nargs="+", required=True,
                    type=parse_labeled_path,
                    help="LABEL:PATH for each condition TAD BED")
    ap.add_argument("--insulation", nargs="*", default=[],
                    type=parse_labeled_path,
                    help="LABEL:PATH for each condition insulation .bm")
    ap.add_argument("--region", required=True, type=parse_region_str,
                    help="Region in format chrom:start-end")
    ap.add_argument("--annotations", default=None,
                    help="Feature annotation BED")
    ap.add_argument("--tolerance", type=int, default=50000,
                    help="Boundary-match tolerance in bp (default 50,000)")
    ap.add_argument("--out", required=True, help="Output prefix")
    args = ap.parse_args()

    chrom, region_start, region_end = args.region
    print(f"Region: {chrom}:{region_start:,}-{region_end:,}", flush=True)
    print(f"Tolerance: {args.tolerance:,} bp", flush=True)

    tads_by_cond = {}
    for cond, path in args.tads:
        tads = load_tads(path, chrom_filter=chrom,
                          start_filter=region_start, end_filter=region_end)
        tads_by_cond[cond] = tads
        print(f"  {cond}: {len(tads)} TADs in region from {Path(path).name}",
              flush=True)

    insulation_by_cond = {}
    for cond, path in args.insulation:
        ins = load_insulation(path, chrom_filter=chrom)
        insulation_by_cond[cond] = ins
        n = len(ins.get(chrom, []))
        print(f"  {cond}: {n} insulation bins on {chrom} from "
              f"{Path(path).name}", flush=True)

    annotations = load_annotations(args.annotations)
    print(f"Loaded {len(annotations)} feature annotations", flush=True)

    # Extract boundaries per condition (restricted to in-window)
    cond_boundaries = {}
    for cond, tads in tads_by_cond.items():
        bounds = sorted(set(extract_boundaries(tads)))
        bounds = [b for b in bounds if region_start <= b <= region_end]
        cond_boundaries[cond] = bounds

    # Match across conditions
    events = match_boundaries(cond_boundaries, args.tolerance)
    print(f"\nFound {len(events)} boundary events:", flush=True)
    by_status = {}
    for ev in events:
        by_status.setdefault(ev["status"], 0)
        by_status[ev["status"]] += 1
    for s, n in sorted(by_status.items()):
        print(f"  {s}: {n}", flush=True)

    # Write boundary calls TSV
    conds = sorted(tads_by_cond.keys())
    rows = []
    for ev in events:
        row = {
            "representative_pos": ev["representative_pos"],
            "status": ev["status"],
            "n_present": ev["n_present"],
            "max_shift_bp": ev["max_shift_bp"],
        }
        for c in conds:
            row[f"pos_{c}"] = ev["by_cond"][c]
        rows.append(row)
    df = pd.DataFrame(rows)
    bnd_tsv = f"{args.out}.boundaries.tsv"
    df.to_csv(bnd_tsv, sep="\t", index=False)
    print(f"\nWrote {bnd_tsv}", flush=True)

    # Pairwise shifts table: for each pair of conditions, what boundaries
    # are present in both but at different positions?
    shift_rows = []
    for i, c1 in enumerate(conds):
        for c2 in conds[i + 1:]:
            for ev in events:
                p1 = ev["by_cond"][c1]
                p2 = ev["by_cond"][c2]
                if p1 is not None and p2 is not None and abs(p1 - p2) > 0:
                    shift_rows.append({
                        "cond_a": c1,
                        "cond_b": c2,
                        "pos_a": p1,
                        "pos_b": p2,
                        "shift_bp": p2 - p1,
                        "abs_shift_bp": abs(p2 - p1),
                    })
    if shift_rows:
        sh_df = pd.DataFrame(shift_rows)
        sh_df = sh_df.sort_values("abs_shift_bp", ascending=False)
        sh_tsv = f"{args.out}.shifts.tsv"
        sh_df.to_csv(sh_tsv, sep="\t", index=False)
        print(f"Wrote {sh_tsv} ({len(sh_df)} boundary-shift records)",
              flush=True)

    # Plot
    pdf_path = f"{args.out}.pdf"
    with PdfPages(pdf_path) as pdf:
        fig = plot_tad_compare(
            tads_by_cond, insulation_by_cond, annotations,
            events, region_start, region_end, chrom, args.tolerance
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
