#!/usr/bin/env python3
"""
consolidate_regional_enrichment.py
==================================
Consolidate per-comparison regional-enrichment TSVs (from
regional_enrichment.py batch) into a single long-format table, and optionally
render a region x comparison heatmap.

Each input TSV has one row per (comparison, region_name) with columns:
  comparison, region_name, chrom, start, end, L_region_bins, L_locus_bins,
  pi, k, n_locus, k_pos, k_neg, p_spatial, p_directional, p_combined, q_combined

Consolidation simply stacks all input rows (they already carry a `comparison`
column) into one table, preserving every column.

Two entry points are provided so the workflow can call one script for both
steps:
  * default / --out            : write the consolidated long table
  * --heatmap-out (optional)   : also render a heatmap from the consolidated
                                 table in the same run

A companion mode (--from-consolidated) renders a heatmap from an
already-consolidated table, so the workflow can run consolidation and heatmap
as two invocations of this one script if preferred.

ASCII-only source; matplotlib only imported when a heatmap is requested.
"""

import argparse
import csv
import sys


def read_tsv(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader), (reader.fieldnames or [])


def consolidate(tsv_paths, out_path):
    """Stack per-comparison TSVs into one long table. Returns (rows, header)."""
    all_rows = []
    header = None
    for p in tsv_paths:
        rows, fields = read_tsv(p)
        if not fields:
            continue
        if header is None:
            header = list(fields)
        all_rows.extend(rows)
    if header is None:
        sys.exit("ERROR: no readable input TSVs (no header found in any).")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for r in all_rows:
            writer.writerow({k: r.get(k, "") for k in header})
    sys.stderr.write(
        "Wrote consolidated table: %s (%d rows from %d files)\n"
        % (out_path, len(all_rows), len(tsv_paths))
    )
    return all_rows, header


def render_heatmap(rows, out_path, comparisons=None, value_col="q_combined"):
    """Render a region x comparison heatmap of -log10(value_col).

    Directionality (net sign of k_pos - k_neg) is overlaid as +/- annotations
    so a reader sees both significance (color) and direction (text).
    """
    import math

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Region order: first appearance; comparison order: as given or first seen.
    regions = []
    comps = []
    cell = {}        # (region, comparison) -> -log10(value)
    sign = {}        # (region, comparison) -> "+"/"-"/""
    for r in rows:
        reg = r.get("region_name", "")
        comp = r.get("comparison", "")
        if not reg or not comp:
            continue
        if reg not in regions:
            regions.append(reg)
        if comp not in comps:
            comps.append(comp)
        try:
            q = float(r.get(value_col, "nan"))
        except ValueError:
            q = float("nan")
        nlp = (-math.log10(q)) if (q == q and q > 0) else 0.0
        cell[(reg, comp)] = nlp
        try:
            kp = float(r.get("k_pos", "0") or 0)
            kn = float(r.get("k_neg", "0") or 0)
        except ValueError:
            kp = kn = 0.0
        sign[(reg, comp)] = "+" if kp > kn else ("-" if kn > kp else "")

    if comparisons:
        comps = [c for c in comparisons if c in comps] or comps

    if not regions or not comps:
        sys.exit("ERROR: no (region, comparison) pairs to plot.")

    mat = np.zeros((len(regions), len(comps)), dtype=float)
    for i, reg in enumerate(regions):
        for j, comp in enumerate(comps):
            mat[i, j] = cell.get((reg, comp), 0.0)

    fig_w = max(4.0, 1.2 * len(comps) + 2.0)
    fig_h = max(3.0, 0.6 * len(regions) + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0.0)
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels(comps, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(regions)))
    ax.set_yticklabels(regions, fontsize=8)
    for i in range(len(regions)):
        for j in range(len(comps)):
            s = sign.get((regions[i], comps[j]), "")
            if s:
                ax.text(j, i, s, ha="center", va="center",
                        color="white", fontsize=9, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("-log10(%s)" % value_col, fontsize=8)
    ax.set_title("Regional enrichment (-log10 %s); +/- = net direction"
                 % value_col, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    sys.stderr.write("Wrote heatmap: %s\n" % out_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsvs", nargs="+",
                    help="Per-comparison enrich_*.tsv files to consolidate.")
    ap.add_argument("--out",
                    help="Output consolidated long-format TSV.")
    ap.add_argument("--heatmap-out",
                    help="If set, also render a heatmap to this path in the "
                         "same run (after consolidation).")
    ap.add_argument("--from-consolidated",
                    help="Render a heatmap from an already-consolidated TSV "
                         "(skip consolidation).")
    ap.add_argument("--comparisons", nargs="*", default=None,
                    help="Optional explicit comparison column order for the "
                         "heatmap.")
    ap.add_argument("--value-col", default="q_combined",
                    help="Column to color the heatmap by (default q_combined).")
    args = ap.parse_args()

    if args.from_consolidated:
        rows, _ = read_tsv(args.from_consolidated)
        out = args.heatmap_out or args.out
        if not out:
            sys.exit("ERROR: --from-consolidated needs --heatmap-out or --out.")
        render_heatmap(rows, out, comparisons=args.comparisons,
                       value_col=args.value_col)
        return

    if not args.tsvs or not args.out:
        sys.exit("ERROR: provide --tsvs and --out (or use --from-consolidated).")

    rows, _ = consolidate(args.tsvs, args.out)
    if args.heatmap_out:
        render_heatmap(rows, args.heatmap_out, comparisons=args.comparisons,
                       value_col=args.value_col)


if __name__ == "__main__":
    main()
