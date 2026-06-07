#!/usr/bin/env python3
"""
diffuse_architecture_top1pct.py

Quantify how "diffuse" or "focal" each track's top differential signal is,
across a locus of interest.

Strategy:
  1. Load a BEARING per-bin differential table (one comparison) with columns
     chrom, start, end, kl_1..kl_6 (per-track KL contributions).
  2. For each track, compute the 99th-percentile threshold of |kl_T| within
     the analysis window (or genome-wide -- selectable).
  3. Tile the locus into non-overlapping 100 kb windows.
  4. Count top-1% bins per track per window.
  5. Compute focal-vs-flank enrichment for a user-specified focal sub-region
     (Fisher's exact, two-sided).

Output:
  - Per-track per-window count TSV (long format).
  - Per-track focal-vs-flank summary TSV.
  - Optional matplotlib figure showing top-1% count per 100 kb per track.

Track order is fixed: ATAC, RNA+, RNA-, CTCF, RAD21 (Cohesin), H3K27ac.

Usage (Tcrb DN-vs-DP, recombination center as focal region):
  python3 diffuse_architecture_top1pct.py \
      --diff diff_DN_vs_DP.tsv \
      --region chr6:40400000-42400000 \
      --focal chr6:41000000-41700000 \
      --window 100000 \
      --pct 99 \
      --threshold-scope window \
      --abs \
      --out-prefix tcrb_DN_vs_DP_diffuse

ASCII-only. No emoji, no unicode.
"""

import argparse
import sys

import numpy as np
import pandas as pd


# Track order is fixed (per project convention).
TRACK_NAMES = ["ATAC", "RNA_plus", "RNA_minus", "CTCF", "RAD21", "H3K27ac"]
KL_COLS = ["kl_1", "kl_2", "kl_3", "kl_4", "kl_5", "kl_6"]


def parse_locus(spec):
    """Parse 'chr6:40400000-42400000' -> ('chr6', 40400000, 42400000)."""
    chrom, rest = spec.split(":")
    start_s, end_s = rest.split("-")
    return chrom, int(start_s.replace(",", "")), int(end_s.replace(",", ""))


def load_diff(path):
    """Load the BEARING differential TSV. Expects header with chrom,start,end,kl_1..kl_6.

    Transparently handles gzip/bgzip-compressed input (detected by the gzip
    magic bytes) so a .tsv or a .tsv.gz both work.
    """
    compression = "infer"
    try:
        with open(path, "rb") as _fh:
            if _fh.read(2) == b"\x1f\x8b":
                compression = "gzip"
    except OSError:
        pass
    df = pd.read_csv(path, sep="\t", compression=compression)
    needed = {"chrom", "start", "end"} | set(KL_COLS)
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit("ERROR: input file missing required columns: %s" % sorted(missing))
    return df


def subset_window(df, chrom, start, end):
    """Restrict df to the genomic window."""
    m = (df["chrom"] == chrom) & (df["start"] < end) & (df["end"] > start)
    return df.loc[m].copy().reset_index(drop=True)


def compute_thresholds(df, pct, use_abs):
    """For each track, return the percentile threshold of |kl_T| (or kl_T)."""
    out = {}
    for col, name in zip(KL_COLS, TRACK_NAMES):
        vals = df[col].to_numpy()
        if use_abs:
            vals = np.abs(vals)
        # Drop NaN if any
        vals = vals[~np.isnan(vals)]
        if vals.size == 0:
            out[name] = np.nan
        else:
            out[name] = float(np.percentile(vals, pct))
    return out


def per_track_top_mask(df, thresholds, use_abs):
    """Boolean mask per track: True where |kl_T| (or kl_T) >= threshold."""
    masks = {}
    for col, name in zip(KL_COLS, TRACK_NAMES):
        v = df[col].to_numpy()
        if use_abs:
            v = np.abs(v)
        masks[name] = v >= thresholds[name]
    return masks


def tile_windows(start, end, window):
    """Yield (w_start, w_end) windows of given size tiling [start, end)."""
    cur = start
    while cur < end:
        nxt = min(cur + window, end)
        yield (cur, nxt)
        cur = nxt


def count_per_window(df, masks, chrom, win_start, win_end):
    """Count, per track, the number of top-1% bins whose midpoint falls in the window."""
    mid = ((df["start"] + df["end"]) // 2).to_numpy()
    in_win = (mid >= win_start) & (mid < win_end)
    total_bins = int(in_win.sum())
    counts = {name: int((masks[name] & in_win).sum()) for name in TRACK_NAMES}
    return counts, total_bins


def load_genes(path, chrom, start, end):
    """
    Load BED gene annotation, returning a list of dicts for records overlapping
    the analysis region. Accepts BED6 or BED9; uses cols 1-4 and 6 (strand).
    """
    if path is None:
        return []
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            c = parts[0]
            try:
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            if c != chrom:
                continue
            if e <= start or s >= end:
                continue
            name = parts[3] if len(parts) > 3 else ""
            strand = parts[5] if len(parts) > 5 else "."
            rows.append({"chrom": c, "start": s, "end": e, "name": name, "strand": strand})
    rows.sort(key=lambda r: r["start"])
    return rows


def render_genes_panel(ax, genes, region_start_mb, region_end_mb,
                       focal_start_mb, focal_end_mb,
                       label_fontsize, max_labels):
    """
    Draw an arrow-box gene track. Boxes are vertically staggered if labels
    would overlap. Strand encoded as left/right arrow shape via a polygon.
    region_*_mb and focal_*_mb are in Mb (matches the other panels' x-axis).
    Gene boxes have a minimum visible width to handle very short genes
    (e.g. J segments at antigen receptor loci that are ~50 bp).
    """
    ax.set_xlim(region_start_mb, region_end_mb)
    ax.set_yticks([])
    ax.set_ylabel("genes", rotation=0, ha="right", va="center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Focal shading consistent with other panels
    if focal_end_mb > focal_start_mb:
        ax.axvspan(focal_start_mb, focal_end_mb, alpha=0.12, color="#cccccc", zorder=0)

    if not genes:
        ax.set_ylim(0, 1.0)
        ax.text((region_start_mb + region_end_mb) / 2.0, 0.5,
                "no gene annotation in region", ha="center", va="center",
                fontsize=label_fontsize, color="#888888")
        return

    skip_labels = len(genes) > max_labels

    region_span_mb = region_end_mb - region_start_mb

    # Minimum visible gene width: 0.4% of region width. For very short genes
    # (J segments etc.) this prevents the box from being a single-pixel hairline.
    min_visible_width = region_span_mb * 0.004

    # Label width heuristic: ~1.2% of region width per character. Used only
    # for staggering, not for the boxes themselves.
    char_width_mb = region_span_mb * 0.012

    # Layout: 6 tiers, alternating to spread out crowded regions.
    # Tier 0 is the bottom (lowest), tier 5 is the top.
    tier_count = 6
    box_height = 0.10
    tier_gap = 0.13
    # y positions for the bottom of each tier's box
    tier_y = [0.05 + i * tier_gap for i in range(tier_count)]
    # The y-limit needs to accommodate the top tier's label
    y_max = tier_y[-1] + box_height + 0.18
    ax.set_ylim(0, y_max)

    # Track the rightmost extent (box + reserved label) per tier
    last_end_per_tier = [-1e18] * tier_count

    for g in genes:
        gs = g["start"] / 1e6
        ge = g["end"] / 1e6
        name = g["name"]
        strand = g["strand"]

        # Enforce minimum visible width by extending around the gene midpoint
        if (ge - gs) < min_visible_width:
            mid = (gs + ge) / 2.0
            box_gs = mid - min_visible_width / 2.0
            box_ge = mid + min_visible_width / 2.0
        else:
            box_gs = gs
            box_ge = ge

        # Pick the lowest tier whose last reserved x-extent is left of this gene's start
        needed_label_width = max(len(name), 1) * char_width_mb if not skip_labels else 0
        # Reserve the wider of (box, box + label) for collision purposes
        reserve_end = max(box_ge, box_gs + needed_label_width)
        tier = 0
        # find first available tier
        for t in range(tier_count):
            if box_gs >= last_end_per_tier[t]:
                tier = t
                break
        else:
            # all tiers occupied; pick the one with the smallest right extent
            tier = min(range(tier_count), key=lambda t: last_end_per_tier[t])
        last_end_per_tier[tier] = reserve_end

        y0 = tier_y[tier]
        # Arrow box: simple rectangle with a triangle tail/head depending on strand.
        gene_len = box_ge - box_gs
        head_width = min(0.006 * region_span_mb, gene_len * 0.4)
        if strand == "+":
            poly_x = [box_gs, box_ge - head_width, box_ge - head_width, box_ge,
                      box_ge - head_width, box_ge - head_width, box_gs]
            poly_y = [y0, y0, y0 - 0.012, y0 + box_height / 2.0,
                      y0 + box_height + 0.012, y0 + box_height, y0 + box_height]
        elif strand == "-":
            poly_x = [box_gs + head_width, box_ge, box_ge, box_gs + head_width,
                      box_gs + head_width, box_gs, box_gs + head_width]
            poly_y = [y0, y0, y0 + box_height, y0 + box_height,
                      y0 + box_height + 0.012, y0 + box_height / 2.0, y0 - 0.012]
        else:
            poly_x = [box_gs, box_ge, box_ge, box_gs]
            poly_y = [y0, y0, y0 + box_height, y0 + box_height]

        ax.fill(poly_x, poly_y, color="#444444", edgecolor="none", zorder=2)

        if not skip_labels and name:
            ax.text((box_gs + box_ge) / 2.0, y0 + box_height + 0.025, name,
                    ha="center", va="bottom", fontsize=label_fontsize,
                    color="#222222", zorder=3, clip_on=False)


def fisher_focal_flank(focal_top, focal_total, flank_top, flank_total):
    """
    Two-sided Fisher's exact for a 2x2 contingency:
        [[focal_top,            focal_total - focal_top],
         [flank_top,            flank_total - flank_top]]
    Returns (odds_ratio, p_value). Uses scipy if present, else gives an OR
    estimate and skips the p-value (printed as NaN).
    """
    try:
        from scipy.stats import fisher_exact
        table = [[focal_top, max(0, focal_total - focal_top)],
                 [flank_top, max(0, flank_total - flank_top)]]
        odds, p = fisher_exact(table, alternative="two-sided")
        return float(odds), float(p)
    except Exception:
        # Fall back to a Haldane-Anscombe corrected OR estimate, no p-value
        a = focal_top + 0.5
        b = max(0, focal_total - focal_top) + 0.5
        c = flank_top + 0.5
        d = max(0, flank_total - flank_top) + 0.5
        return float((a * d) / (b * c)), float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diff", required=True,
                    help="BEARING per-bin differential TSV with chrom,start,end,kl_1..kl_6")
    ap.add_argument("--region", required=True,
                    help="Analysis window, e.g. chr6:40400000-42400000")
    ap.add_argument("--focal", required=True,
                    help="Focal sub-region (e.g. recombination center) chr6:41000000-41700000")
    ap.add_argument("--window", type=int, default=100000,
                    help="Tiling window size in bp (default 100000)")
    ap.add_argument("--pct", type=float, default=99.0,
                    help="Percentile threshold (default 99 = top 1%%)")
    ap.add_argument("--threshold-scope", choices=["window", "genome"], default="window",
                    help="Compute percentile threshold within the analysis window (default) "
                         "or genome-wide (whole diff table). Window-scope is more conservative "
                         "for locus-focused claims.")
    ap.add_argument("--abs", action="store_true", default=True,
                    help="Use |kl_T| rather than signed kl_T (default True)")
    ap.add_argument("--signed", dest="abs", action="store_false",
                    help="Use signed kl_T instead of absolute value")
    ap.add_argument("--out-prefix", required=True,
                    help="Prefix for output files (.per_window.tsv and .focal_vs_flank.tsv)")
    ap.add_argument("--figure", default=None,
                    help="Optional path for matplotlib figure (png/pdf). If omitted, no figure.")
    ap.add_argument("--genes", default=None,
                    help="Optional BED file (BED6 or BED9) for gene annotation. Columns used: "
                         "chrom, start, end, name, _, strand. Records overlapping the analysis "
                         "region are drawn as a track at the bottom of the figure.")
    ap.add_argument("--gene-label-fontsize", type=float, default=6.5,
                    help="Font size for gene labels (default 6.5). Reduce for dense regions.")
    ap.add_argument("--max-gene-labels", type=int, default=60,
                    help="If region has more than this many overlapping genes, draw boxes only "
                         "and skip labels to avoid overprinting (default 60).")
    args = ap.parse_args()

    # Load and subset
    df_all = load_diff(args.diff)
    region_chrom, region_start, region_end = parse_locus(args.region)
    focal_chrom, focal_start, focal_end = parse_locus(args.focal)
    if focal_chrom != region_chrom:
        raise SystemExit("ERROR: focal chrom (%s) != region chrom (%s)" % (focal_chrom, region_chrom))
    if not (region_start <= focal_start and focal_end <= region_end):
        sys.stderr.write("WARN: focal region is not entirely inside the analysis window\n")

    df_region = subset_window(df_all, region_chrom, region_start, region_end)
    if df_region.empty:
        raise SystemExit("ERROR: no bins in region %s" % args.region)
    sys.stderr.write("loaded %d bins in region %s\n" % (len(df_region), args.region))

    # Compute thresholds
    if args.threshold_scope == "window":
        thr_df = df_region
    else:
        thr_df = df_all
    thresholds = compute_thresholds(thr_df, args.pct, args.abs)
    sys.stderr.write("track thresholds (p%g, abs=%s, scope=%s):\n" %
                     (args.pct, args.abs, args.threshold_scope))
    for n in TRACK_NAMES:
        sys.stderr.write("  %-10s %.6f\n" % (n, thresholds[n]))

    # Compute top masks on the region
    masks = per_track_top_mask(df_region, thresholds, args.abs)

    # Tile windows and count
    rows = []
    for w_start, w_end in tile_windows(region_start, region_end, args.window):
        counts, total = count_per_window(df_region, masks, region_chrom, w_start, w_end)
        in_focal = (w_start < focal_end) and (w_end > focal_start)
        row = {
            "chrom": region_chrom,
            "win_start": w_start,
            "win_end": w_end,
            "in_focal": int(in_focal),
            "total_bins": total,
        }
        for n in TRACK_NAMES:
            row["top_%s" % n] = counts[n]
        rows.append(row)
    per_win = pd.DataFrame(rows)
    per_win_path = args.out_prefix + ".per_window.tsv"
    per_win.to_csv(per_win_path, sep="\t", index=False)
    sys.stderr.write("wrote %s (%d windows)\n" % (per_win_path, len(per_win)))

    # Focal-vs-flank summary
    focal_mask = per_win["in_focal"] == 1
    flank_mask = ~focal_mask
    rows2 = []
    for n in TRACK_NAMES:
        col = "top_%s" % n
        focal_top = int(per_win.loc[focal_mask, col].sum())
        focal_total = int(per_win.loc[focal_mask, "total_bins"].sum())
        flank_top = int(per_win.loc[flank_mask, col].sum())
        flank_total = int(per_win.loc[flank_mask, "total_bins"].sum())
        focal_rate = focal_top / focal_total if focal_total > 0 else float("nan")
        flank_rate = flank_top / flank_total if flank_total > 0 else float("nan")
        odds, p = fisher_focal_flank(focal_top, focal_total, flank_top, flank_total)
        rows2.append({
            "track": n,
            "focal_top": focal_top,
            "focal_total_bins": focal_total,
            "focal_top1pct_rate": focal_rate,
            "flank_top": flank_top,
            "flank_total_bins": flank_total,
            "flank_top1pct_rate": flank_rate,
            "focal_vs_flank_OR": odds,
            "fisher_p_two_sided": p,
        })
    summary = pd.DataFrame(rows2)
    summary_path = args.out_prefix + ".focal_vs_flank.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    sys.stderr.write("wrote %s\n" % summary_path)

    # Print a short readable summary
    sys.stderr.write("\nfocal-vs-flank summary (OR > 1 = focal-enriched, OR ~ 1 = diffuse):\n")
    sys.stderr.write("  %-10s  %12s  %12s  %8s  %10s\n" %
                     ("track", "focal_rate", "flank_rate", "OR", "p"))
    for _, r in summary.iterrows():
        sys.stderr.write("  %-10s  %12.4f  %12.4f  %8.3f  %10.3g\n" %
                         (r["track"], r["focal_top1pct_rate"],
                          r["flank_top1pct_rate"], r["focal_vs_flank_OR"],
                          r["fisher_p_two_sided"]))

    # Optional figure
    if args.figure is not None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.gridspec import GridSpec
            colors = {
                "ATAC":      "#b39ddb",
                "RNA_plus":  "#6495ed",
                "RNA_minus": "#1a3a8f",
                "CTCF":      "#ff2200",
                "RAD21":     "#8b0000",
                "H3K27ac":   "#00c864",
            }
            # Optionally load gene annotation
            genes = load_genes(args.genes, region_chrom, region_start, region_end) \
                if args.genes else []
            has_genes_panel = args.genes is not None
            sys.stderr.write("gene records overlapping region: %d\n" % len(genes))

            # Layout: 6 track panels (height 1) + optional gene panel (height 2.2)
            n_tracks = len(TRACK_NAMES)
            if has_genes_panel:
                height_ratios = [1.0] * n_tracks + [2.2]
                fig_height = 1.1 * n_tracks + 2.2
                n_rows = n_tracks + 1
            else:
                height_ratios = [1.0] * n_tracks
                fig_height = 1.1 * n_tracks
                n_rows = n_tracks
            fig = plt.figure(figsize=(10, fig_height))
            gs = GridSpec(n_rows, 1, height_ratios=height_ratios,
                          hspace=0.15, figure=fig)
            track_axes = [fig.add_subplot(gs[i, 0]) for i in range(n_tracks)]
            gene_ax = fig.add_subplot(gs[n_tracks, 0], sharex=track_axes[0]) \
                if has_genes_panel else None
            # Share x across track axes manually
            for ax in track_axes[1:]:
                ax.sharex(track_axes[0])

            mids_mb = (per_win["win_start"] + per_win["win_end"]) / 2.0 / 1e6
            for ax, n in zip(track_axes, TRACK_NAMES):
                ax.bar(mids_mb, per_win["top_%s" % n],
                       width=args.window / 1e6 * 0.9,
                       color=colors.get(n, "#666666"))
                ax.axvspan(focal_start / 1e6, focal_end / 1e6,
                           alpha=0.12, color="#cccccc", zorder=0)
                ax.set_ylabel(n, rotation=0, ha="right", va="center")
                ax.tick_params(left=True, labelleft=True)
                # Hide x tick labels on all but the last visible axis
                ax.tick_params(labelbottom=False)

            # Gene panel
            if gene_ax is not None:
                render_genes_panel(
                    gene_ax,
                    genes,
                    region_start / 1e6, region_end / 1e6,
                    focal_start / 1e6, focal_end / 1e6,
                    label_fontsize=args.gene_label_fontsize,
                    max_labels=args.max_gene_labels,
                )
                gene_ax.set_xlabel("%s position (Mb)" % region_chrom)
            else:
                track_axes[-1].set_xlabel("%s position (Mb)" % region_chrom)
                track_axes[-1].tick_params(labelbottom=True)

            # Title on top axis
            track_axes[0].set_title(
                "Top %g%% |kl| count per %d kb window (focal=grey)"
                % (100 - args.pct, args.window // 1000))

            fig.tight_layout()
            fig.savefig(args.figure, dpi=150)
            sys.stderr.write("wrote %s\n" % args.figure)
        except Exception as e:
            sys.stderr.write("figure not produced: %s\n" % e)

    sys.stderr.write("done.\n")


if __name__ == "__main__":
    main()
