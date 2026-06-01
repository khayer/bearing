#!/usr/bin/env python3
"""
visualize_bearing_tracks.py
===========================
Track-stack view of the synthetic BEARING data, in the spirit of a genome
browser / pyGenomeTracks panel, but self-contained (matplotlib only).

For a chosen genomic window it stacks:
  * the six input signal tracks (ATAC, RNA+, RNA-, CTCF, RAD21, H3K27ac),
    each read straight from its BigWig, in the established track colors;
  * shaded spans marking the planted ground-truth blocks (from the blocks BED),
    labeled by block_id;
  * (optional) the per-bin total BEARING score and the per-track BEARING score
    stack, read from a qcat.bgz produced by bigwig_to_qcat.py, so you can see
    the planted signal AND where BEARING puts score, side by side.

This is the most direct sanity check that the simulator plants signal where the
truth table says, and that BEARING lights up there and stays quiet in null bins.

USAGE
    # signal tracks + planted blocks only
    python3 visualize_bearing_tracks.py \
        --prefix bearing_sim \
        --region chrSim:90000-260000 \
        --out tracks_named_blocks.png

    # also overlay the BEARING score (run the scorer first, or pass --qcat)
    python3 visualize_bearing_tracks.py \
        --prefix bearing_sim \
        --qcat bearing_sim.qcat.bgz \
        --region chrSim:90000-260000 \
        --out tracks_with_score.png

    # default region = the two named blocks (auto-located from the blocks BED)

DEPENDENCIES
    numpy, matplotlib, pyBigWig; pysam only if --qcat is given.
"""

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Fixed six-track order and the established BEARING track colors (matching the
# first six categories of bigwig_to_qcat.py: RAD21 == "Cohesin").
TRACK_NAMES = ["ATAC", "RNA+", "RNA-", "CTCF", "RAD21", "H3K27ac"]
TRACK_TOKENS = ["atac", "rnaplus", "rnaminus", "ctcf", "rad21", "h3k27ac"]
TRACK_COLORS = ["#be92e0", "#6495ed", "#1a3a8f", "#ff2200", "#8b0000", "#00e676"]
N_TRACKS = 6


def parse_region(region, chrom_default=None):
    """Parse 'chrom:start-end' (or 'chrom' for whole chrom). Returns (chrom, s, e)
    with s,e possibly None when only a chrom is given."""
    if region is None:
        return chrom_default, None, None
    if ":" not in region:
        return region, None, None
    chrom, span = region.split(":", 1)
    span = span.replace(",", "")
    s, e = span.split("-")
    return chrom, int(s), int(e)


def load_blocks_bed(path):
    """Return list of (chrom, start, end, name) from the planted-blocks BED."""
    out = []
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            out.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else ""))
    return out


def default_region_from_blocks(blocks, pad=20000):
    """Window spanning the two named blocks, with padding. Falls back to the
    first block if the named ones are absent."""
    named = [b for b in blocks if "named_" in b[3]]
    use = named if named else blocks[:1]
    if not use:
        return None, None, None
    chrom = use[0][0]
    s = min(b[1] for b in use) - pad
    e = max(b[2] for b in use) + pad
    return chrom, max(0, s), e


def read_bigwig_region(path, chrom, start, end, nbins):
    import pyBigWig
    bw = pyBigWig.open(path)
    # clamp end to chrom length
    clen = bw.chroms().get(chrom)
    if clen is None:
        bw.close()
        raise SystemExit("ERROR: chrom %s not in BigWig %s" % (chrom, path))
    end = min(end, clen)
    vals = bw.stats(chrom, start, end, type="mean", nBins=nbins)
    bw.close()
    return np.array([v if v is not None else 0.0 for v in vals], dtype=float)


def parse_qcat_region(qcat_bgz, chrom, start, end):
    """Return (centers, total_score, per_track_score [n x 6]) for bins in window."""
    import pysam
    import json
    tb = pysam.TabixFile(qcat_bgz)
    centers, totals, pertrack = [], [], []
    for line in tb.fetch(chrom, start, end):
        parts = line.rstrip("\n").split("\t")
        bs, be = int(parts[1]), int(parts[2])
        col = parts[3]
        i = col.find("qcat:")
        if i < 0:
            continue
        j = i + len("qcat:")
        depth, k_end = 0, None
        for k in range(j, len(col)):
            if col[k] == "[":
                depth += 1
            elif col[k] == "]":
                depth -= 1
                if depth == 0:
                    k_end = k + 1
                    break
        if k_end is None:
            continue
        pairs = json.loads(col[j:k_end])
        vec = np.zeros(N_TRACKS)
        for score, state in pairs:
            vec[int(state) - 1] = float(score)
        centers.append((bs + be) / 2.0)
        totals.append(vec.sum())
        pertrack.append(vec)
    if not centers:
        return (np.array([]), np.array([]), np.zeros((0, N_TRACKS)))
    return (np.array(centers), np.array(totals), np.vstack(pertrack))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Track-stack view of synthetic BEARING tracks with planted "
                    "blocks (and optional BEARING score overlay).")
    ap.add_argument("--prefix", required=True,
                    help="simulator prefix (expects <prefix>_<token>.bw and "
                         "<prefix>_blocks.bed)")
    ap.add_argument("--region", default=None,
                    help="chrom:start-end (default: window around the two "
                         "named blocks)")
    ap.add_argument("--qcat", default=None,
                    help="optional qcat.bgz to overlay the BEARING score")
    ap.add_argument("--nbins", type=int, default=1000,
                    help="number of display bins across the window (default 1000)")
    ap.add_argument("--out", default=None,
                    help="output PNG (default: <prefix>_tracks.png)")
    args = ap.parse_args(argv)

    blocks_bed = "%s_blocks.bed" % args.prefix
    blocks = load_blocks_bed(blocks_bed)

    bw_paths = ["%s_%s.bw" % (args.prefix, t) for t in TRACK_TOKENS]
    for p in bw_paths:
        if not os.path.isfile(p):
            sys.exit("ERROR: missing BigWig %s" % p)

    chrom, start, end = parse_region(args.region)
    if start is None:
        chrom, start, end = default_region_from_blocks(blocks)
        if chrom is None:
            sys.exit("ERROR: no --region given and no blocks BED to derive one")
    have_score = args.qcat is not None and os.path.isfile(args.qcat)

    # Decide the set of windows to draw as columns.
    #   - explicit --region  -> single window (chrom, start, end)
    #   - otherwise          -> column 1 = whole-chromosome overview,
    #                           columns 2.. = zoom-ins on the named blocks
    if args.region is not None:
        windows = [(chrom, start, end, "%s:%s-%s"
                    % (chrom, f"{start:,}", f"{end:,}"))]
    else:
        # whole chromosome overview
        chrom_len = max(b[2] for b in blocks) if blocks else end
        # try to read the true chrom length from a chrom.sizes if present
        cs = "%s.chrom.sizes" % args.prefix
        if os.path.isfile(cs):
            with open(cs) as fh:
                for line in fh:
                    c, L = line.split()[:2]
                    if c == chrom:
                        chrom_len = int(L)
                        break
        windows = [(chrom, 0, chrom_len, "whole chromosome (overview)")]
        # zoom-ins: the named blocks, padded
        pad = 15000
        named = [b for b in blocks if "named_" in b[3]]
        for (_, bs, be, name) in named:
            windows.append((chrom, max(0, bs - pad), be + pad,
                            "zoom: %s" % name.replace("named_", "")))

    n_signal_rows = N_TRACKS
    n_extra = 2 if have_score else 0
    n_rows = n_signal_rows + n_extra
    n_cols = len(windows)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.0 * n_cols, 1.05 * n_rows),
                             squeeze=False)

    bw_cache = {}

    def render_window(col, chrom, start, end, col_title):
        nbins = args.nbins
        xs = np.linspace(start, end, nbins)
        win_blocks = [b for b in blocks
                      if b[0] == chrom and b[2] > start and b[1] < end]

        def shade(ax):
            for (_, bs, be, _name) in win_blocks:
                ax.axvspan(max(bs, start), min(be, end), color="0.85", zorder=0)

        def label_top(ax):
            # only label in overview if not too crowded, and always in zooms
            if len(win_blocks) > 12:
                return
            for (_, bs, be, name) in win_blocks:
                xc = (max(bs, start) + min(be, end)) / 2.0
                ax.annotate(name, xy=(xc, 1.04),
                            xycoords=("data", "axes fraction"),
                            rotation=90, va="bottom", ha="center", fontsize=4.5,
                            color="0.35", annotation_clip=False)

        # six signal tracks
        for ti in range(N_TRACKS):
            ax = axes[ti][col]
            key = (bw_paths[ti], chrom, start, end, nbins)
            if key not in bw_cache:
                bw_cache[key] = read_bigwig_region(bw_paths[ti], chrom,
                                                   start, end, nbins)
            vals = bw_cache[key]
            shade(ax)
            if ti == 0:
                label_top(ax)
            ax.fill_between(xs, 0, vals, color=TRACK_COLORS[ti], lw=0)
            if col == 0:
                ax.set_ylabel(TRACK_NAMES[ti], rotation=0, ha="right",
                              va="center", fontsize=8)
            ax.set_yticks([])
            for sp in ("top", "right", "left"):
                ax.spines[sp].set_visible(False)

        # optional score rows
        if have_score:
            centers, totals, pertrack = parse_qcat_region(args.qcat, chrom,
                                                          start, end)
            ax_tot = axes[N_TRACKS][col]
            shade(ax_tot)
            if len(centers):
                ax_tot.fill_between(centers, 0, totals, color="0.2", lw=0)
            if col == 0:
                ax_tot.set_ylabel("BEARING\nscore", rotation=0, ha="right",
                                  va="center", fontsize=8)
            ax_tot.set_yticks([])
            for sp in ("top", "right", "left"):
                ax_tot.spines[sp].set_visible(False)

            ax_stk = axes[N_TRACKS + 1][col]
            shade(ax_stk)
            if len(centers):
                bottom = np.zeros(len(centers))
                width = (centers[1] - centers[0]) if len(centers) > 1 else \
                    (end - start) / max(nbins, 1)
                for ti in range(N_TRACKS):
                    ax_stk.bar(centers, pertrack[:, ti], bottom=bottom,
                               width=width, color=TRACK_COLORS[ti], lw=0,
                               label=TRACK_NAMES[ti] if col == 0 else None)
                    bottom += pertrack[:, ti]
            if col == 0:
                ax_stk.set_ylabel("per-track\nscore", rotation=0, ha="right",
                                  va="center", fontsize=8)
                ax_stk.legend(ncol=6, fontsize=6, loc="upper center",
                              bbox_to_anchor=(0.5, -0.45))
            ax_stk.set_yticks([])
            for sp in ("top", "right", "left"):
                ax_stk.spines[sp].set_visible(False)

        axes[0][col].set_title(col_title, fontsize=9)
        axes[-1][col].set_xlabel("%s position (bp)" % chrom, fontsize=8)
        axes[-1][col].set_xlim(start, end)
        for r in range(n_rows):
            axes[r][col].tick_params(axis="x", labelsize=6)
        return len(win_blocks)

    total_blocks = 0
    for col, (c, s, e, t) in enumerate(windows):
        total_blocks += render_window(col, c, s, e, t)

    title = "Synthetic BEARING tracks  (%s)" % chrom
    if have_score:
        title += "   (+ BEARING score)"
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])

    out = args.out or ("%s_tracks.png" % args.prefix)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote %s  (%d window(s): %s%s)"
          % (out, len(windows), ", ".join(w[3] for w in windows),
             ", + BEARING score" if have_score else ""))


if __name__ == "__main__":
    main()
