#!/usr/bin/env python3
"""
plot_hic_superbins.py -- render BEARING super-bin Hi-C at a locus, aligned to the
BEARING 1D track and an antigen-receptor (AgR) annotation, so the 1D chromatin
change and the 3D contact change can be read on the same coordinates.

Two heatmap modes:
  single : one cooler -> width-normalized contact DENSITY  count / (w_i * w_j)
           (log color scale). Use to see "how the cool file looks" at a locus.
  ratio  : two coolers A,B -> log2( density_A / density_B ), masking any pixel
           that is zero in either condition. Because A and B share one grid, the
           w_i*w_j factors cancel exactly, so this equals log2(count_A/count_B);
           the densities are still formed explicitly per the requested recipe and
           so the code stays correct if grids ever differ.

Panels (top to bottom, shared genomic x-axis):
  1. AgR annotation track  (--annot AgR_mm10.bed; boxes coloured V/D/J/C/other)
  2. BEARING 1D track      (--bearing-track *.qcat.bgz | *.bedgraph | *.bw)
  3. Hi-C heatmap          (square contact map over the region)

Requires: cooler, numpy, pandas, matplotlib. qcat.bgz reading uses pysam if a
.tbi is present, else a gzip stream. bigwig needs pyBigWig (optional). ASCII only.
"""
import argparse
import gzip
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.patches import Rectangle


# ---------------------------------------------------------------- region utils
def parse_region(s):
    chrom, rng = s.split(":")
    a, b = (int(x.replace(",", "")) for x in rng.split("-"))
    return chrom, a, b


# ------------------------------------------------------------- cooler fetching
def fetch_region(cool_path, region):
    import cooler
    c = cooler.Cooler(cool_path)
    bins = c.bins().fetch(region).reset_index(drop=True)
    mat = c.matrix(balance=False, sparse=False).fetch(region).astype(np.float64)
    starts = bins["start"].to_numpy()
    ends = bins["end"].to_numpy()
    width = bins["width"].to_numpy().astype(np.float64) if "width" in bins \
        else (ends - starts).astype(np.float64)
    edges = np.append(starts, ends[-1]).astype(np.float64)
    return edges, mat, width, starts, ends


# ---------------------------------------------------------- BEARING 1D readers
def read_qcat(path, chrom, a, b):
    """Return (starts, ends, total_score) over [a,b). Sums KL scores per bin."""
    rows = []

    def handle(line):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or f[0] != chrom:
            return None
        s, e = int(f[1]), int(f[2])
        if e <= a or s >= b:
            return "skip"
        j = f[3].find("qcat:")
        tot = 0.0
        if j >= 0:
            try:
                pairs = json.loads(f[3][j + 5:])
                tot = float(sum(p[0] for p in pairs))
            except Exception:
                tot = 0.0
        rows.append((s, e, tot))
        return "ok"

    used_tabix = False
    if os.path.exists(path + ".tbi"):
        try:
            import pysam
            tb = pysam.TabixFile(path)
            for line in tb.fetch(chrom, a, b):
                handle(line + "\n")
            used_tabix = True
        except Exception:
            used_tabix = False
    if not used_tabix:
        opener = gzip.open
        entered = False
        with opener(path, "rt") as fh:
            for line in fh:
                r = handle(line)
                if r == "ok":
                    entered = True
                elif entered and line.split("\t", 1)[0] == chrom and r == "skip":
                    # past the region on this chrom -> stop (file is sorted)
                    s = int(line.split("\t")[1])
                    if s >= b:
                        break
    if not rows:
        return None
    rows.sort()
    return (np.array([r[0] for r in rows]),
            np.array([r[1] for r in rows]),
            np.array([r[2] for r in rows]))


def read_bedgraph(path, chrom, a, b):
    op = gzip.open if path.endswith(".gz") else open
    s_, e_, v_ = [], [], []
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line[0] in "#tb":
                continue
            f = line.split()
            if f[0] != chrom:
                continue
            s, e = int(f[1]), int(f[2])
            if e <= a or s >= b:
                continue
            s_.append(s); e_.append(e); v_.append(float(f[3]))
    if not s_:
        return None
    return np.array(s_), np.array(e_), np.array(v_)


def read_bigwig(path, chrom, a, b, nbins=600):
    import pyBigWig
    bw = pyBigWig.open(path)
    vals = np.array(bw.stats(chrom, a, b, nBins=nbins, type="mean"), dtype=np.float64)
    bw.close()
    vals = np.nan_to_num(vals)
    edges = np.linspace(a, b, nbins + 1)
    return edges[:-1], edges[1:], vals


def read_bearing_track(path, chrom, a, b):
    pl = path.lower()
    if pl.endswith(".qcat.bgz") or pl.endswith(".qcat.gz"):
        return read_qcat(path, chrom, a, b)
    if pl.endswith(".bw") or pl.endswith(".bigwig"):
        return read_bigwig(path, chrom, a, b)
    return read_bedgraph(path, chrom, a, b)


# --------------------------------------------------------------- AgR annotation
def classify_segment(name):
    n = (name or "").lower()
    if re.search(r"trbv|ighv|igkv|trav|tcrv|[^a-z]v\d|^v\d", n):
        return "V", "#2c6fbb"
    if re.search(r"trbd|ighd|trad|[^a-z]d\d", n):
        return "D", "#1a9850"
    if re.search(r"trbj|ighj|igkj|traj|[^a-z]j\d", n):
        return "J", "#f08c00"
    if re.search(r"trbc|ighc|igkc|trac|constant|c1|c2|cmu|cgamma", n):
        return "C", "#7b3fa0"
    if re.search(r"enh|eb|emu|3'rr|3rr|cbe|ctcf|pax|silenc", n):
        return "reg", "#d62728"
    return "other", "#888888"


def read_annot(path, chrom, a, b):
    op = gzip.open if path.endswith(".gz") else open
    feats = []
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if f[0] != chrom:
                continue
            s, e = int(f[1]), int(f[2])
            if e <= a or s >= b:
                continue
            name = f[3] if len(f) > 3 else ""
            feats.append((s, e, name))
    return feats


# ------------------------------------------------------------------- plotting
def draw_annot(ax, feats, a, b):
    cats = {}
    for s, e, name in feats:
        cls, col = classify_segment(name)
        cats[cls] = col
        ax.add_patch(Rectangle((s, 0.15), max(e - s, (b - a) * 0.0008), 0.7,
                               facecolor=col, edgecolor="none"))
    ax.set_ylim(0, 1)
    ax.set_xlim(a, b)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.set_ylabel("AgR", rotation=0, ha="right", va="center", fontsize=9)
    if cats:
        handles = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="none")
                   for c in cats.values()]
        ax.legend(handles, list(cats.keys()), ncol=len(cats), fontsize=7,
                  loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False,
                  handlelength=1.0, columnspacing=1.0)


def draw_bearing(ax, track, a, b):
    s, e, v = track
    widths = e - s
    pos = v >= 0
    ax.bar(s[pos], v[pos], width=widths[pos], align="edge",
           color="#11806a", edgecolor="none")
    if (~pos).any():
        ax.bar(s[~pos], v[~pos], width=widths[~pos], align="edge",
               color="#b5670f", edgecolor="none")
    ax.axhline(0, color="#999999", lw=0.5)
    ax.set_xlim(a, b)
    ax.set_ylabel("BEARING\n1D", rotation=0, ha="right", va="center", fontsize=9)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)


def mb(x, _):
    return "%.2f" % (x / 1e6)


def nan_smooth(a, sigma):
    """nan-aware Gaussian smooth: smooths valid data and its support mask, so
    isolated masked pixels get filled but large masked regions stay masked."""
    from scipy.ndimage import gaussian_filter
    valid = np.isfinite(a)
    num = gaussian_filter(np.where(valid, a, 0.0), sigma=sigma, mode="nearest")
    den = gaussian_filter(valid.astype(float), sigma=sigma, mode="nearest")
    out = np.full(a.shape, np.nan, dtype=float)
    ok = den > 0.15
    out[ok] = num[ok] / den[ok]
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cool", help="single-condition cooler (single mode)")
    ap.add_argument("--cool-a", help="condition A cooler (ratio mode)")
    ap.add_argument("--cool-b", help="condition B cooler (ratio mode)")
    ap.add_argument("--region", required=True, help="chr6:40400000-42400000")
    ap.add_argument("--annot", default=None, help="AgR_mm10.bed")
    ap.add_argument("--bearing-track", default=None,
                    help="BEARING 1D: *.qcat.bgz | *.bedgraph | *.bw "
                         "(use a *_diff_*.qcat.bgz to match a ratio plot)")
    ap.add_argument("--vmax", type=float, default=2.0,
                    help="ratio color limit, log2 units (default 2 = 4x)")
    ap.add_argument("--min-count", type=int, default=0,
                    help="mask any pixel below this raw count in either "
                         "condition before ratio/density (try 5 to drop "
                         "1-vs-1 / 2-vs-1 sampling noise)")
    ap.add_argument("--smooth", type=float, default=0.0,
                    help="nan-aware Gaussian smoothing sigma in bins (e.g. 1.0) "
                         "to read domain-level change over salt-and-pepper; 0=off")
    ap.add_argument("--cmap", default=None)
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    chrom, a, b = parse_region(args.region)
    mode = "ratio" if (args.cool_a and args.cool_b) else "single"
    if mode == "single" and not args.cool:
        sys.exit("[ERROR] single mode needs --cool, ratio mode needs --cool-a and --cool-b")

    if mode == "single":
        edges, mat, w, starts, ends = fetch_region(args.cool, args.region)
        dens = mat / np.outer(w, w)
        dens[mat < max(args.min_count, 1)] = np.nan
        cmap = plt.get_cmap(args.cmap or "YlOrRd").copy()
        cmap.set_bad("#f4f4f4")
        finite = dens[np.isfinite(dens)]
        vmin = np.percentile(finite, 5) if finite.size else 1e-9
        vmax = np.percentile(finite, 99.5) if finite.size else 1.0
        norm = LogNorm(vmin=max(vmin, vmax * 1e-4), vmax=vmax)
        hm_label = "contact density  (counts / bp^2)"
    else:
        edges, matA, w, starts, ends = fetch_region(args.cool_a, args.region)
        _, matB, _, _, _ = fetch_region(args.cool_b, args.region)
        if matA.shape != matB.shape:
            sys.exit("[ERROR] A and B matrices differ in shape; not the same grid")
        densA = matA / np.outer(w, w)
        densB = matB / np.outer(w, w)
        thr = max(args.min_count, 1)
        mask = (matA < thr) | (matB < thr)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.log2(densA / densB)
        ratio[mask] = np.nan
        cmap = plt.get_cmap(args.cmap or "RdBu_r").copy()
        cmap.set_bad("#d9d9d9")
        norm = TwoSlopeNorm(vmin=-args.vmax, vcenter=0.0, vmax=args.vmax)
        mat = ratio
        an = os.path.basename(args.cool_a).split(".")[0]
        bn = os.path.basename(args.cool_b).split(".")[0]
        hm_label = "log2( %s / %s ) density" % (an, bn)

    # panel rows
    rows, hr = [], []
    if args.annot:
        rows.append("annot"); hr.append(0.5)
    if args.bearing_track:
        rows.append("bearing"); hr.append(1.1)
    rows.append("hm"); hr.append(7.0)

    fig = plt.figure(figsize=(8.5, 1.1 * sum(hr) + 1.2))
    gs = fig.add_gridspec(len(rows), 2, width_ratios=[40, 1],
                          height_ratios=hr, hspace=0.12, wspace=0.03)
    axes = {}
    base = None
    for i, r in enumerate(rows):
        ax = fig.add_subplot(gs[i, 0], sharex=base) if base is not None \
            else fig.add_subplot(gs[i, 0])
        base = ax if base is None else base
        axes[r] = ax

    if "annot" in axes:
        draw_annot(axes["annot"], read_annot(args.annot, chrom, a, b), a, b)
        plt.setp(axes["annot"].get_xticklabels(), visible=False)
    if "bearing" in axes:
        tr = read_bearing_track(args.bearing_track, chrom, a, b)
        if tr is None:
            axes["bearing"].text(0.5, 0.5, "no BEARING signal in region",
                                 ha="center", va="center",
                                 transform=axes["bearing"].transAxes, fontsize=8)
            axes["bearing"].set_xlim(a, b)
        else:
            draw_bearing(axes["bearing"], tr, a, b)
        plt.setp(axes["bearing"].get_xticklabels(), visible=False)

    hm = axes["hm"]
    grid = dens if mode == "single" else mat
    if args.smooth > 0:
        grid = nan_smooth(grid, args.smooth)
    pcm = hm.pcolormesh(edges, edges, grid, norm=norm, cmap=cmap, shading="flat")
    hm.set_xlim(a, b); hm.set_ylim(a, b)
    hm.invert_yaxis()
    hm.set_aspect("equal")
    hm.xaxis.set_major_formatter(plt.FuncFormatter(mb))
    hm.yaxis.set_major_formatter(plt.FuncFormatter(mb))
    hm.set_xlabel("%s position (Mb)" % chrom)
    cax = fig.add_subplot(gs[len(rows) - 1, 1])
    cb = fig.colorbar(pcm, cax=cax)
    cb.set_label(hm_label, fontsize=8)

    ttl = args.title or ("%s  %s" % (mode, args.region))
    fig.suptitle(ttl, fontsize=11, y=0.995)
    fig.savefig(args.out, dpi=170, bbox_inches="tight")
    sys.stderr.write("[plot] %s -> %s\n" % (mode, args.out))


if __name__ == "__main__":
    main()
