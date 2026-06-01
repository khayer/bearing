#!/usr/bin/env python3
"""
compartment_resolution_grid.py

Render an A/B compartment-eigenvector GRID: rows = resolutions, columns =
conditions, one PC1 track per cell, over one region. Closes the Supp S7B layout
gap (hic/compartment_analysis.py emits per-resolution figures; this arranges the
resolution x condition composite from the PC1 bigwigs directly).

PC1 sign convention matches hic/compartment_analysis.py: positive = A (red),
negative = B (blue). Optionally orient by gene-density landmarks passed as
--landmarks so A/B is consistent across cells.

Usage:
  python3 compartment_resolution_grid.py \\
      --pc1 DN=DN_pc1_{res}.bw DP=DP_pc1_{res}.bw EbKO=... ProB=... S3T3=... \\
      --resolutions 25000 50000 100000 \\
      --region chr6:39000000-43500000 \\
      --landmarks Cntnap2:chr6:46400000:B Gpnmb:chr6:48900000:A Tcrb:chr6:41000000:A \\
      --out suppS7_panelB_compartment_res.png

Each --pc1 path may contain a literal {res} placeholder substituted per row.
ASCII-only.
"""

import argparse
import sys

import numpy as np


def parse_region(s):
    chrom, rest = s.split(":")
    start, end = rest.replace(",", "").split("-")
    return chrom, int(start), int(end)


def read_bigwig_binned(path, chrom, start, end, nbins=600):
    import pyBigWig
    bw = pyBigWig.open(path)
    try:
        vals = bw.stats(chrom, start, end, type="mean", nBins=nbins)
    finally:
        bw.close()
    return np.array([np.nan if v is None else float(v) for v in vals])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pc1", nargs="+", required=True,
                    help="COND=path entries; path may contain {res}")
    ap.add_argument("--resolutions", nargs="+", type=int, required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--landmarks", nargs="*", default=[],
                    help="NAME:chrom:pos:AB landmark annotations")
    ap.add_argument("--nbins", type=int, default=600)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    conds, tmpls = [], []
    for entry in args.pc1:
        c, p = entry.split("=", 1)
        conds.append(c); tmpls.append(p)
    chrom, start, end = parse_region(args.region)

    landmarks = []
    for lm in args.landmarks:
        parts = lm.split(":")
        if len(parts) >= 4:
            landmarks.append((parts[0], parts[1], int(parts[2]), parts[3]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nrows = len(args.resolutions)
    ncols = len(conds)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(2.6 * ncols, 1.8 * nrows),
                             squeeze=False, sharex=True)
    x = np.linspace(start, end, args.nbins)

    for ri, res in enumerate(args.resolutions):
        for ci, (cond, tmpl) in enumerate(zip(conds, tmpls)):
            ax = axes[ri][ci]
            path = tmpl.replace("{res}", str(res))
            try:
                pc1 = read_bigwig_binned(path, chrom, start, end, args.nbins)
            except Exception as e:
                ax.text(0.5, 0.5, "read failed\n%s" % e, ha="center",
                        va="center", fontsize=6, transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
                continue
            pos = np.where(pc1 > 0, pc1, 0.0)
            neg = np.where(pc1 < 0, pc1, 0.0)
            ax.fill_between(x, 0, pos, color="#c0392b", linewidth=0)  # A
            ax.fill_between(x, 0, neg, color="#0b3d91", linewidth=0)  # B
            ax.axhline(0, color="#888", linewidth=0.5)
            for _, lchrom, lpos, ab in landmarks:
                if lchrom == chrom and start <= lpos <= end:
                    ax.axvline(lpos, color="#444", linewidth=0.6, linestyle=":")
            if ri == 0:
                ax.set_title(cond, fontsize=9)
            if ci == 0:
                ax.set_ylabel("%d kb" % (res // 1000), fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])

    # landmark legend row labels along the bottom
    if landmarks:
        names = ", ".join("%s=%s" % (n, ab) for n, _, _, ab in landmarks)
        fig.text(0.5, 0.005, "landmarks: " + names, ha="center", fontsize=8)
    fig.suptitle("A (red) / B (blue) compartments, %s" % args.region, fontsize=11)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    sys.stderr.write("wrote %s\n" % args.out)


if __name__ == "__main__":
    main()
