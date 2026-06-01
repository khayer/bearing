#!/usr/bin/env python3
"""
oe_contact_grid.py

Render a fixed-scale observed/expected (O/E) contact-heatmap GRID across
conditions for one region. Closes the Fig 3B / Fig 4A layout gap: the analysis
(O/E values, self-contact) is done by hic/tcrb_contact_isolation.py; this script
only arranges per-condition O/E matrices into a tidy matrix with a shared color
scale and one colorbar.

O/E is computed exactly as in hic/tcrb_contact_isolation.py: balanced matrix
divided by the per-diagonal (distance) expected, then log2.

Usage:
  python3 oe_contact_grid.py \\
      --cool DN=DN_25kb.cool DP=DP_25kb.cool EbKO=... ProB=... S3T3=... \\
      --region chr6:40400000-42400000 --resolution 25000 \\
      --vmin -2 --vmax 2 --ncols 5 \\
      --annotate-oe-region chr6:41000000-41600000 \\
      --out fig3_panelB_oe_grid.png

ASCII-only.
"""

import argparse
import sys

import numpy as np


def parse_region(s):
    chrom, rest = s.split(":")
    start, end = rest.replace(",", "").split("-")
    return chrom, int(start), int(end)


def fetch_balanced(path, chrom, start, end, balance=True):
    import cooler
    c = cooler.Cooler(path)
    region = "%s:%d-%d" % (chrom, start, end)
    mat = c.matrix(balance=balance, sparse=False).fetch(region)
    return np.asarray(mat, dtype=float)


def oe_normalize(mat):
    """Divide each diagonal by its mean, then log2. NaN-safe."""
    n = mat.shape[0]
    oe = np.full_like(mat, np.nan)
    for k in range(n):
        diag_idx = (np.arange(n - k), np.arange(k, n))
        vals = mat[diag_idx]
        m = np.nanmean(vals)
        if m and np.isfinite(m) and m > 0:
            oe[diag_idx] = vals / m
            # mirror
            oe[(np.arange(k, n), np.arange(n - k))] = vals / m
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(oe)


def mean_oe_block(oe, region, sub, resolution, region_start):
    """Mean O/E (linear, not log2) within a sub-block, for annotation."""
    s = (sub[1] - region_start) // resolution
    e = (sub[2] - region_start) // resolution
    s = max(0, s); e = min(oe.shape[0], e)
    if e <= s:
        return np.nan
    block = oe[s:e, s:e]
    return float(np.nanmean(np.power(2.0, block)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cool", nargs="+", required=True,
                    help="COND=path.cool entries (order = grid order)")
    ap.add_argument("--region", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--vmin", type=float, default=-2.0)
    ap.add_argument("--vmax", type=float, default=2.0)
    ap.add_argument("--cmap", default="RdBu_r")
    ap.add_argument("--ncols", type=int, default=5)
    ap.add_argument("--no-balance", action="store_true")
    ap.add_argument("--annotate-oe-region", default=None,
                    help="sub-region chr:start-end; overlay its mean linear O/E")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    conds, paths = [], []
    for entry in args.cool:
        if "=" not in entry:
            sys.exit("ERROR: --cool entries must be COND=path, got %r" % entry)
        c, p = entry.split("=", 1)
        conds.append(c); paths.append(p)

    chrom, start, end = parse_region(args.region)
    sub = parse_region(args.annotate_oe_region) if args.annotate_oe_region else None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(conds)
    ncols = min(args.ncols, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 3.2 * nrows),
                             squeeze=False)
    im = None
    extent = [start, end, end, start]
    for i, (cond, path) in enumerate(zip(conds, paths)):
        ax = axes[i // ncols][i % ncols]
        try:
            mat = fetch_balanced(path, chrom, start, end,
                                 balance=not args.no_balance)
            oe = oe_normalize(mat)
        except Exception as e:
            ax.text(0.5, 0.5, "%s\n(read failed)\n%s" % (cond, e),
                    ha="center", va="center", fontsize=7, transform=ax.transAxes)
            ax.axis("off")
            continue
        im = ax.imshow(oe, cmap=args.cmap, vmin=args.vmin, vmax=args.vmax,
                       extent=extent, aspect="equal", interpolation="none")
        title = cond
        if sub is not None:
            moe = mean_oe_block(oe, (chrom, start, end), sub,
                                args.resolution, start)
            title = "%s  (O/E %.2f)" % (cond, moe)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])

    # hide unused axes
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    if im is not None:
        cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label("log2 O/E", fontsize=9)
    fig.suptitle("O/E contact, %s (%d kb)" % (args.region, args.resolution // 1000),
                 fontsize=11)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    sys.stderr.write("wrote %s\n" % args.out)


if __name__ == "__main__":
    main()
