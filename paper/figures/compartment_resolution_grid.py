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
import os
import sys

import numpy as np

# Reuse the canonical gene-density orientation from the compartment-analysis
# module so the two scripts can never drift. hic/ is not a package, so add the
# repo root to sys.path and import by module name.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_HIC_DIR = os.path.join(_REPO_ROOT, "hic")
if _HIC_DIR not in sys.path:
    sys.path.insert(0, _HIC_DIR)
try:
    from compartment_analysis import (load_gene_counts_per_bin,
                                       orient_pc1_by_gene_density)
    _HAVE_GENE_ORIENT = True
except Exception:
    _HAVE_GENE_ORIENT = False


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


def orient_pc1_genes(pc1, chrom, start, end, gtf_path):
    """Orient PC1 sign by gene density over the plotted region (canonical
    method, matching hic/compartment_analysis.py): correlate per-bin gene
    counts against PC1 and flip if the correlation is negative. Returns
    (oriented_pc1, rho) or (None, nan) if it could not be computed.
    """
    if not (_HAVE_GENE_ORIENT and gtf_path and os.path.exists(gtf_path)):
        return None, float("nan")
    nb = len(pc1)
    edges = np.linspace(start, end, nb + 1).astype(int)
    bin_starts = edges[:-1]
    bin_ends = edges[1:]
    try:
        gene_counts = load_gene_counts_per_bin(gtf_path, chrom,
                                               bin_starts, bin_ends)
        import pandas as pd
        df = pd.DataFrame({"value": pc1})
        flip, rho = orient_pc1_by_gene_density(df, np.asarray(gene_counts))
    except Exception:
        return None, float("nan")
    return (-pc1 if flip else pc1), rho


def orient_pc1(pc1, x, start, end, chrom, landmarks):
    """Fallback orientation: flip the PC1 sign so the A-labelled landmarks come
    out positive (red). Used only when gene-density orientation is unavailable
    (no GTF). The compartment eigenvector sign is arbitrary per computation;
    this anchors it to known biology.

    Returns the (possibly sign-flipped) PC1 array.
    """
    if not landmarks:
        return pc1
    votes = 0.0
    for _, lchrom, lpos, ab in landmarks:
        if lchrom != chrom or not (start <= lpos <= end):
            continue
        idx = int(np.clip(round((lpos - start) / (end - start) * (len(x) - 1)),
                          0, len(x) - 1))
        v = pc1[idx]
        if v is None or np.isnan(v):
            continue
        want_pos = (ab.upper() == "A")
        if want_pos:
            votes += 1.0 if v > 0 else -1.0
        else:
            votes += 1.0 if v < 0 else -1.0
    return -pc1 if votes < 0 else pc1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pc1", nargs="+", required=True,
                    help="COND=path entries; path may contain {res}")
    ap.add_argument("--resolutions", nargs="+", type=int, required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--landmarks", nargs="*", default=[],
                    help="NAME:chrom:pos:AB landmark annotations")
    ap.add_argument("--nbins", type=int, default=600)
    ap.add_argument("--gtf", default=None,
                    help="GTF for gene-density PC1 orientation (canonical "
                         "method; matches hic/compartment_analysis.py). If "
                         "omitted, falls back to --landmarks anchoring.")
    ap.add_argument("--orient-region", default=None,
                    help="Region used to DECIDE the sign flip via gene density "
                         "(default: whole chromosome of --region). Orientation "
                         "is more stable over a wide region; the decided sign "
                         "is then applied to the plotted --region.")
    ap.add_argument("--orient-nbins", type=int, default=2000,
                    help="Bins for the orientation-decision region.")
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

    # Layout: rows = conditions, columns = resolutions. All panels share the
    # same genomic x-axis, so stacking conditions vertically aligns them and a
    # single Tcrb gridline reads straight down every condition.
    nrows = len(conds)
    ncols = len(args.resolutions)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(2.8 * ncols, 1.5 * nrows),
                             squeeze=False, sharex=True)
    x = np.linspace(start, end, args.nbins)

    # Tcrb gets a distinct marker; other landmarks stay faint.
    def is_tcrb(name):
        return name.lower().startswith("tcrb") or name.lower().startswith("trcb")

    # Orientation-decision region: prefer a wide window (whole chromosome span
    # if not given) so the gene-density correlation is stable; the resulting
    # per-(cond,res) sign is then applied to the plotted region.
    if args.orient_region:
        ochrom, ostart, oend = parse_region(args.orient_region)
    else:
        ochrom, ostart, oend = chrom, start, end

    def decide_flip(tmpl, res):
        """Return -1.0 or +1.0 to multiply the plotted PC1 by, decided over the
        wide orient-region. Gene density preferred; landmark anchor fallback."""
        path = tmpl.replace("{res}", str(res))
        try:
            opc1 = read_bigwig_binned(path, ochrom, ostart, oend,
                                      args.orient_nbins)
        except Exception:
            return 1.0
        if args.gtf:
            oriented, rho = orient_pc1_genes(opc1, ochrom, ostart, oend,
                                             args.gtf)
            if oriented is not None:
                # oriented == -opc1 when a flip happened
                return -1.0 if (np.nansum(np.abs(oriented + opc1)) <
                                np.nansum(np.abs(oriented - opc1))) else 1.0
        # landmark fallback over the orient-region
        ox = np.linspace(ostart, oend, len(opc1))
        anchored = orient_pc1(opc1, ox, ostart, oend, ochrom, landmarks)
        return -1.0 if (np.nansum(np.abs(anchored + opc1)) <
                        np.nansum(np.abs(anchored - opc1))) else 1.0

    for ri, (cond, tmpl) in enumerate(zip(conds, tmpls)):
        for ci, res in enumerate(args.resolutions):
            ax = axes[ri][ci]
            path = tmpl.replace("{res}", str(res))
            try:
                pc1 = read_bigwig_binned(path, chrom, start, end, args.nbins)
            except Exception as e:
                ax.text(0.5, 0.5, "read failed\n%s" % e, ha="center",
                        va="center", fontsize=6, transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
                continue
            # Anchor the eigenvector sign before plotting. The flip is decided
            # over the wide orient-region (gene density preferred), then applied
            # here so all panels share a consistent A=positive convention.
            pc1 = pc1 * decide_flip(tmpl, res)
            pos = np.where(pc1 > 0, pc1, 0.0)
            neg = np.where(pc1 < 0, pc1, 0.0)
            ax.fill_between(x, 0, pos, color="#c0392b", linewidth=0)  # A
            ax.fill_between(x, 0, neg, color="#0b3d91", linewidth=0)  # B
            ax.axhline(0, color="#888", linewidth=0.5)
            for name, lchrom, lpos, ab in landmarks:
                if lchrom == chrom and start <= lpos <= end:
                    if is_tcrb(name):
                        ax.axvline(lpos, color="#111", linewidth=1.4,
                                   linestyle="-", zorder=5)
                    else:
                        ax.axvline(lpos, color="#999", linewidth=0.5,
                                   linestyle=":")
            # Label the Tcrb line once, on the top row only.
            if ri == 0:
                ax.set_title("%d kb" % (res // 1000), fontsize=9)
                for name, lchrom, lpos, ab in landmarks:
                    if is_tcrb(name) and lchrom == chrom and start <= lpos <= end:
                        ax.annotate("Tcrb", xy=(lpos, 0.96),
                                    xycoords=("data", "axes fraction"),
                                    xytext=(2, 0), textcoords="offset points",
                                    ha="left", va="top", fontsize=7,
                                    color="#111", zorder=6)
            if ci == 0:
                ax.set_ylabel(cond, fontsize=9)
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
