#!/usr/bin/env python3
"""
sanity_check_tad_boundaries.py

Sanity-check TAD boundary calls and insulation values at two disputed
positions on chr6 in the BEARING wide-locus analysis.

Disputed claim 1: DP lacks the chr6:41,562,500 boundary (Eb-Trbv31 region)
                  that is present in DN/EbKO/ProB/S3T3.
Disputed claim 2: chr6:41,275,000 boundary is present in DN/EbKO/ProB
                  but absent in DP/S3T3 (Vb-cluster-internal).

For each disputed position the script reports:
  (a) every boundary call within +/-100 kb in each of the 5 conditions
  (b) the raw insulation score at the disputed position plus neighbors
      at +/-25, 50, 75, 100 kb
  (c) local maxima (boundary-like peaks) within +/-150 kb that may have
      been sub-threshold in the original boundary call
  (d) a small zoom-in PDF of the insulation tracks for visual check

Adjust BOUNDARY_PATTERN and INSULATION_PATTERN to your file naming.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Configure here
# ---------------------------------------------------------------------
DEFAULT_TAD_DIR = (
    "/mnt/isilon/bassing_lab/projects/HiC_V31_NT_allele/"
    "data/endpoints/HiC_explorer_mm10/07tad"
)
CONDITIONS = ["DN", "DP", "EbKO", "ProB", "S3T3"]

# File-name patterns - edit if your naming differs
BOUNDARY_PATTERN = "{cond}_boundaries.bed"
INSULATION_PATTERN = "{cond}_tad_score.bm"

# Which column of the .bm file holds the insulation score to use.
# HiCExplorer .bm is bedgraph-like with one score column per window size.
# 0-based index into the score columns (i.e. 4th column = index 0).
INSULATION_SCORE_COL = 0

# Disputed positions: (chrom, position_bp, label)
DISPUTED = [
    ("chr6", 41562500, "Eb_Trbv31_boundary"),
    ("chr6", 41275000, "Vb_cluster_internal"),
]

# Tolerance values to test (bp)
TOLERANCES = [25000, 50000, 75000]


# ---------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------
def read_boundaries(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df = df.iloc[:, :3].copy()
    df.columns = ["chrom", "start", "end"]
    df["mid"] = (df["start"] + df["end"]) // 2
    return df


def read_insulation(path, score_col=0):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    # cols: chrom, start, end, score_1, [score_2, ...]
    score_idx = 3 + score_col
    out = df.iloc[:, [0, 1, 2, score_idx]].copy()
    out.columns = ["chrom", "start", "end", "score"]
    out["mid"] = (out["start"] + out["end"]) // 2
    return out


# ---------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------
def boundaries_near(df, chrom, pos, window):
    sub = df[(df["chrom"] == chrom)
             & (df["mid"] >= pos - window)
             & (df["mid"] <= pos + window)]
    return sub.sort_values("mid")


def insulation_profile(df, chrom, pos, window=100000, step=25000):
    offsets = list(range(-window, window + 1, step))
    out = []
    for off in offsets:
        target = pos + off
        sub = df[(df["chrom"] == chrom)
                 & (df["start"] <= target)
                 & (df["end"] > target)]
        score = float(sub.iloc[0]["score"]) if len(sub) else np.nan
        out.append((off, target, score))
    return out


def local_extrema(df, chrom, center, half_width=150000):
    sub = df[(df["chrom"] == chrom)
             & (df["mid"] >= center - half_width)
             & (df["mid"] <= center + half_width)].copy()
    sub = sub.sort_values("mid").reset_index(drop=True)
    s = sub["score"].values
    found = []
    for i in range(1, len(s) - 1):
        if np.isnan(s[i-1]) or np.isnan(s[i]) or np.isnan(s[i+1]):
            continue
        if s[i] > s[i-1] and s[i] > s[i+1]:
            found.append(("max", int(sub.iloc[i]["mid"]), float(s[i])))
        elif s[i] < s[i-1] and s[i] < s[i+1]:
            found.append(("min", int(sub.iloc[i]["mid"]), float(s[i])))
    return found


def classify_at_tolerance(bounds_by_cond, chrom, pos, tol):
    """Return per-condition (present/absent, nearest_offset)."""
    result = {}
    for c, df in bounds_by_cond.items():
        near = boundaries_near(df, chrom, pos, tol)
        if len(near) == 0:
            result[c] = ("absent", None)
        else:
            offs = (near["mid"] - pos).abs()
            i = offs.idxmin()
            result[c] = ("present", int(near.loc[i, "mid"]) - pos)
    return result


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_zoom(insul_by_cond, bounds_by_cond, chrom, center, half_width,
              out_pdf, title):
    fig, axes = plt.subplots(len(CONDITIONS), 1,
                             figsize=(9, 1.6 * len(CONDITIONS)),
                             sharex=True)
    if len(CONDITIONS) == 1:
        axes = [axes]
    lo = center - half_width
    hi = center + half_width
    for ax, c in zip(axes, CONDITIONS):
        ins = insul_by_cond[c]
        sub = ins[(ins["chrom"] == chrom)
                  & (ins["end"] >= lo) & (ins["start"] <= hi)]
        ax.plot(sub["mid"], sub["score"], color="black", lw=1.2)
        ax.fill_between(sub["mid"], 0, sub["score"], alpha=0.18)
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.axvline(center, color="red", lw=1.0, alpha=0.7)
        bnd = bounds_by_cond[c]
        bsub = bnd[(bnd["chrom"] == chrom)
                   & (bnd["mid"] >= lo) & (bnd["mid"] <= hi)]
        for _, b in bsub.iterrows():
            ax.axvline(b["mid"], color="blue", lw=0.8, alpha=0.5)
        ax.set_ylabel(c, rotation=0, ha="right", va="center")
        ax.set_xlim(lo, hi)
    axes[-1].set_xlabel("{} position (bp)".format(chrom))
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tad-dir", default=DEFAULT_TAD_DIR)
    ap.add_argument("--out-tsv", default="tad_sanity_report.tsv")
    ap.add_argument("--out-pdf-prefix", default="tad_sanity_zoom")
    ap.add_argument("--zoom-half-width", type=int, default=300000,
                    help="bp on each side of disputed pos in zoom PDF")
    args = ap.parse_args()

    # Load
    bounds = {}
    insul = {}
    for c in CONDITIONS:
        bp = os.path.join(args.tad_dir, BOUNDARY_PATTERN.format(cond=c))
        ip = os.path.join(args.tad_dir, INSULATION_PATTERN.format(cond=c))
        for p in (bp, ip):
            if not os.path.exists(p):
                sys.exit("MISSING: {}".format(p))
        bounds[c] = read_boundaries(bp)
        insul[c] = read_insulation(ip, INSULATION_SCORE_COL)
        print("Loaded {:<6} boundaries={} insulation_bins={}"
              .format(c, len(bounds[c]), len(insul[c])))

    rows = []
    for chrom, pos, label in DISPUTED:
        print()
        print("=" * 72)
        print("DISPUTED: {} at {} {:,}".format(label, chrom, pos))
        print("=" * 72)

        # (a) boundaries within +/-100 kb
        print()
        print("(a) Boundaries within +/-100 kb of disputed position")
        for c in CONDITIONS:
            near = boundaries_near(bounds[c], chrom, pos, 100000)
            if len(near) == 0:
                print("  {:<6} NONE within +/-100 kb".format(c))
                rows.append(dict(label=label, condition=c,
                                 type="boundary_call",
                                 offset_bp=None, abs_pos=None,
                                 score=None, note="no boundary +/-100 kb"))
            else:
                for _, b in near.iterrows():
                    off = int(b["mid"]) - pos
                    print("  {:<6} boundary at {:>12,} (offset {:+8,} bp)"
                          .format(c, int(b["mid"]), off))
                    rows.append(dict(label=label, condition=c,
                                     type="boundary_call",
                                     offset_bp=off, abs_pos=int(b["mid"]),
                                     score=None, note=""))

        # (b) tolerance sensitivity
        print()
        print("(b) Classification under different tolerances")
        header = "  cond   " + "  ".join(
            ["{:>5}kb".format(t // 1000) for t in TOLERANCES])
        print(header)
        for c in CONDITIONS:
            cells = []
            for tol in TOLERANCES:
                res = classify_at_tolerance({c: bounds[c]}, chrom, pos, tol)
                state, off = res[c]
                if state == "absent":
                    cells.append("absent ")
                else:
                    cells.append("{:+5,}".format(off))
            print("  {:<6}".format(c), "  ".join("{:>7}".format(x)
                                                 for x in cells))

        # (c) insulation profile at +/-100 kb (25 kb steps)
        print()
        print("(c) Insulation score profile at +/-100 kb (25 kb steps)")
        steps = list(range(-100000, 100001, 25000))
        hdr = "  cond   " + " ".join("{:>+5}kb".format(s // 1000)
                                      for s in steps)
        print(hdr)
        for c in CONDITIONS:
            prof = insulation_profile(insul[c], chrom, pos, 100000, 25000)
            cells = []
            for off, target, score in prof:
                cells.append(" {:>7.3f}".format(score)
                             if not np.isnan(score) else "    NA  ")
                rows.append(dict(label=label, condition=c,
                                 type="insulation",
                                 offset_bp=off, abs_pos=target,
                                 score=None if np.isnan(score) else score,
                                 note=""))
            print("  {:<6}".format(c), " ".join(cells))

        # (d) local extrema (sub-threshold boundary candidates)
        print()
        print("(d) Local extrema in insulation track within +/-150 kb")
        for c in CONDITIONS:
            ext = local_extrema(insul[c], chrom, pos, 150000)
            if not ext:
                print("  {:<6} no extrema found".format(c))
                continue
            for kind, mid, score in ext:
                off = mid - pos
                tag = "MAX(boundary?)" if kind == "max" else "MIN(TAD core)"
                print("  {:<6} {} at {:>12,} (offset {:+8,} bp) "
                      "score={:.3f}".format(c, tag, mid, off, score))
                rows.append(dict(label=label, condition=c,
                                 type="extremum_" + kind,
                                 offset_bp=off, abs_pos=mid,
                                 score=score, note=""))

        # (e) zoom PDF
        out_pdf = "{}_{}.pdf".format(args.out_pdf_prefix, label)
        plot_zoom(insul, bounds, chrom, pos, args.zoom_half_width,
                  out_pdf,
                  "{}  {}:{:,}  (red=disputed, blue=boundary calls)"
                  .format(label, chrom, pos))
        print()
        print("  zoom PDF -> {}".format(out_pdf))

    out = pd.DataFrame(rows)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print()
    print("Report TSV -> {}".format(args.out_tsv))


if __name__ == "__main__":
    main()
