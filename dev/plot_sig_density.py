#!/usr/bin/env python3
"""Per-region x per-sample FDR-significant bin DENSITY heatmap.

Reads significant_bins_summary.tsv (region_name, region, <sample>_sig,
<sample>_total ...) and renders a heatmap of significant bins per 1000 scorable
bins, so region size does not dominate. Conditions are shown in canonical order
(DN, EbKO, DP, ProB, S3T3) with replicates side by side; regions are grouped by
locus category. Cells are annotated with the raw significant-bin count.

ASCII-only.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

# canonical condition order; replicates appended in file order
COND_ORDER = ["DN", "EbKO", "DP", "ProB", "S3T3"]
# canonical comparison order for diff-mode summaries
CMP_ORDER = ["DN_vs_EbKO", "DN_vs_DP", "DN_vs_ProB", "DN_vs_S3T3"]

# region grouping for readability (label -> region_names, in display order)
REGION_GROUPS = [
    ("Tcrb locus", ["wide", "new_wide", "tcrb", "vcluster", "v1_wide", "v1",
                     "tcrb_v7_v21", "tcrb_v11_v15", "tcrb_v21_try4",
                     "dj_cluster", "rc", "trcb_up"]),
    ("Ig loci", ["igk", "igh", "igh_wide"]),
    ("Recomb.", ["rag", "rag1_rag2_zoom"]),
    ("Housekeep.", ["gapdh", "actb", "emc7"]),
    ("Other", ["cd4", "cd8a", "wapl", "Col1a1"]),
]


def order_samples(columns):
    cols = [c[:-4] for c in columns if c.endswith("_sig")]
    is_diff = any("_vs_" in c for c in cols)
    if is_diff:
        ordered = [c for c in CMP_ORDER if c in cols]
        ordered.extend(sorted(c for c in cols if c not in ordered))
        return ordered, True
    ordered = []
    for cond in COND_ORDER:
        reps = sorted([s for s in cols if s.rsplit("_", 1)[0] == cond])
        ordered.extend(reps)
    ordered.extend([s for s in cols if s not in ordered])
    return ordered, False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per", type=float, default=1000.0,
                    help="Express density per N scorable bins (default 1000).")
    args = ap.parse_args()

    df = pd.read_csv(args.summary, sep="\t").set_index("region_name")
    samples, is_diff = order_samples(df.columns)

    # build ordered region list from groups (only those present), then any extras
    ordered_regions, group_bounds, group_labels = [], [], []
    for label, names in REGION_GROUPS:
        present = [n for n in names if n in df.index]
        if not present:
            continue
        group_bounds.append((len(ordered_regions), len(ordered_regions) + len(present)))
        group_labels.append(label)
        ordered_regions.extend(present)
    extras = [r for r in df.index if r not in ordered_regions]
    if extras:
        group_bounds.append((len(ordered_regions), len(ordered_regions) + len(extras)))
        group_labels.append("(unclassified)")
        ordered_regions.extend(extras)

    nrow, ncol = len(ordered_regions), len(samples)
    dens = np.zeros((nrow, ncol))
    sig = np.zeros((nrow, ncol), dtype=int)
    for i, reg in enumerate(ordered_regions):
        for j, s in enumerate(samples):
            ns = int(df.loc[reg, "%s_sig" % s])
            nt = int(df.loc[reg, "%s_total" % s])
            sig[i, j] = ns
            dens[i, j] = (args.per * ns / nt) if nt > 0 else 0.0

    # colormap: white -> deep blue (sequential); sqrt norm so the large AgR
    # loci (low density, high raw count) remain visible alongside small
    # high-density regions like rag1_rag2_zoom.
    pos = dens[dens > 0]
    vmax = np.percentile(pos, 98) if pos.size else 1.0
    vmax = max(vmax, 1.0)
    cmap = LinearSegmentedColormap.from_list(
        "wb", ["#ffffff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"])
    norm = PowerNorm(gamma=0.5, vmin=0, vmax=vmax)

    fig_h = 0.34 * nrow + 2.2
    fig_w = 0.62 * ncol + 4.6
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
    im = ax.imshow(dens, aspect="auto", cmap=cmap, norm=norm)

    # annotate raw significant counts
    for i in range(nrow):
        for j in range(ncol):
            if sig[i, j] > 0:
                shade = norm(dens[i, j])
                col = "white" if shade > 0.55 else "#222222"
                ax.text(j, i, str(sig[i, j]), ha="center", va="center",
                        fontsize=6, color=col)

    ax.set_xticks(range(ncol))
    ax.set_xticklabels(samples, rotation=90, fontsize=7)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(ordered_regions, fontsize=7)

    # condition separators + top labels (per-sample mode groups reps by
    # condition; diff mode shows each comparison as its own column).
    conds = [s.rsplit("_", 1)[0] for s in samples]
    if not is_diff:
        for j in range(1, ncol):
            if conds[j] != conds[j - 1]:
                ax.axvline(j - 0.5, color="black", lw=1.2)
        for cond in COND_ORDER:
            idx = [j for j, c in enumerate(conds) if c == cond]
            if idx:
                ax.text(np.mean(idx), -0.75, cond, ha="center", va="bottom",
                        fontsize=9, fontweight="bold")

    # region group separators + labels (placed in the LEFT margin, clear of grid)
    for (a, b), label in zip(group_bounds, group_labels):
        if a > 0:
            ax.axhline(a - 0.5, color="black", lw=1.2)
        ax.text(-3.2, (a + b - 1) / 2.0, label, rotation=90,
                ha="center", va="center", fontsize=7.5, color="#333333",
                fontweight="bold")

    ax.set_xlim(-0.5, ncol - 0.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cb.set_label("FDR-significant bins per %d scorable bins (sqrt scale)"
                 % int(args.per), fontsize=8)
    ax.set_title("BEARING significant-bin density by region and sample\n"
                 "(cell number = raw FDR-significant bin count)",
                 fontsize=9, pad=24)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print("Wrote %s (%d regions x %d samples)" % (args.out, nrow, ncol))


if __name__ == "__main__":
    main()
