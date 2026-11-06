#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : make_figS12_elemclass.py
# Version  : 1.0.0
# Date     : 2026-11-06
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Supplementary Figure S12 (element-resolution regional enrichment at Tcrb)
# from the per-class outputs of the Snakemake rule element_class_enrichment
# (regional/elemclass_<class>_<cmp>.tsv). Joint Benjamini-Hochberg over every
# element x comparison row (165 x 4 = 660 in the manuscript), as in
# consolidate_regional_enrichment.py --joint-bh.
#   (A) number of elements per class with joint q < 0.05, per comparison
#   (B) every element significant in at least one comparison: -log10 joint q
#       (capped at 20), direction (^ DN-enriched, v other-enriched) and
#       k/n (below-threshold bins / scored bins); outline = q < 0.05
# Also writes <out>.source.tsv (all rows with q_joint) and prints a summary
# for the legend. Usage:
#   python make_figS12_elemclass.py --regional-dir results_prior1_n100_l1/regional \
#       --out S12_elemclass_l1
# ASCII only.
import argparse
import csv
import glob
import math
import os
import re

CLASSES = ["CBE", "RSS", "Vbeta_functional", "Vbeta_pseudogene", "promoter", "lncRNA",
           "DJC_recomb_center"]
CLABEL = {"CBE": "CBE (CTCF)", "RSS": "RSS", "Vbeta_functional": "Vb functional",
          "Vbeta_pseudogene": "Vb pseudogene", "promoter": "promoter", "lncRNA": "lncRNA",
          "DJC_recomb_center": "Db-Jb-Cb / RC"}
CMPS = ["DN_vs_DP", "DN_vs_EbKO", "DN_vs_ProB", "DN_vs_S3T3"]
CMPL = {"DN_vs_DP": "vs DP", "DN_vs_EbKO": "vs EbKO", "DN_vs_ProB": "vs Pro-B",
        "DN_vs_S3T3": "vs 3T3"}


def bh(p):
    n = len(p)
    o = sorted(range(n), key=lambda i: p[i])
    q = [0.0] * n
    prev = 1.0
    for r, i in enumerate(reversed(o)):
        k = n - r
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regional-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    a = ap.parse_args()
    rows = []
    for f in sorted(glob.glob(os.path.join(a.regional_dir, "elemclass_*_DN_vs_*.tsv"))):
        m = re.match(r"elemclass_(.+)_(DN_vs_[A-Za-z0-9]+)\.tsv$", os.path.basename(f))
        if not m or m.group(1) not in CLASSES or m.group(2) not in CMPS:
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            r["element_class"] = m.group(1)
            r["comparison"] = m.group(2)
            rows.append(r)
    if not rows:
        raise SystemExit("no elemclass_<class>_<cmp>.tsv files in %s" % a.regional_dir)
    q = bh([float(r["p_combined"]) for r in rows])
    for r, v in zip(rows, q):
        r["q_joint"] = v
    print("rows (element x comparison): %d" % len(rows))
    hdr = ["element_class", "region_name", "comparison", "L_region_bins", "k", "k_pos", "k_neg",
           "p_spatial", "p_directional", "p_combined", "q_joint", "g_dir"]
    with open(a.out + ".source.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(h, "")) for h in hdr) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    sig = [r for r in rows if r["q_joint"] < a.alpha]
    fig = plt.figure(figsize=(12, 4.6))
    axA = fig.add_axes([0.11, 0.18, 0.30, 0.70])
    cnt = [[sum(1 for r in sig if r["element_class"] == c and r["comparison"] == k)
            for k in CMPS] for c in CLASSES]
    vmax = max(1, max(max(x) for x in cnt))
    axA.imshow(cnt, cmap="Blues", vmin=0, vmax=vmax * 1.3, aspect="auto")
    for i, c in enumerate(CLASSES):
        for j in range(len(CMPS)):
            axA.text(j, i, str(cnt[i][j]), ha="center", va="center", fontsize=9,
                     color="white" if cnt[i][j] > vmax * 0.6 else "black")
            if cnt[i][j]:
                axA.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#e75480", lw=2))
    axA.set_yticks(range(len(CLASSES)))
    axA.set_yticklabels([CLABEL[c] for c in CLASSES], fontsize=8)
    axA.set_xticks(range(len(CMPS)))
    axA.set_xticklabels([CMPL[k] for k in CMPS], fontsize=8)
    axA.set_title("A  Significant elements per class (joint q < %.2f)" % a.alpha, fontsize=9,
                  loc="left")

    axB = fig.add_axes([0.58, 0.18, 0.34, 0.70])
    names = []
    for r in sig:
        if r["region_name"] not in names:
            names.append(r["region_name"])
    if not names:
        axB.text(0.5, 0.5, "no element reaches joint q < %.2f" % a.alpha, ha="center")
        axB.axis("off")
    else:
        idx = {(r["region_name"], r["comparison"]): r for r in rows}
        mat = []
        for n in names:
            line = []
            for k in CMPS:
                r = idx.get((n, k))
                line.append(min(20.0, -math.log10(max(r["q_joint"], 1e-300))) if r else 0.0)
            mat.append(line)
        axB.imshow(mat, cmap="Blues", vmin=0, vmax=20, aspect="auto")
        for i, n in enumerate(names):
            for j, k in enumerate(CMPS):
                r = idx.get((n, k))
                if r is None:
                    continue
                if r["q_joint"] < a.alpha:
                    kp, kn = int(float(r["k_pos"])), int(float(r["k_neg"]))
                    arrow = "^" if kp >= kn else "v"
                    txt = "%.0e\n%s %s/%s" % (r["q_joint"], arrow, r["k"], r["L_region_bins"])
                    axB.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#e75480",
                                            lw=2))
                else:
                    txt = "n.s."
                axB.text(j, i, txt, ha="center", va="center", fontsize=7,
                         color="white" if mat[i][j] > 12 else "black")
        axB.set_yticks(range(len(names)))
        axB.set_yticklabels(names, fontsize=8)
        axB.set_xticks(range(len(CMPS)))
        axB.set_xticklabels([CMPL[k] for k in CMPS], fontsize=8)
    axB.set_title("B  Significant elements: -log10 joint q (cap 20), direction, k/n", fontsize=9,
                  loc="left")
    fig.text(0.11, 0.03, "^ = DN-enriched (k_pos >= k_neg); k/n = below-threshold bins / scored "
             "bins; joint BH over %d element x comparison tests." % len(rows), fontsize=7.5)
    fig.savefig(a.out + ".pdf")
    fig.savefig(a.out + ".png", dpi=250)
    print("significant element x comparison pairs: %d" % len(sig))
    for r in sorted(sig, key=lambda x: x["q_joint"]):
        print("  %-28s %-11s q=%.2g  k=%s/%s  k_pos=%s k_neg=%s  class=%s" % (
            r["region_name"], r["comparison"], r["q_joint"], r["k"], r["L_region_bins"],
            r["k_pos"], r["k_neg"], r["element_class"]))
    print("wrote %s.pdf/.png/.source.tsv" % a.out)


if __name__ == "__main__":
    main()
