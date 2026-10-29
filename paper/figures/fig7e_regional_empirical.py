#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : paper/figures/fig7e_regional_empirical.py
# Version  : 1.1.0
# Date     : 2026-10-28
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
#
# 1.1.0 (2026-10-28): q shown with 2 significant digits (0.0075, not 0.007);
#   a q below --alpha that would round up to it is printed with more digits
#   (0.04997, not 0.05); emp p at the null floor (emp_count 0) printed as
#   "<= floor"; emp p < 0.01 with 2 significant digits (0.0030); less empty space above and between the panels.
"""
Figure 7E: region x comparison heatmap of the CALIBRATED regional test.

Replaces the panel_E_heatmap_<locus>.pdf of consolidate_regional_enrichment.py,
which colours the NOMINAL q_combined (not usable as a regional call; see the
contrast-matched calibration, Methods M.7 / Table S11). This panel reads the
regional_empirical_<locus>.tsv tables (dev/regional_empirical_table.py) and
shows, per cell:

  colour  -log10(q_testable)   (BH over testable regions; primary family)
  text    emp p, q_testable, q_all, and k_pos/k_neg (bins at p < threshold,
          first-named condition higher / second higher)
  outline q_testable < --alpha (default 0.05): a regional call
  hatched not testable (fewer than --min-bins tested bins; emp p undefined)

Input  : one or more regional_empirical_<locus>.tsv (one panel per file)
Output : --out-prefix.{pdf,png} and --out-prefix.source.tsv (the plotted cells)

Style  : thesis/manuscript house look (purple on grey neutrals, Arial).
ASCII only. Reads real tables; fabricates nothing.
"""

import argparse
import csv
import math
import os
import sys

P_DARK, P_MID, P_LIGHT, P_PALE = "#4B2E83", "#8A6FBF", "#C9B8E6", "#EFE9F7"
G_DARK, G_MID, G_LIGHT, G_PALE = "#3A3A3A", "#7A7A7A", "#C8C8C8", "#F2F2F2"
INK = "#222222"

REGION_LABELS = {
    "Trbv1": "Trbv1",
    "Vbeta_cluster_main": "V-beta cluster",
    "NEGCTL_inter_cluster_silent": "Neg. control (silent)",
    "DJ_RC": "DJ recombination center",
    "Ebeta_Trbv31_CBE3": "E-beta / Trbv31 / CBE3",
    "NEGCTL_3flanking_silent": "Neg. control (3' flank)",
    "CH_JH_iEmu": "CH / JH / iEmu",
    "DH_cluster": "DH cluster",
    "VH_proximal": "VH proximal",
    "VH_distal": "VH distal",
}
LOCUS_LABELS = {"tcrb": "Tcrb", "igh": "Igh"}


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fmt_p(x):
    if x is None:
        return "-"
    if x < 0.001:
        return "%.1e" % x
    if x < 0.01:
        return "%#.2g" % x
    return "%.2f" % x


def fmt_q(x, alpha=0.05):
    """q value, 2 significant digits; never print a value below alpha as alpha
    or above (a call must not read as q 0.05)."""
    if x is None:
        return "-"
    if x < 0.001:
        return "%.1e" % x
    if x >= 0.1:
        return "%.2f" % x
    s = "%#.2g" % x
    if x < alpha <= float(s):
        s = "%#.4g" % x
    return s


def fmt_emp(x, count):
    """emp p; at the floor (0 null draws as extreme) print it as an upper bound."""
    s = fmt_p(x)
    try:
        if int(float(count)) == 0:
            return "<= " + s
    except (TypeError, ValueError):
        pass
    return s


def cmp_label(c, ref_first=True):
    a, _, b = c.partition("_vs_")
    b = {"S3T3": "3T3", "ProB": "Pro-B", "EbKO": "E-beta KO"}.get(b, b)
    a = {"S3T3": "3T3", "ProB": "Pro-B", "EbKO": "E-beta KO"}.get(a, a)
    return "%s vs\n%s" % (a, b)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", nargs="+", required=True,
                    help="regional_empirical_<locus>.tsv files, one panel each")
    ap.add_argument("--comparisons", nargs="*", default=None,
                    help="column order (default: order of first appearance)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--vmax", type=float, default=3.0,
                    help="colour scale cap for -log10(q_testable) (default 3)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle

    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                                             "DejaVu Sans"],
                         "font.size": 9, "text.color": INK,
                         "axes.labelcolor": INK, "xtick.color": INK,
                         "ytick.color": INK})
    cmap = LinearSegmentedColormap.from_list(
        "bearing_purple", ["#FFFFFF", P_PALE, P_LIGHT, P_MID, P_DARK])

    panels = []
    for path in a.tables:
        rows = read_rows(path)
        if not rows:
            sys.exit("empty table: %s" % path)
        locus = rows[0].get("locus") or os.path.basename(path)
        cmps = a.comparisons or []
        if not cmps:
            for r in rows:
                if r["comparison"] not in cmps:
                    cmps.append(r["comparison"])
        regions = []
        for r in rows:
            if r["region"] not in regions:
                regions.append(r["region"])
        cell = {(r["region"], r["comparison"]): r for r in rows}
        panels.append((locus, regions, list(cmps), cell))

    heights = [max(2.2, 0.62 * len(p[1]) + 1.2) for p in panels]
    ncols = max(len(p[2]) for p in panels)
    fig_w = 2.3 + 1.55 * ncols
    fig_h = sum(heights) + 1.3
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(len(panels), 1, height_ratios=heights, hspace=0.22)

    src = []
    im = None
    for pi, (locus, regions, cmps, cell) in enumerate(panels):
        ax = fig.add_subplot(gs[pi, 0])
        nr, nc = len(regions), len(cmps)
        ax.set_xlim(0, nc)
        ax.set_ylim(nr, 0)
        for i, reg in enumerate(regions):
            for j, c in enumerate(cmps):
                r = cell.get((reg, c))
                if r is None:
                    ax.add_patch(Rectangle((j, i), 1, 1, facecolor="white",
                                           edgecolor=G_LIGHT, lw=0.8))
                    ax.text(j + 0.5, i + 0.5, "absent", ha="center", va="center",
                            fontsize=7, color=G_MID)
                    continue
                testable = r.get("testable", "yes") == "yes"
                emp = fnum(r.get("emp_p"))
                qt = fnum(r.get("q_testable"))
                qa = fnum(r.get("q_all"))
                kp, kn = r.get("k_pos", ""), r.get("k_neg", "")
                nb = r.get("tested_bins", "")
                if not testable or emp is None:
                    ax.add_patch(Rectangle((j, i), 1, 1, facecolor=G_PALE,
                                           edgecolor=G_LIGHT, lw=0.8, hatch="////"))
                    ax.text(j + 0.5, i + 0.5, "not testable\n(%s bins)" % nb,
                            ha="center", va="center", fontsize=7, color=G_DARK)
                    src.append([locus, reg, c, nb, "no", "", "", "", kp, kn])
                    continue
                v = min(-math.log10(max(qt, 1e-300)), a.vmax)
                col = cmap(v / a.vmax)
                call = qt < a.alpha
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor=col,
                                       edgecolor=P_DARK if call else "white",
                                       lw=2.2 if call else 0.8))
                dark = v / a.vmax > 0.55
                tc = "white" if dark else INK
                ax.text(j + 0.5, i + 0.30, "q %s" % fmt_q(qt, a.alpha), ha="center",
                        va="center", fontsize=9, fontweight="bold" if call else "normal",
                        color=tc)
                ax.text(j + 0.5, i + 0.55, "emp p %s" % fmt_emp(emp, r.get("emp_count")),
                        ha="center",
                        va="center", fontsize=7, color=tc)
                ax.text(j + 0.5, i + 0.78, "k %s/%s  (q_all %s)" % (kp, kn, fmt_q(qa, a.alpha)),
                        ha="center", va="center", fontsize=6.5, color=tc)
                src.append([locus, reg, c, nb, "yes", "%.6g" % emp, "%.6g" % qt,
                            "" if qa is None else "%.6g" % qa, kp, kn])
        ax.set_xticks([j + 0.5 for j in range(nc)])
        ax.set_xticklabels([cmp_label(c) for c in cmps], fontsize=8.5)
        ax.xaxis.tick_top()
        ax.set_yticks([i + 0.5 for i in range(nr)])
        ax.set_yticklabels([REGION_LABELS.get(r, r) for r in regions], fontsize=8.5)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(LOCUS_LABELS.get(locus, locus), loc="left", fontsize=11,
                     fontweight="bold", color=P_DARK, pad=26)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, a.vmax))
    cax = fig.add_axes([0.25, 0.55 / fig_h, 0.5, 0.13 / fig_h])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("-log10 q_testable (contrast-matched empirical null; capped at %g)"
                 % a.vmax, fontsize=8, color=G_DARK)
    cb.ax.tick_params(labelsize=7, color=G_MID)
    cb.ax.xaxis.set_label_position("top")
    cb.outline.set_edgecolor(G_LIGHT)
    fig.text(0.02, 0.012,
             "Outlined: regional call (q_testable < %g). k = bins at p < 0.05, "
             "first condition higher / second higher. Hatched: < 5 tested bins. "
             "<=: 0 of the matched background draws as extreme (null floor)." % a.alpha,
             fontsize=7, color=G_MID)
    fig.subplots_adjust(left=0.22, right=0.98, top=1 - 0.75 / fig_h, bottom=1.15 / fig_h)

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    fig.savefig(a.out_prefix + ".pdf")
    fig.savefig(a.out_prefix + ".png", dpi=300)
    with open(a.out_prefix + ".source.tsv", "w") as fh:
        fh.write("\t".join(["locus", "region", "comparison", "tested_bins", "testable",
                            "emp_p", "q_testable", "q_all", "k_pos", "k_neg"]) + "\n")
        for r in src:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("wrote %s.{pdf,png,source.tsv} (%d cells)" % (a.out_prefix, len(src)))


if __name__ == "__main__":
    main()
