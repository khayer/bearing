#!/usr/bin/env python3
# regen_fig7_panelB.py  (v2 - correct qcat parser)
# Figure 7 Panel B: pooled empirical permutation null for the DN-vs-EbKO
# differential, from the frozen N=100 production run, with the rejected
# Gamma fit overlaid.
#
# qcat.bgz format (headerless, tab-separated):
#   chrom  start  end  id:<n>,qcat:[[score,track],[score,track],... x6]
# The 6 scores are the per-track differential (clamped-KL) contributions at
# that bin. The per-bin statistic below defaults to sum of |score| across the
# 6 tracks (L1 differential magnitude). If bearing_pvalue.py thresholds a
# different per-bin statistic, set --stat to match (sum_abs | sum | max_abs).
#
#   conda activate bearing
#   python regen_fig7_panelB.py --comp DN_vs_EbKO --out fig7_panelB.pdf
#   # faster preview on 1/10 of bins:  --sample 10
# ASCII only. Reads real data; fabricates nothing.

import argparse, glob, gzip, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCORE_RE = re.compile(r"\[\s*(-?[0-9.eE+]+)\s*,\s*\d+\s*\]")

def bin_stats_from_file(path, stat, sample):
    # yield the per-bin statistic for every (sampled) bin in one qcat file
    out = []
    with gzip.open(path, "rt") as fh:
        for ln, line in enumerate(fh):
            if sample > 1 and (ln % sample): 
                continue
            tab = line.rstrip("\n").split("\t")
            if len(tab) < 4:
                continue
            scores = SCORE_RE.findall(tab[3])
            if not scores:
                continue
            v = np.fromiter((float(s) for s in scores), dtype=float)
            if stat == "sum_abs":
                out.append(np.abs(v).sum())
            elif stat == "sum":
                out.append(v.sum())
            elif stat == "max_abs":
                out.append(np.abs(v).max())
            else:
                sys.exit("unknown --stat %s" % stat)
    return np.asarray(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="workflow/results")
    ap.add_argument("--comp", default="DN_vs_EbKO")
    ap.add_argument("--stat", default="sum_abs",
                    choices=["sum_abs", "sum", "max_abs"])
    ap.add_argument("--sample", type=int, default=1,
                    help="use every Nth bin for speed (1 = all bins)")
    ap.add_argument("--nbins", type=int, default=60)
    ap.add_argument("--out", default="fig7_panelB.pdf")
    a = ap.parse_args()

    perm_files = sorted(glob.glob(os.path.join(
        a.results, "perm", "perm*", "diff_comparison",
        "diff_%s.qcat.bgz" % a.comp)))
    if not perm_files:
        sys.exit("no perm files under %s/perm/perm*/diff_comparison/" % a.results)
    print("pooling %d permutation files for %s (stat=%s, sample=1/%d)"
          % (len(perm_files), a.comp, a.stat, a.sample))

    # pass 1: first file sets the histogram range from its nonzero values
    first = bin_stats_from_file(perm_files[0], a.stat, a.sample)
    nz = first[first > 0]
    if nz.size == 0:
        sys.exit("all-zero first file; check --stat / inputs")
    lo, hi = np.percentile(nz, 0.1), first.max()
    lo = max(lo, 1e-6)
    edges = np.logspace(np.log10(lo), np.log10(hi * 1.05), a.nbins)
    print("  range: [%.4g, %.4g]  (first-file median nonzero=%.4g)"
          % (lo, hi, np.median(nz)))

    # accumulate counts across all files (memory-safe)
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    total = 0; total_nz = 0
    pooled_for_moments = []   # keep a subsample for the Gamma fit only
    for i, pf in enumerate(perm_files):
        v = first if i == 0 else bin_stats_from_file(pf, a.stat, a.sample)
        c, _ = np.histogram(v[v > 0], bins=edges)
        counts += c
        total += v.size; total_nz += int((v > 0).sum())
        if len(pooled_for_moments) < 5_000_000:
            pooled_for_moments.append(v[v > 0][::5])
        if (i + 1) % 20 == 0:
            print("  ...%d/%d files" % (i + 1, len(perm_files)))
    pooled_for_moments = np.concatenate(pooled_for_moments)
    print("  pooled: %d bins (%d nonzero)  median=%.4g"
          % (total, total_nz, np.median(pooled_for_moments)))

    # observed (unpermuted) diff, if present
    obs = None
    for cand in ("diff_comparison/diff_%s.qcat.bgz" % a.comp,
                 "diff_%s.qcat.bgz" % a.comp):
        p = os.path.join(a.results, cand)
        if os.path.exists(p):
            ov = bin_stats_from_file(p, a.stat, a.sample)
            obs = ov[ov > 0]
            print("  observed from %s: median=%.4g" % (p, np.median(obs)))
            break

    # density from accumulated counts
    widths = np.diff(edges)
    dens = counts / (counts.sum() * widths)
    centers = np.sqrt(edges[:-1] * edges[1:])

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    ax.bar(centers, dens, width=widths, align="center", color="#4c72b0",
           alpha=0.75, label="pooled permutation null (N=%d)" % len(perm_files))
    try:
        from scipy.stats import gamma
        m, var = pooled_for_moments.mean(), pooled_for_moments.var()
        shape = m * m / var; scale = var / m
        xs = np.logspace(np.log10(lo), np.log10(hi), 300)
        ax.plot(xs, gamma.pdf(xs, shape, scale=scale), "r--", lw=1.5,
                label="Gamma fit (rejected)")
    except Exception as e:
        print("  (Gamma overlay skipped: %s)" % e)
    if obs is not None and obs.size:
        ax.axvline(np.median(obs), color="k", ls=":", lw=1.2,
                   label="observed median")
    ax.set_xscale("log")
    ax.set_xlabel("per-bin differential statistic (%s of 6 tracks)" % a.stat)
    ax.set_ylabel("density")
    ax.set_title("Figure 7B: pooled empirical null, %s (N=%d perms)"
                 % (a.comp, len(perm_files)))
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote %s" % a.out)

if __name__ == "__main__":
    main()
