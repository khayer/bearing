#!/usr/bin/env python3
# regen_fig7_panelB.py
# Regenerate Figure 7 Panel B: pooled empirical permutation null for the
# DN-vs-EbKO differential, from the frozen N=100 production run.
# Panel B shows why the parametric Gamma null was rejected: the pooled
# circular-shift null is plotted as a histogram with the (failed) Gamma fit
# overlaid for contrast.
#
# Run from the repo root (paths are relative to workflow/results).
# ASCII only. Reads real data; fabricates nothing.
#
#   conda activate bearing
#   python regen_fig7_panelB.py --comp DN_vs_EbKO --out fig7_panelB.pdf

import argparse, glob, gzip, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_scores_from_qcat(path, score_col=None):
    # perm diff qcats are bgzipped TSVs. Return the per-bin differential
    # statistic column as a float array. Adjust score_col if the header
    # name differs in your build (printed on first file for confirmation).
    vals = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if score_col is None:
            # pick the differential-score column by common names
            for cand in ("diff_kl", "dkl", "score", "stat", "kl", "d_stat"):
                if cand in header:
                    score_col = cand; break
        if score_col is None or score_col not in header:
            sys.exit("could not find score column in %s; header=%s "
                     "-- pass --score-col" % (path, header))
        ci = header.index(score_col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try: vals.append(float(f[ci]))
            except (ValueError, IndexError): pass
    return np.asarray(vals), score_col

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="workflow/results")
    ap.add_argument("--comp", default="DN_vs_EbKO")
    ap.add_argument("--score-col", default=None)
    ap.add_argument("--out", default="fig7_panelB.pdf")
    a = ap.parse_args()

    perm_glob = os.path.join(a.results, "perm", "perm*",
                             "diff_comparison", "diff_%s.qcat.bgz" % a.comp)
    perm_files = sorted(glob.glob(perm_glob))
    if not perm_files:
        sys.exit("no perm files matched %s" % perm_glob)
    print("pooling %d permutation files for %s" % (len(perm_files), a.comp))

    pooled = []
    score_col = a.score_col
    for i, pf in enumerate(perm_files):
        v, score_col = read_scores_from_qcat(pf, score_col)
        pooled.append(v)
        if i == 0:
            print("  using score column: %s" % score_col)
    pooled = np.concatenate(pooled)
    pooled = pooled[np.isfinite(pooled)]
    print("  pooled null: n=%d  min=%.4g  median=%.4g  max=%.4g"
          % (pooled.size, pooled.min(), np.median(pooled), pooled.max()))

    # observed per-bin differential statistic (same column) for the real data
    obs_path = os.path.join(a.results, "pvalue", "diff_%s.stats.tsv" % a.comp)
    obs = None
    if os.path.exists(obs_path):
        with open(obs_path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            if score_col in head:
                ci = head.index(score_col)
                obs = np.array([float(x.split("\t")[ci])
                                for x in fh
                                if _num(x.split("\t"), ci)])
        if obs is not None:
            obs = obs[np.isfinite(obs)]
            print("  observed: n=%d  median=%.4g  max=%.4g"
                  % (obs.size, np.median(obs), obs.max()))

    # plot: pooled empirical null histogram + failed Gamma fit overlay
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    lo = max(pooled[pooled > 0].min(), 1e-6)
    bins = np.logspace(np.log10(lo), np.log10(pooled.max()), 60)
    ax.hist(pooled, bins=bins, density=True, color="#4c72b0",
            alpha=0.75, label="pooled permutation null (N=%d)" % len(perm_files))
    # failed Gamma fit (method-of-moments) -- shown to illustrate poor fit
    try:
        from scipy.stats import gamma
        p = pooled[pooled > 0]
        m, v = p.mean(), p.var()
        shape = m * m / v; scale = v / m
        xs = np.logspace(np.log10(lo), np.log10(pooled.max()), 200)
        ax.plot(xs, gamma.pdf(xs, shape, scale=scale), "r--", lw=1.5,
                label="Gamma fit (rejected)")
    except Exception as e:
        print("  (scipy Gamma overlay skipped: %s)" % e)
    if obs is not None and obs.size:
        ax.axvline(np.median(obs), color="k", ls=":", lw=1.2,
                   label="observed median")
    ax.set_xscale("log")
    ax.set_xlabel("per-bin differential statistic (%s)" % score_col)
    ax.set_ylabel("density")
    ax.set_title("Figure 7B: pooled empirical null, %s (N=%d perms)"
                 % (a.comp, len(perm_files)))
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote %s" % a.out)

def _num(fields, ci):
    try: float(fields[ci]); return True
    except Exception: return False

if __name__ == "__main__":
    main()
