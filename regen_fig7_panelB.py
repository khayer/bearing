#!/usr/bin/env python3
# regen_fig7_panelB.py  (v3 - matches bearing_pvalue.py exactly)
# Figure 7 Panel B: pooled empirical permutation null for the DN-vs-EbKO
# differential, N=100 production run, with the rejected Gamma fit overlaid.
#
# Per-bin statistic MATCHES bearing_pvalue.py diff_mode:
#   score_total = sum of the 6 per-track KL scores (signed)
#   tested      = abs(score_total)                     <-- the null statistic
#   bins kept if abs(score_total) >= min_signal (0.1)  <-- same filter as the test
# (default --stat abs_sum reproduces this; sum_abs/sum/max_abs are diagnostics.)
#
#   conda activate bearing
#   python regen_fig7_panelB.py --comp DN_vs_EbKO --out fig7_panelB.pdf
#   # fast preview:  --sample 10
# ASCII only. Reads real data; fabricates nothing.

import argparse, glob, gzip, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCORE_RE = re.compile(r"\[\s*(-?[0-9.eE+]+)\s*,\s*\d+\s*\]")

def bin_stats_from_file(path, stat, min_signal, sample):
    # per-bin statistic for (sampled) bins, filtered like the test (abs(sum)>=min_signal)
    out = []
    with gzip.open(path, "rt") as fh:
        for ln, line in enumerate(fh):
            if sample > 1 and (ln % sample):
                continue
            if line.startswith("#") or not line.strip():
                continue
            tab = line.rstrip("\n").split("\t")
            if len(tab) < 4:
                continue
            sc = SCORE_RE.findall(tab[3])
            if not sc:
                continue
            v = np.fromiter((float(s) for s in sc), dtype=float)
            signed_sum = v.sum()
            if abs(signed_sum) < min_signal:      # same filter as bearing_pvalue diff_mode
                continue
            if stat == "abs_sum":
                out.append(abs(signed_sum))        # == tested in bearing_pvalue.py
            elif stat == "sum":
                out.append(signed_sum)
            elif stat == "sum_abs":
                out.append(np.abs(v).sum())
            elif stat == "max_abs":
                out.append(np.abs(v).max())
            else:
                sys.exit("unknown --stat %s" % stat)
    return np.asarray(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="workflow/results")
    ap.add_argument("--comp", default="DN_vs_EbKO")
    ap.add_argument("--stat", default="abs_sum",
                    choices=["abs_sum", "sum", "sum_abs", "max_abs"])
    ap.add_argument("--min-signal", type=float, default=0.1,
                    help="match bearing_pvalue.py diff_mode filter (default 0.1)")
    ap.add_argument("--sample", type=int, default=1,
                    help="use every Nth line for speed (1 = all)")
    ap.add_argument("--nbins", type=int, default=60)
    ap.add_argument("--out", default="fig7_panelB.pdf")
    a = ap.parse_args()

    perm_files = sorted(glob.glob(os.path.join(
        a.results, "perm", "perm*", "diff_comparison",
        "diff_%s.qcat.bgz" % a.comp)))
    if not perm_files:
        sys.exit("no perm files under %s/perm/perm*/diff_comparison/" % a.results)
    print("pooling %d perms for %s (stat=%s, min_signal=%.3g, sample=1/%d)"
          % (len(perm_files), a.comp, a.stat, a.min_signal, a.sample))

    first = bin_stats_from_file(perm_files[0], a.stat, a.min_signal, a.sample)
    if first.size == 0:
        sys.exit("no bins pass the filter in the first file; check --stat/--min-signal")
    lo = max(first.min(), a.min_signal * 0.9, 1e-6)
    hi = first.max() * 1.05
    edges = np.logspace(np.log10(lo), np.log10(hi), a.nbins)

    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    moments = []
    for i, pf in enumerate(perm_files):
        v = first if i == 0 else bin_stats_from_file(pf, a.stat, a.min_signal, a.sample)
        counts += np.histogram(v, bins=edges)[0]
        if sum(m.size for m in moments) < 5_000_000:
            moments.append(v[::5])
        if (i + 1) % 20 == 0:
            print("  ...%d/%d files" % (i + 1, len(perm_files)))
    moments = np.concatenate(moments)
    print("  pooled null: n=%d  median=%.4g  max=%.4g"
          % (counts.sum(), np.median(moments), hi))

    # observed (unpermuted) diff -- lives under compare/
    obs = None
    for cand in ("compare/diff_%s.qcat.bgz" % a.comp,
                 "diff_comparison/diff_%s.qcat.bgz" % a.comp):
        p = os.path.join(a.results, cand)
        if os.path.exists(p):
            obs = bin_stats_from_file(p, a.stat, a.min_signal, a.sample)
            print("  observed from %s: n=%d median=%.4g" % (p, obs.size, np.median(obs)))
            break

    widths = np.diff(edges)
    dens = counts / (counts.sum() * widths)
    centers = np.sqrt(edges[:-1] * edges[1:])

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    ax.bar(centers, dens, width=widths, align="center", color="#4c72b0",
           alpha=0.75, label="pooled permutation null (N=%d)" % len(perm_files))
    try:
        from scipy.stats import gamma
        m, var = moments.mean(), moments.var()
        shape = m * m / var; scale = var / m
        xs = np.logspace(np.log10(lo), np.log10(hi), 300)
        ax.plot(xs, gamma.pdf(xs, shape, scale=scale), "r--", lw=1.5,
                label="Gamma fit (rejected)")
    except Exception as e:
        print("  (Gamma overlay skipped: %s)" % e)
    if obs is not None and obs.size:
        ax.axvline(np.median(obs), color="k", ls=":", lw=1.2, label="observed median")
    ax.set_xscale("log")
    ax.set_xlabel("per-bin |differential| (abs sum of 6 tracks)")
    ax.set_ylabel("density")
    ax.set_title("Empirical null vs Gamma: %s" % a.comp.replace("_vs_", " vs "))
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote %s" % a.out)

if __name__ == "__main__":
    main()
