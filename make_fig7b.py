#!/usr/bin/env python3
# make_fig7b.py -- Figure 7B: pooled empirical permutation null for one differential comparison
# Authors: Katharina E. Hayer (CHOP/Drexel) with Claude (Anthropic), co-created
# Version: 1.1   Date: 2026-09-17  (v1.1: caches the histogram in <out-prefix>.npz; --from-npz re-renders without streaming; bounded log axis)
#
# Streams the N permuted diff qcat files (never holds them all in memory), builds a fixed-bin
# histogram of |diff score| above the production min-signal floor, overlays the 2-component
# Gamma-mixture fit (the parametric alternative), marks the observed FDR-threshold score, and
# writes the legend numbers (total null bins, min achievable p, Gamma false-call fraction) to JSON.
#
# Usage (same inputs as the production bearing_pvalue.py --diff call):
#   python make_fig7b.py --null-qcat perm_nulls/perm*/diff/diff_DN_vs_EbKO.qcat.bgz \
#       --observed-stats workflow/results/pvalue/diff_DN_vs_EbKO.stats.tsv \
#       --label "DN - EbKO" --min-signal 0.5 --out-prefix Fig7B_DN_vs_EbKO
#
import argparse, json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bearing_pvalue import parse_qcat, fit_gamma_mixture  # same repo, same parser

ap = argparse.ArgumentParser()
ap.add_argument("--null-qcat", nargs="+", required=True)
ap.add_argument("--observed-stats", required=True, help="<prefix>.stats.tsv from bearing_pvalue.py --diff")
ap.add_argument("--label", default="DN - EbKO")
ap.add_argument("--min-signal", type=float, default=0.5, help="production floor on |score| (KL: 0.5)")
ap.add_argument("--xmax", type=float, default=6.0)
ap.add_argument("--nbins", type=int, default=220)
ap.add_argument("--gamma-subsample", type=int, default=5_000_000)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--out-prefix", required=True)
ap.add_argument("--from-npz", default=None, help="re-render from a cached <prefix>.npz instead of streaming the null files")
a = ap.parse_args()

edges = np.linspace(a.min_signal, a.xmax, a.nbins + 1)
counts = np.zeros(a.nbins, dtype=np.int64)
total = 0; overflow = 0; per_file = []
if a.from_npz:
    z = np.load(a.from_npz, allow_pickle=True)
    edges, counts, total, overflow = z["edges"], z["counts"], int(z["total"]), int(z["overflow"])
    reservoir = list(z["reservoir"]); per_file = [tuple(t) for t in z["per_file"]]
    a.null_qcat = [f for f, _ in per_file]
rng = np.random.default_rng(a.seed)
if not a.from_npz:
    reservoir = []  # uniform subsample of null scores for the Gamma fit and the false-call fraction
seen = 0
for path in ([] if a.from_npz else a.null_qcat):
    n = 0
    buf = []
    for _, _, _, score, _ in parse_qcat(path, min_signal=a.min_signal, diff_mode=True):
        s = abs(score)
        if s < a.min_signal:
            continue
        buf.append(s); n += 1
        if len(buf) >= 2_000_000:
            arr = np.asarray(buf); buf = []
            h, _ = np.histogram(arr, bins=edges); counts += h; overflow += int((arr > a.xmax).sum())
            for v in arr[rng.random(len(arr)) < (a.gamma_subsample / 7.5e8)]:
                reservoir.append(v)
    if buf:
        arr = np.asarray(buf)
        h, _ = np.histogram(arr, bins=edges); counts += h; overflow += int((arr > a.xmax).sum())
        for v in arr[rng.random(len(arr)) < (a.gamma_subsample / 7.5e8)]:
            reservoir.append(v)
    per_file.append((path, n)); total += n
    print(f"{path}: {n:,} null bins", file=sys.stderr)
np.savez(a.out_prefix + ".npz", edges=edges, counts=counts, total=total, overflow=overflow, reservoir=np.asarray(reservoir), per_file=np.asarray(per_file, dtype=object))
if any(int(n) == 0 for _, n in per_file):
    sys.exit("ERROR: a null file contributed zero bins; refusing to draw a partial null")
min_p = 1.0 / (total + 1)
sub = np.asarray(reservoir)
print(f"total null bins {total:,}; min achievable p {min_p:.2e}; Gamma subsample {len(sub):,}", file=sys.stderr)

# parametric alternative on the same null scores
fitted, k_bg, th_bg, cap = fit_gamma_mixture(sub)
gamma_p = 1.0 - fitted.cdf(sub)
false_call = float((gamma_p < 0.05).mean())

# observed: FDR threshold score and largest observed |score|
import csv
obs = []; fdr_scores = []
with open(a.observed_stats) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    scol = "bearing_score" if "bearing_score" in rd.fieldnames else [c for c in rd.fieldnames if "score" in c][0]
    sigcol = [c for c in rd.fieldnames if c.startswith("significant")]
    for row in rd:
        try:
            s = abs(float(row[scol]))
        except ValueError:
            continue
        if s >= a.min_signal:
            obs.append(s)
            if sigcol and row[sigcol[0]] in ("True", "1", "true"):
                fdr_scores.append(s)
obs = np.asarray(obs)
fdr_thr = float(min(fdr_scores)) if fdr_scores else None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
centers = 0.5 * (edges[1:] + edges[:-1]); width = edges[1] - edges[0]
dens = counts / (total * width)
fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=200)
ax.bar(centers, dens, width=width, color="#b9c8ea", edgecolor="none", label=f"pooled empirical null, N = {len(a.null_qcat)} permutations")
xs = np.linspace(a.min_signal, a.xmax, 400)
ax.plot(xs, fitted.pdf(xs) / max(1e-12, 1.0 - fitted.cdf(a.min_signal)), "--", color="#555555", lw=1.2,
        label=f"2-component Gamma mixture (calls {false_call*100:.0f}% of null at p < 0.05)")
if len(obs):
    ho, _ = np.histogram(obs, bins=edges); ax.step(edges[:-1], ho / (len(obs) * width), where="post", color="#b00020", lw=1.0, label=f"observed {a.label}")
if fdr_thr:
    ax.axvline(fdr_thr, color="#b00020", lw=1.0, ls=":"); ax.text(fdr_thr, 5, " FDR 0.05", color="#b00020", fontsize=6.5)
ax.set_yscale("log"); ax.set_ylim(0.3 / (total * width), 20); ax.set_xlabel("|differential BEARING score|", fontsize=8); ax.set_ylabel("density", fontsize=8)
ax.tick_params(labelsize=7); ax.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1)
ax.set_title(f"Pooled empirical null, {a.label}\n{total:,.0f} null bins from {len(a.null_qcat)} permutations; min p = {min_p:.1e}", fontsize=7.5)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(a.out_prefix + ".png", dpi=300); fig.savefig(a.out_prefix + ".pdf")
json.dump({"label": a.label, "n_perms": len(a.null_qcat), "total_null_bins": total, "min_achievable_p": min_p,
           "min_signal": a.min_signal, "overflow_above_xmax": overflow, "gamma_k_bg": k_bg, "gamma_theta_bg": th_bg,
           "gamma_false_call_fraction_at_0.05": false_call, "fdr_threshold_score": fdr_thr,
           "observed_bins_above_floor": int(len(obs)), "per_file_counts": per_file},
          open(a.out_prefix + ".json", "w"), indent=2)
print("wrote", a.out_prefix + ".png/.pdf/.json")
