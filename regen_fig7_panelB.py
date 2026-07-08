#!/usr/bin/env python3
# regen_fig7_panelB.py  (v5 - honest null panel)
# Figure 7 Panel B: pooled empirical permutation null for a differential,
# from the frozen N=100 production run.
#
# WHY NO SINGLE-GAMMA "REJECTED" OVERLAY:
#   BEARING's parametric option is a 2-component Gamma MIXTURE (fit_gamma_mixture
#   in bearing_pvalue.py) that must separate a low-score background/null component
#   from a high-score signal component. It fails on these unimodal replicate-
#   differential distributions (the mixture collapses; pi_bg degenerates), which
#   is why the empirical permutation null is used. A single Gamma actually fits
#   the pooled distribution well (KS D~0.03), so a single-Gamma "rejected" panel
#   is self-contradictory -- do not use it.
#
# Per-bin statistic MATCHES bearing_pvalue.py diff_mode: tested = abs(sum of the
# 6 per-track KL scores), bins kept if abs(sum) >= min_signal (0.1).
#
#   conda activate bearing
#   python regen_fig7_panelB.py --comp DN_vs_EbKO --out fig7_panelB.pdf
#   # demonstrate the parametric failure (runs the pipeline's own mixture EM):
#   python regen_fig7_panelB.py --comp DN_vs_EbKO --show-mixture --out fig7_panelB.pdf
# ASCII only. Reads real data; fabricates nothing.

import argparse, glob, gzip, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCORE_RE = re.compile(r"\[\s*(-?[0-9.eE+]+)\s*,\s*\d+\s*\]")

def bin_stats(path, min_signal, sample):
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
            s = sum(float(x) for x in sc)
            if abs(s) >= min_signal:
                out.append(abs(s))
    return np.asarray(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="workflow/results")
    ap.add_argument("--comp", default="DN_vs_EbKO")
    ap.add_argument("--min-signal", type=float, default=0.1)
    ap.add_argument("--sample", type=int, default=1)
    ap.add_argument("--nbins", type=int, default=60)
    ap.add_argument("--show-mixture", action="store_true",
                    help="run bearing_pvalue.fit_gamma_mixture on the OBSERVED "
                         "scores and overlay it, to show the parametric failure")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default="fig7_panelB.pdf")
    a = ap.parse_args()

    perm_files = sorted(glob.glob(os.path.join(
        a.results, "perm", "perm*", "diff_comparison",
        "diff_%s.qcat.bgz" % a.comp)))
    if not perm_files:
        sys.exit("no perm files under %s/perm/perm*/diff_comparison/" % a.results)
    print("pooling %d perms for %s (min_signal=%.3g, sample=1/%d)"
          % (len(perm_files), a.comp, a.min_signal, a.sample))

    first = bin_stats(perm_files[0], a.min_signal, a.sample)
    if first.size == 0:
        sys.exit("no bins pass the filter; check inputs/--min-signal")
    lo = max(first.min(), a.min_signal * 0.9, 1e-6)
    hi = first.max() * 1.05
    edges = np.logspace(np.log10(lo), np.log10(hi), a.nbins)
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    for i, pf in enumerate(perm_files):
        v = first if i == 0 else bin_stats(pf, a.min_signal, a.sample)
        counts += np.histogram(v, bins=edges)[0]
        if (i + 1) % 20 == 0:
            print("  ...%d/%d" % (i + 1, len(perm_files)))
    print("  pooled null: n=%d" % counts.sum())

    # observed
    obs = None
    for cand in ("compare/diff_%s.qcat.bgz" % a.comp,
                 "diff_comparison/diff_%s.qcat.bgz" % a.comp):
        p = os.path.join(a.results, cand)
        if os.path.exists(p):
            obs = bin_stats(p, a.min_signal, a.sample)
            print("  observed from %s: n=%d median=%.4g" % (p, obs.size, np.median(obs)))
            break

    widths = np.diff(edges)
    dens = counts / (counts.sum() * widths)
    centers = np.sqrt(edges[:-1] * edges[1:])

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    ax.bar(centers, dens, width=widths, align="center", color="#4c72b0", alpha=0.75)
    handles = [Patch(facecolor="#4c72b0", alpha=0.75,
                     label="pooled permutation null (N=%d)" % len(perm_files))]
    note = ""

    if a.show_mixture:
        if obs is None:
            print("  --show-mixture needs the observed diff; not found, skipping")
        else:
            try:
                sys.path.insert(0, os.getcwd())
                from bearing_pvalue import fit_gamma_mixture
                samp = obs if obs.size <= 2_000_000 else \
                       np.random.default_rng(42).choice(obs, 2_000_000, replace=False)
                # capture the pi_bg the pipeline reports
                import io, contextlib
                buf = io.StringIO()
                with contextlib.redirect_stderr(buf):
                    fitted_bg, k_bg, th_bg, _ = fit_gamma_mixture(samp)
                log = buf.getvalue()
                pi = re.search(r"pi=([0-9.]+)", log)
                pib = float(pi.group(1)) if pi else float("nan")
                conv = "converged" in log and "reached max_iter" not in log
                from scipy.stats import gamma as gd
                xs = np.logspace(np.log10(lo), np.log10(hi), 300)
                ax.plot(xs, gd.pdf(xs, a=k_bg, scale=th_bg), "r--", lw=1.5)
                handles.append(Line2D([0], [0], color="r", ls="--", lw=1.5,
                                      label="mixture-EM background (parametric null)"))
                bad = (not conv) or (pib < 0.3 or pib > 0.97) or np.isnan(pib)
                note = ("parametric mixture-EM FAILED\n"
                        "(pi_bg=%.2f, %s): no separable\nbackground -> empirical null used"
                        % (pib, "did not converge" if not conv else "collapsed")) if bad \
                       else "mixture-EM pi_bg=%.2f" % pib
                print("  mixture-EM: pi_bg=%.3f converged=%s -> %s"
                      % (pib, conv, "FAILED" if bad else "ok"))
            except Exception as e:
                print("  --show-mixture failed: %s" % e)

    if obs is not None and obs.size:
        ax.axvline(np.median(obs), color="k", ls=":", lw=1.2)
        handles.append(Line2D([0], [0], color="k", ls=":", lw=1.2, label="observed median"))

    ax.set_xscale("log")
    ax.set_xlabel("per-bin |differential| (abs sum of 6 tracks)")
    ax.set_ylabel("density")
    ax.set_title(a.title if a.title else
                 "Pooled empirical null: %s" % a.comp.replace("_vs_", " vs "))
    leg = ax.legend(handles=handles, fontsize=8, frameon=False, loc="upper right",
                    title="null: |differential| >= %.2g (min_signal)" % a.min_signal)
    leg.get_title().set_fontsize(7)
    if note:
        ax.text(0.03, 0.03, note, transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote %s" % a.out)

if __name__ == "__main__":
    main()
