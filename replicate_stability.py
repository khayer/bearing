#!/usr/bin/env python3
"""
replicate_stability.py
======================
Leave-one-replicate-out stability of the BEARING differential.

WHY
The production p-values come from the permutation null; replicates are averaged
at the score level before the differential is taken. A reviewer asked whether
that averaging hides replicate disagreement -- i.e. whether one replicate drives
a call. With two biological replicates per condition, dropping one leaves a
single-replicate differential; the four ways of pairing one A-replicate with one
B-replicate are the complete leave-one-out set. If the biology is real, all four
pairings should (i) recover the same regions and (ii) attribute the change to the
same track, and the per-bin differential sign should be concordant across them.

WHAT (three readouts, per key contrast)
1. Regional reproducibility: top-1% fold-enrichment + hypergeometric p for each
   pre-specified region, in each of the four single-replicate pairings.
2. Per-track attribution stability: the per-track share of the absolute
   differential within each region, in each pairing, with the dominant track and
   the coefficient of variation of its share across pairings.
3. Per-bin sign concordance: fraction of bins on which all four pairings agree on
   the sign of the differential, genome-wide and within each region.

HOW (reuses verified code)
Imports load_matrix, read_sheet, read_chrom_sizes, compute_Q_all_bins,
load_regions, region_bin_mask, hyper_sf and import_production_scorer from
baseline_comparison.py (whose per-sample reconstruction is validated against the
production qcat, Spearman 0.9987). Each of the up-to-four replicate samples is
loaded once; the per-track KL score is computed per sample against that sample's
own genome-wide Q, exactly as production scores each sample. The single-replicate
differential per track is klA[:,i] - klB[:,i], which sums to bearing_kl(A) -
bearing_kl(B); the tested differential is the absolute summed value, matching the
production diff_mode statistic. Because each per-sample per-track KL is faithful
to production, so is its per-track difference -- no averaged-differential
reproduction is required.

USAGE
  conda activate bearing
  python replicate_stability.py \\
      --sheet workflow/config/samples.tsv \\
      --chrom-sizes workflow/resources/mm10.chrom.sizes \\
      --regions regions_manuscript.tsv \\
      --cond-a "DN:DN_rep1,DN_rep2" \\
      --cond-b "ProB:ProB_rep1,ProB_rep2" \\
      --track-names ATAC,RNA+,RNA-,CTCF,Cohesin,H3K27ac \\
      --out replicate_stability_DN_vs_ProB.tsv

  # then DN-vs-S3T3 (tcrb) and DN-vs-DP similarly.

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import os
import sys

import numpy as np


def load_baseline(repo_dir):
    """baseline_comparison.py lives in the repo root (where it is run from)."""
    sys.path.insert(0, repo_dir)
    try:
        import baseline_comparison as bc
    except Exception as e:
        sys.exit("could not import baseline_comparison.py from %s: %s\n"
                 "(run from the repo root, or pass --repo pointing at it)" % (repo_dir, e))
    return bc


def per_track_kl(R, kl_fn, prob_fn, Q):
    """bins x tracks production per-track KL (unsummed). Same call as
    baseline_comparison.statistics uses, minus the .sum(axis=1)."""
    P = prob_fn(R.astype(np.float64))
    res = kl_fn(P, Q, raw_signal_matrix=R, min_signal=0.0, score_method="kl")
    arr = res[0] if isinstance(res, tuple) else res
    return np.asarray(arr, dtype=np.float64)


def parse_cond(spec):
    """'DN:DN_rep1,DN_rep2' -> ('DN', ['DN_rep1','DN_rep2'])"""
    label, samples = spec.split(":")
    return label, [s.strip() for s in samples.split(",") if s.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", default="workflow/config/samples.tsv")
    ap.add_argument("--chrom-sizes", default="workflow/resources/mm10.chrom.sizes")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--cond-a", required=True, help="LABEL:sample1,sample2")
    ap.add_argument("--cond-b", required=True, help="LABEL:sample1,sample2")
    ap.add_argument("--track-names", default="ATAC,RNA+,RNA-,CTCF,Cohesin,H3K27ac")
    ap.add_argument("--data-dir", default=".")
    ap.add_argument("--repo", default=".")
    ap.add_argument("--min-signal", type=float, default=0.1)
    ap.add_argument("--top-frac", type=float, default=0.01)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    bc = load_baseline(os.path.abspath(a.repo))
    kl_fn, prob_fn = bc.import_production_scorer(os.path.abspath(a.repo))
    cs = bc.read_chrom_sizes(a.chrom_sizes, [])
    tnames = [t.strip() for t in a.track_names.split(",")]

    la, sa = parse_cond(a.cond_a)
    lb, sb = parse_cond(a.cond_b)
    all_samples = sa + sb
    print("condition A = %s %s ; condition B = %s %s" % (la, sa, lb, sb))

    # ---- load each sample once: per-track KL (genome-wide) + keep mask --------
    ktrack = {}   # sample -> bins x tracks KL
    keep = {}     # sample -> bool mask
    index = None
    for name in all_samples:
        print("\nloading %s" % name)
        bw = bc.read_sheet(a.sheet, name)
        R, idx = bc.load_matrix(bw, cs, a.data_dir)
        if index is None:
            index = idx
        Q = bc.compute_Q_all_bins(R, prob_fn)
        print("  Q =", np.round(Q, 4))
        ktrack[name] = per_track_kl(R, kl_fn, prob_fn, Q)
        keep[name] = R.sum(axis=1) >= a.min_signal
        del R

    # Common bins: pass the min_signal mask in ALL four samples. baseline_comparison
    # uses the 2-sample both-mask per pairing; we intersect across all four so the
    # pairings share one bin set (required for per-bin sign concordance and for a
    # like-for-like fold-enrichment comparison). This is marginally more conservative
    # -- slightly fewer bins -- and does not bias any single pairing.
    common = np.ones(ktrack[all_samples[0]].shape[0], dtype=bool)
    for name in all_samples:
        common &= keep[name]
    n_common = int(common.sum())
    print("\nbins passing in all %d samples: %d" % (len(all_samples), n_common))
    for name in all_samples:
        ktrack[name] = ktrack[name][common]

    # region masks on the common-bin coordinate system
    regions = bc.load_regions(a.regions)
    rmasks = [(lab, bc.region_bin_mask(index, c, s, e, common))
              for _, c, s, e, lab in regions]

    # ---- the four single-replicate pairings ----------------------------------
    pairings = [(ai, bj) for ai in sa for bj in sb]
    print("\n%d leave-one-out pairings: %s" % (len(pairings), pairings))

    diff_sum = {}       # pairing -> summed differential (signed) per bin
    diff_track = {}     # pairing -> per-track differential (signed) bins x tracks
    for (ai, bj) in pairings:
        dt = ktrack[ai] - ktrack[bj]           # per-track differential, signed
        diff_track[(ai, bj)] = dt
        diff_sum[(ai, bj)] = dt.sum(axis=1)    # = bearing_kl(A) - bearing_kl(B)

    N = n_common
    n_top = max(1, int(round(a.top_frac * N)))

    rows = []

    # ---- readout 1: regional reproducibility ---------------------------------
    print("\n=== 1. REGIONAL REPRODUCIBILITY (top-%.0f%% fold-enrichment; p = hypergeometric) ==="
          % (100 * a.top_frac))
    for lab, rmask in rmasks:
        K = int(rmask.sum())
        if K < 20:
            continue
        print("\n  %-26s (%d region bins / %d genome-wide)" % (lab, K, N))
        print("    %-22s %8s %10s %12s" % ("pairing", "hits", "fold", "p"))
        folds = []
        for (ai, bj) in pairings:
            absd = np.abs(diff_sum[(ai, bj)])
            top = np.argpartition(-absd, n_top - 1)[:n_top]
            sel = np.zeros(N, dtype=bool); sel[top] = True
            k = int((sel & rmask).sum())
            exp = n_top * K / float(N)
            fold = k / exp if exp else float("nan")
            p = bc.hyper_sf(k, N, K, n_top) if k > 0 else 1.0
            folds.append(fold)
            print("    %-22s %8d %10.2f %12.3g" % ("%s vs %s" % (ai, bj), k, fold, p))
            rows.append(("reproducibility", lab, "%s_vs_%s" % (ai, bj), "fold", fold, p))
        cv = (np.std(folds) / np.mean(folds)) if np.mean(folds) else float("nan")
        allsig = all(f > 1.5 for f in folds)
        print("    -> fold CV across pairings = %.2f ; enriched (>1.5x) in all four: %s"
              % (cv, "YES" if allsig else "no"))

    # ---- readout 2: per-track attribution stability --------------------------
    print("\n=== 2. PER-TRACK ATTRIBUTION STABILITY (share of |differential| within region) ===")
    for lab, rmask in rmasks:
        K = int(rmask.sum())
        if K < 20:
            continue
        print("\n  %-26s" % lab)
        print("    %-22s %s" % ("pairing", "  ".join("%-8s" % t for t in tnames)))
        dom = []
        topshares = []
        for (ai, bj) in pairings:
            contrib = np.abs(diff_track[(ai, bj)][rmask]).sum(axis=0)   # per track
            tot = contrib.sum()
            share = contrib / tot if tot else np.zeros_like(contrib)
            di = int(np.argmax(share)); dom.append(di); topshares.append(share[di])
            print("    %-22s %s" % ("%s vs %s" % (ai, bj),
                  "  ".join("%7.1f%%" % (100 * x) for x in share)))
            for ti, x in enumerate(share):
                rows.append(("attribution", lab, "%s_vs_%s" % (ai, bj), tnames[ti], x, ""))
        same = len(set(dom)) == 1
        cvtop = (np.std(topshares) / np.mean(topshares)) if np.mean(topshares) else float("nan")
        print("    -> dominant track: %s ; same in all four: %s ; top-share CV = %.2f"
              % ("/".join(sorted(set(tnames[d] for d in dom))),
                 "YES" if same else "NO", cvtop))

    # ---- readout 3: per-bin sign concordance ---------------------------------
    print("\n=== 3. PER-BIN SIGN CONCORDANCE across the four pairings ===")
    signs = np.vstack([np.sign(diff_sum[p]) for p in pairings])   # 4 x bins
    # a bin is concordant if all four non-zero signs agree
    nz = np.all(signs != 0, axis=0)
    allpos = np.all(signs > 0, axis=0)
    allneg = np.all(signs < 0, axis=0)
    concord = allpos | allneg
    gw = float(concord[nz].mean()) if nz.sum() else float("nan")
    print("  genome-wide (bins non-zero in all four): %.3f concordant" % gw)
    rows.append(("sign_concordance", "genome_wide", "all_pairings", "frac", gw, ""))
    print("  %-26s %10s %12s" % ("region", "n_nz", "concordant"))
    for lab, rmask in rmasks:
        m = rmask & nz
        if m.sum() < 10:
            print("  %-26s %10d   (too few)" % (lab, int(m.sum()))); continue
        c = float(concord[m].mean())
        print("  %-26s %10d %11.3f" % (lab, int(m.sum()), c))
        rows.append(("sign_concordance", lab, "all_pairings", "frac", c, ""))

    with open(a.out, "w") as fh:
        fh.write("readout\tregion\tpairing\tmetric\tvalue\tp\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%s\t%.6g\t%s\n" % r)
    print("\nwrote %s" % a.out)
    print("\nInterpretation: a real call is enriched in all four pairings (readout 1),")
    print("attributes to the same dominant track in all four (readout 2), and has high")
    print("per-bin sign concordance in-region (readout 3). Disagreement would indicate")
    print("the averaged call depends on which replicate is used.")


if __name__ == "__main__":
    main()
