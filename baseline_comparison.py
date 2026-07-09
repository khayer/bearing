#!/usr/bin/env python3
# baseline_comparison.py
#
# AI Reviewer question: "is BEARING simply rediscovering what a simpler
# statistic would show?"  This script computes, for the SAME bins and the SAME
# min_signal mask, the BEARING score and a set of simple baselines a typical
# analyst would try, then reports how strongly they agree (rank correlation)
# and how much their top-scoring bins overlap (Jaccard of the top 1%).
#
# Statistics computed per bin (each a scalar summary of the 6-track vector):
#   bearing_kl     clamped per-track KL, sum over tracks   <- production score
#                  score_i = max(P_i * log2(P_i/Q_i), 0)
#   unclamped_kl   full KL(P||Q) = sum_i P_i * log2(P_i/Q_i)   (no clamp)
#   bearing_jsd    the pipeline's bounded JSD variant (score_method="jsd")
#   jsd_full       full Jensen-Shannon divergence JSD(P||Q), all tracks
#   neg_entropy    -H(P), i.e. composition "peakedness"
#   sum_z          sum_i z_i of raw signal (z from genome-wide mean/sd per track)
#   max_z          max_i |z_i|
#
# The clamped score is recomputed from the qcat's own `raw:` field and checked
# against the qcat's stored scores -- a built-in correctness assertion.
#
# Uses the PRODUCTION scorer (imported from bigwig_to_qcat), so bearing_kl and
# bearing_jsd are exactly what the pipeline computes, not a paraphrase.
#
#   conda activate bearing
#   # single sample, 1/20 of bins (fast):
#   python baseline_comparison.py --qcat workflow/results/DN_rep1.qcat.bgz --sample 20
#   # differential between two samples (the contrast the paper actually makes):
#   python baseline_comparison.py --qcat workflow/results/DN_rep1.qcat.bgz \
#          --vs workflow/results/ProB_rep1.qcat.bgz --sample 20
#
# ASCII only. Reads real data; fabricates nothing.

import argparse
import gzip
import os
import re
import sys

import numpy as np

RAW_RE = re.compile(r"raw:\[([^\]]*)\]")
QCAT_RE = re.compile(r"qcat:\[(.*?)\](?:,raw:|\s*$)", re.S)
SCORE_RE = re.compile(r"\[\s*(-?[0-9.eE+]+)\s*,\s*\d+\s*\]")


def import_production_scorer(repo_root):
    """Import the real scorer so baselines are compared against production code."""
    sys.path.insert(0, repo_root)
    try:
        from bigwig_to_qcat import kl_scores_per_bin, signals_to_prob, PSEUDOCOUNT
    except Exception as e:
        sys.exit("could not import bigwig_to_qcat from %s: %s\n"
                 "run from the repo root, or pass --repo" % (repo_root, e))
    return kl_scores_per_bin, signals_to_prob, PSEUDOCOUNT


def parse_qcat(path, sample):
    """Yield (raw_vector, stored_score_vector) for every (sampled) bin."""
    n_noraw = 0
    with gzip.open(path, "rt") as fh:
        for ln, line in enumerate(fh):
            if sample > 1 and (ln % sample):
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            col = parts[3]
            mraw = RAW_RE.search(col)
            if not mraw:
                n_noraw += 1
                continue
            raw = np.fromstring(mraw.group(1), sep=",")
            scores = np.array([float(s) for s in SCORE_RE.findall(col)])
            if raw.size == 0 or scores.size == 0:
                continue
            yield raw, scores
    if n_noraw:
        print("  WARNING: %d bins had no raw: field (skipped). If this is most of "
              "the file, this qcat predates the raw: field -- rescore or use "
              "bigWigs." % n_noraw, file=sys.stderr)


def collect(path, sample, min_signal):
    """One pass: accumulate per-track sums for Q and z-scores, keep the sampled rows."""
    raws, stored = [], []
    for raw, sc in parse_qcat(path, sample):
        raws.append(raw)
        stored.append(sc)
    if not raws:
        sys.exit("no usable bins parsed from %s (is the raw: field present?)" % path)
    R = np.vstack(raws)
    S = np.vstack(stored)
    # production mask: total RAW signal across tracks below min_signal -> zeroed
    keep = R.sum(axis=1) >= min_signal
    print("  %s: %d sampled bins, %d pass min_signal=%.3g (%.1f%% masked)"
          % (os.path.basename(path), R.shape[0], int(keep.sum()),
             100.0 * (1.0 - keep.mean())))
    return R[keep], S[keep]


def statistics(R, kl_scores_per_bin, signals_to_prob, Q=None):
    """Return dict of per-bin scalar statistics, plus the Q used."""
    P = signals_to_prob(R)
    if Q is None:
        Q = P.mean(axis=0)              # Q = mean of P across bins (definition)
    logPQ = np.log2(P / Q[None, :])

    out = {}
    # production scorer (clamped KL, and the pipeline's bounded jsd variant)
    kl = kl_scores_per_bin(P, Q, raw_signal_matrix=R, min_signal=0.0,
                           score_method="kl")
    out["bearing_kl"] = np.asarray(kl).sum(axis=1)
    try:
        jsd = kl_scores_per_bin(P, Q, raw_signal_matrix=R, min_signal=0.0,
                                score_method="jsd")
        out["bearing_jsd"] = np.asarray(jsd).sum(axis=1)
    except Exception as e:
        print("  (bearing_jsd unavailable: %s)" % e, file=sys.stderr)

    out["unclamped_kl"] = (P * logPQ).sum(axis=1)

    M = 0.5 * (P + Q[None, :])
    jsd_full = 0.5 * ((P * np.log2(P / M)).sum(axis=1)
                      + (Q[None, :] * np.log2(Q[None, :] / M)).sum(axis=1))
    out["jsd_full"] = jsd_full

    out["neg_entropy"] = (P * np.log2(P)).sum(axis=1)   # = -H(P)

    mu = R.mean(axis=0)
    sd = R.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (R - mu[None, :]) / sd[None, :]
    out["sum_z"] = Z.sum(axis=1)
    out["max_z"] = np.abs(Z).max(axis=1)
    return out, Q


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


def top_jaccard(a, b, frac=0.01):
    k = max(1, int(round(frac * a.size)))
    ta = set(np.argpartition(-a, k - 1)[:k].tolist())
    tb = set(np.argpartition(-b, k - 1)[:k].tolist())
    return len(ta & tb) / float(len(ta | tb))


def report(stats, ref="bearing_kl", frac=0.01):
    keys = [k for k in stats if k != ref]
    print("\n%-14s vs %s" % ("baseline", ref))
    print("%-14s %10s %14s" % ("", "Spearman", "top-%.0f%% Jaccard" % (100 * frac)))
    for k in keys:
        print("%-14s %+10.3f %14.3f"
              % (k, spearman(stats[ref], stats[k]), top_jaccard(stats[ref], stats[k], frac)))
    print("\nReading the table:")
    print("  A baseline with BOTH rho ~ 1 and Jaccard ~ 1 reproduces the BEARING")
    print("  ranking -> that baseline is sufficient and BEARING adds nothing on")
    print("  this axis. Low Jaccard with high rho means they agree globally but")
    print("  disagree on exactly the bins that get called -- which is the part")
    print("  that matters for the regional tests.")
    print("  NOTE: per-bin scalar agreement is NOT the whole story. BEARING's")
    print("  claim is per-track ATTRIBUTION; none of these baselines attribute.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qcat", required=True)
    ap.add_argument("--vs", default=None,
                    help="second sample: compare DIFFERENTIALS instead of raw scores")
    ap.add_argument("--repo", default=".", help="repo root (for bigwig_to_qcat import)")
    ap.add_argument("--min-signal", type=float, default=0.1)
    ap.add_argument("--sample", type=int, default=20, help="use every Nth bin")
    ap.add_argument("--top-frac", type=float, default=0.01)
    a = ap.parse_args()

    kl_fn, prob_fn, _ = import_production_scorer(os.path.abspath(a.repo))

    print("collecting %s" % a.qcat)
    R1, S1 = collect(a.qcat, a.sample, a.min_signal)
    st1, Q1 = statistics(R1, kl_fn, prob_fn)

    # correctness check: recomputed clamped score vs the score stored in the qcat
    stored = S1.sum(axis=1)
    rho = spearman(stored, st1["bearing_kl"])
    md = float(np.max(np.abs(stored - st1["bearing_kl"])))
    print("  self-check vs stored qcat scores: Spearman=%.4f  max|diff|=%.3g" % (rho, md))
    if rho < 0.99:
        print("  WARNING: recomputed clamped score does not track the stored score.")
        print("  Q here is the mean of P over SAMPLED bins; production Q is genome-wide.")
        print("  Re-run with --sample 1, or pass the cached Q, before trusting the table.")

    if a.vs is None:
        report(st1, frac=a.top_frac)
        return

    print("\ncollecting %s" % a.vs)
    R2, _ = collect(a.vs, a.sample, a.min_signal)
    n = min(R1.shape[0], R2.shape[0])
    if R1.shape[0] != R2.shape[0]:
        print("  NOTE: different bin counts after masking (%d vs %d); comparing the "
              "first %d. For a bin-exact differential use the compare/ diff qcats."
              % (R1.shape[0], R2.shape[0], n))
    st2, _ = statistics(R2[:n], kl_fn, prob_fn, Q=Q1)

    diff = {k: np.abs(st1[k][:n] - st2[k]) for k in st1 if k in st2}
    print("\nDIFFERENTIAL (|sample1 - sample2|) baselines:")
    report(diff, frac=a.top_frac)


if __name__ == "__main__":
    main()
