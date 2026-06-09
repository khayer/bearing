#!/usr/bin/env python3
"""Diagnose why a per-sample BEARING p-value track is flat (no FDR-significant
bins) when sibling samples are not.

Given several samples' per-sample stats tables (pvalue/<sample>.stats.tsv), this
compares, for each sample:

  * bearing_score distribution (mean, key quantiles, max)
  * raw p-value distribution (min, count and fraction with p < 0.05)
  * BH-adjusted distribution (min pval_adj_bh, count FDR-significant)
  * the p-value of each sample's single TOP-scoring bin (if even the strongest
    bin is non-significant, the null is too wide or scores are globally low)
  * the n_perms-implied p-value floor (if the test is per-bin permutation, the
    smallest attainable p is 1/(n_perms+1); pooled nulls go finer)

The goal is to localize a flat track to one of: (a) genuinely low/flat scores,
(b) an over-wide or miscalibrated null (high scores still get high p), or
(c) a p-value-floor / multiple-testing effect.

Optionally (--null-qcat) reads a sample's permutation-null qcat.bgz to report
the null score distribution directly (requires pysam).

ASCII-only. Pure stdlib + numpy.
"""
import argparse
import math
import os
import sys
import numpy as np


STAT_COLS_NEEDED = ("chrom", "start", "end", "bearing_score", "pval")


def load_stats(path):
    """Return dict of column arrays from a stats.tsv. Missing optional columns
    are returned as None."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        score_i = idx.get("bearing_score")
        pval_i = idx.get("pval")
        padj_i = idx.get("pval_adj_bh")
        sig_col = next((h for h in header if h.lower().startswith("significant_fdr")),
                       None)
        sig_i = idx.get(sig_col) if sig_col else None
        scores, pvals, padjs, sigs = [], [], [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if score_i is None or score_i >= len(f):
                continue
            try:
                scores.append(float(f[score_i]))
            except ValueError:
                continue
            pvals.append(_safe_float(f, pval_i))
            padjs.append(_safe_float(f, padj_i))
            if sig_i is not None and sig_i < len(f):
                sigs.append(f[sig_i].strip().lower() in ("1", "true", "yes", "t"))
            else:
                sigs.append(False)
    return {
        "score": np.array(scores, dtype=float),
        "pval": np.array(pvals, dtype=float),
        "padj": np.array(padjs, dtype=float),
        "sig": np.array(sigs, dtype=bool),
        "sig_col": sig_col,
    }


def _safe_float(fields, i):
    if i is None or i >= len(fields):
        return np.nan
    try:
        return float(fields[i])
    except ValueError:
        return np.nan


def fmt_q(a):
    a = a[np.isfinite(a)]
    if a.size == 0:
        return "n/a"
    qs = np.percentile(a, [0, 50, 90, 99, 100])
    return ("min=%.3g  med=%.3g  p90=%.3g  p99=%.3g  max=%.3g"
            % (qs[0], qs[1], qs[2], qs[3], qs[4]))


def report_sample(name, d, n_perms):
    score = d["score"]
    pval = d["pval"]
    padj = d["padj"]
    n = score.size
    print("=" * 72)
    print("Sample: %s   (%d scorable bins)" % (name, n))
    if n == 0:
        print("  no bins.")
        return
    print("  bearing_score : %s" % fmt_q(score))
    fin_p = pval[np.isfinite(pval)]
    if fin_p.size:
        n_raw = int(np.sum(fin_p < 0.05))
        print("  raw p-value   : min=%.4g   p<0.05: %d (%.2f%%)"
              % (np.min(fin_p), n_raw, 100.0 * n_raw / fin_p.size))
    fin_q = padj[np.isfinite(padj)]
    if fin_q.size:
        n_fdr = int(np.sum(fin_q < 0.05))
        print("  BH pval_adj   : min=%.4g   q<0.05: %d (%.3f%%)"
              % (np.min(fin_q), n_fdr, 100.0 * n_fdr / fin_q.size))
    n_sig = int(np.sum(d["sig"]))
    print("  FDR-significant bins (%s): %d" % (d["sig_col"] or "?", n_sig))

    # the single strongest bin: is it even close to significant?
    top = int(np.nanargmax(score))
    print("  TOP-scoring bin: score=%.3g  p=%.4g  q=%.4g"
          % (score[top],
             pval[top] if top < pval.size else float("nan"),
             padj[top] if top < padj.size else float("nan")))

    if n_perms:
        floor = 1.0 / (n_perms + 1)
        print("  per-bin perm p floor (n_perms=%d) = %.3g  "
              "(raw p never below this if the null is per-bin)" % (n_perms, floor))


def cross_compare(samples):
    print("=" * 72)
    print("CROSS-SAMPLE READOUT")
    # rank samples by median score and by FDR-sig count
    rows = []
    for name, d in samples:
        sc = d["score"][np.isfinite(d["score"])]
        rows.append((name, np.median(sc) if sc.size else float("nan"),
                     np.max(sc) if sc.size else float("nan"),
                     int(np.sum(d["sig"]))))
    print("%-14s  %10s  %10s  %10s" % ("sample", "med_score", "max_score", "FDR_sig"))
    for name, med, mx, nsig in rows:
        print("%-14s  %10.3g  %10.3g  %10d" % (name, med, mx, nsig))
    print("")
    print("Interpretation guide:")
    print("  * If the flat sample's med/max score is comparable to a working")
    print("    sample but FDR_sig=0 and its TOP bin has a large p -> the null")
    print("    is too wide / miscalibrated for that sample (calibration).")
    print("  * If its scores are globally lower -> low effect size (data/biology).")
    print("  * If raw p has many <0.05 but BH q has none -> multiple-testing /")
    print("    p-floor effect (compare n_perms and pooled-vs-per-bin null).")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pvalue-dir", required=True)
    ap.add_argument("--samples", nargs="+", required=True,
                    help="Sample names (e.g. ProB_rep1 ProB_rep2 DN_rep1).")
    ap.add_argument("--n-perms", type=int, default=None,
                    help="n_perms used, to show the per-bin p-value floor.")
    ap.add_argument("--null-qcat", default=None,
                    help="Optional permutation-null qcat.bgz for the flat sample "
                         "to report the null score distribution directly.")
    args = ap.parse_args()

    loaded = []
    for s in args.samples:
        path = os.path.join(args.pvalue_dir, "%s.stats.tsv" % s)
        if not os.path.exists(path):
            print("[WARN] missing %s" % path, file=sys.stderr)
            continue
        loaded.append((s, load_stats(path)))
    if not loaded:
        print("[ERROR] no stats files loaded", file=sys.stderr)
        sys.exit(1)

    for name, d in loaded:
        report_sample(name, d, args.n_perms)
    cross_compare(loaded)

    if args.null_qcat:
        print("=" * 72)
        print("Null score distribution from %s" % args.null_qcat)
        try:
            import pysam
            scores = []
            with pysam.TabixFile(args.null_qcat) as tb:
                for row in tb.fetch():
                    parts = row.split("\t")
                    # qcat rows: chrom start end <state scores...>; sum |kl| as
                    # the per-bin null score proxy if a total is not stored.
                    try:
                        vals = [float(x) for x in parts[3:]]
                        scores.append(sum(abs(v) for v in vals))
                    except ValueError:
                        continue
            arr = np.array(scores, dtype=float)
            print("  null per-bin score: %s" % fmt_q(arr))
        except ImportError:
            print("  pysam not available; skipping null read.")
        except Exception as e:
            print("  could not read null qcat: %s" % e)


if __name__ == "__main__":
    main()
