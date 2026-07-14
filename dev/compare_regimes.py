#!/usr/bin/env python3
"""Summarize significance-set concordance of several BEARING scoring regimes
against a single reference run, in one table.

Each regime is compared to the reference the same way compare_default_vs_qnorm.py
does pairwise: per differential comparison, the Jaccard overlap of the
FDR-significant bin SETS (unsigned and direction-matched) and the Spearman
correlation of per-bin BEARING scores over shared bins. This wrapper just runs
that across every regime and collects it into one table plus a per-regime
summary, so the robustness suite is a single command instead of N pairwise runs.

Read the seed regime's Jaccard as the permutation-noise FLOOR: the scoring is
deterministic, so a same-settings different-seed run differs only by Monte-Carlo
sampling of the null. Metric / normalization / n_perms regimes are robust if
their Jaccard is at or above that floor.

Example (perm10 suite vs the seed42 reference, restricted to Tcrb):
    compare_regimes.py \\
      --ref     seed42=workflow/results_seed42/pvalue \\
      --regime  seed444=workflow/results_seed444/pvalue \\
      --regime  jsd=workflow/results_jsd/pvalue \\
      --regime  qnorm=workflow/results_qnorm/pvalue \\
      --regime  perm100=workflow/results/pvalue \\
      --region  chr6:40000000-42000000 \\
      --out     regime_concordance_tcrb.tsv

Reuses load_stats / sig_keyset / jaccard / spearman / parse_region from
compare_default_vs_qnorm.py. No third-party deps beyond numpy. ASCII-only.
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

from compare_default_vs_qnorm import (
    load_stats, sig_keyset, jaccard, spearman, parse_region,
)


def _split_label(spec):
    """'label=path' -> (label, path); 'path' -> (basename, path)."""
    if "=" in spec:
        label, path = spec.split("=", 1)
        return label, path
    return os.path.basename(os.path.normpath(spec)), spec


def concordance_pair(ref_dir, reg_dir, region=None):
    """Per-comparison concordance of reg_dir vs ref_dir. Returns list of
    (cmp, sig_ref, sig_reg, jacc, jacc_dir, rho, n_shared)."""
    files = sorted(glob.glob(os.path.join(ref_dir, "diff_*.stats.tsv")))
    cmps = [os.path.basename(p)[len("diff_"):-len(".stats.tsv")] for p in files]
    out = []
    for cmp in cmps:
        pr = os.path.join(ref_dir, "diff_%s.stats.tsv" % cmp)
        pg = os.path.join(reg_dir, "diff_%s.stats.tsv" % cmp)
        if not os.path.exists(pg):
            continue
        dr = load_stats(pr, region=region)
        dg = load_stats(pg, region=region)
        if dr is None or dg is None:
            continue
        sr, sg = sig_keyset(dr), sig_keyset(dg)
        srd, sgd = sig_keyset(dr, signed=True), sig_keyset(dg, signed=True)
        j = jaccard(sr, sg)
        jd = jaccard(srd, sgd)
        md = {k: v for k, v in zip(dr["key"], dr["score"])}
        shared = [(md[k], v) for k, v in zip(dg["key"], dg["score"]) if k in md]
        if shared:
            a = np.array([x for x, _ in shared])
            b = np.array([y for _, y in shared])
            rho = spearman(a, b)
        else:
            rho = float("nan")
        out.append((cmp, int(dr["sig"].sum()), int(dg["sig"].sum()),
                    j, jd, rho, len(shared)))
    return out


def _median(xs):
    xs = [x for x in xs if x == x]  # drop NaN
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True,
                    help="reference run: label=path to its pvalue dir.")
    ap.add_argument("--regime", action="append", required=True,
                    help="regime to compare vs ref: label=path (repeatable).")
    ap.add_argument("--region", default=None,
                    help="restrict to a region: 'chr6:40000000-42000000' or 'chr6'. "
                         "Omit for genome-wide.")
    ap.add_argument("--out", default=None, help="optional TSV path for the full table.")
    args = ap.parse_args()

    ref_label, ref_dir = _split_label(args.ref)
    if not os.path.isdir(ref_dir):
        sys.exit("[ERROR] reference dir not found: %s" % ref_dir)
    region = parse_region(args.region)
    rtxt = ("genome-wide" if region is None
            else (region[0] if region[1] is None else "%s:%d-%d" % region))

    print("=" * 86)
    print("Regime concordance vs reference '%s'  |  scope: %s" % (ref_label, rtxt))
    print("-" * 86)
    print("%-12s %-16s %8s %8s %8s %8s %8s %8s"
          % ("regime", "comparison", "sig_ref", "sig_reg", "jacc", "jacc_dir",
             "rho_sc", "n_shared"))

    all_rows = []          # (regime, cmp, sig_ref, sig_reg, jacc, jacc_dir, rho, n_shared)
    summary = []           # (regime, median_jacc, median_jacc_dir, median_rho, n_cmp)
    for spec in args.regime:
        label, path = _split_label(spec)
        if not os.path.isdir(path):
            print("%-12s   (dir not found: %s)" % (label, path))
            continue
        pairs = concordance_pair(ref_dir, path, region=region)
        for (cmp, sr, sg, j, jd, rho, ns) in pairs:
            print("%-12s %-16s %8d %8d %8.3f %8.3f %8.3f %8d"
                  % (label, cmp, sr, sg, j, jd, rho, ns))
            all_rows.append((label, cmp, sr, sg, j, jd, rho, ns))
        if pairs:
            summary.append((label,
                            _median([p[3] for p in pairs]),     # median jacc
                            _median([p[4] for p in pairs]),     # median jacc_dir
                            _median([p[5] for p in pairs]),     # median rho
                            len(pairs)))
        print("-" * 86)

    print("\nPER-REGIME SUMMARY (median across comparisons):")
    print("%-12s %10s %10s %10s %6s" %
          ("regime", "med_jacc", "med_jaccdir", "med_rho", "n_cmp"))
    for label, mj, mjd, mr, nc in summary:
        print("%-12s %10.3f %10.3f %10.3f %6d" % (label, mj, mjd, mr, nc))
    print("\nRead metric/normalization/n_perms regimes against the seed regime's")
    print("median Jaccard -- that is the permutation-noise floor. At or above it")
    print("=> robust; well below it => the regime changes the significant SET.")

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["regime", "comparison", "sig_ref", "sig_regime",
                        "jaccard", "jaccard_dir", "rho_score", "n_shared", "scope"])
            for (label, cmp, sr, sg, j, jd, rho, ns) in all_rows:
                w.writerow([label, cmp, sr, sg, "%.4f" % j, "%.4f" % jd,
                            "%.4f" % rho, ns, rtxt])
        print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
