#!/usr/bin/env python3
"""Multi-regime robustness summary for BEARING differential calls.

Compares each scoring/sampling regime against ONE reference run and emits a
single table, so the regime axes (permutation seed, permutation depth, metric,
normalization) can be read side by side instead of as separate pairwise reports.

It reuses compare_default_vs_qnorm.py's per-comparison math verbatim (FDR-
significant bin-set Jaccard, direction-aware Jaccard, Spearman of per-bin scores
over shared bins), so numbers match that tool exactly. The reference is normally
the perm10 default (results_seed42/pvalue); include the matched seed-variant
(e.g. seed444) as one regime -- its Jaccard is the permutation-noise FLOOR, and
the metric/normalization/perm-depth Jaccards should be read against it.

Usage:
  robustness_summary.py \
    --ref seed42=workflow/results_seed42/pvalue \
    --regime seed444=workflow/results_seed444/pvalue \
    --regime jsd=workflow/results_jsd/pvalue \
    --regime qnorm=workflow/results_qnorm/pvalue \
    --regime perm100=workflow/results/pvalue \
    [--region chr6:40000000-42000000] \
    --out robustness_summary.tsv

Each DIR is a run's pvalue/ folder (the one holding diff_<cmp>.stats.tsv).
The reference defines the comparison set; regimes are matched by comparison name.
ASCII only.
"""
import argparse
import glob
import os
import sys

import numpy as np

# Reuse the reference tool's exact per-comparison math.
from compare_default_vs_qnorm import (
    load_stats, sig_keyset, jaccard, spearman, parse_region,
)


def parse_labeled_dir(s):
    if "=" not in s:
        sys.exit("[ERROR] expected LABEL=DIR, got: %s" % s)
    label, d = s.split("=", 1)
    label, d = label.strip(), d.strip()
    if not label or not d:
        sys.exit("[ERROR] empty label or dir in: %s" % s)
    return label, d


def comparisons_in(pvalue_dir):
    files = sorted(glob.glob(os.path.join(pvalue_dir, "diff_*.stats.tsv")))
    return [os.path.basename(p)[len("diff_"):-len(".stats.tsv")] for p in files]


def concordance_one(ref_dir, reg_dir, cmp, region):
    """Per-comparison concordance of reg vs ref. Returns a dict or None."""
    pr = os.path.join(ref_dir, "diff_%s.stats.tsv" % cmp)
    pg = os.path.join(reg_dir, "diff_%s.stats.tsv" % cmp)
    if not (os.path.exists(pr) and os.path.exists(pg)):
        return None
    dr, dg = load_stats(pr, region=region), load_stats(pg, region=region)
    if dr is None or dg is None:
        return None
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
    return {
        "cmp": cmp,
        "ref_sig": int(dr["sig"].sum()),
        "reg_sig": int(dg["sig"].sum()),
        "jacc": j,
        "jacc_dir": jd,
        "rho": rho,
        "n_shared": len(shared),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, metavar="LABEL=DIR",
                    help="reference run's pvalue/ dir (defines the comparison set).")
    ap.add_argument("--regime", required=True, action="append", metavar="LABEL=DIR",
                    help="a regime's pvalue/ dir; repeatable.")
    ap.add_argument("--region", default=None,
                    help="restrict to a region ('chr6:40000000-42000000' or 'chr6'); "
                         "omit for genome-wide.")
    ap.add_argument("--out", default=None, help="per-(regime,comparison) TSV path.")
    args = ap.parse_args()

    ref_label, ref_dir = parse_labeled_dir(args.ref)
    regimes = [parse_labeled_dir(s) for s in args.regime]
    region = parse_region(args.region)

    if not os.path.isdir(ref_dir):
        sys.exit("[ERROR] reference dir not found: %s" % ref_dir)
    cmps = comparisons_in(ref_dir)
    if not cmps:
        sys.exit("[ERROR] no diff_*.stats.tsv in reference dir: %s" % ref_dir)

    rtxt = "genome-wide" if region is None else (
        region[0] if region[1] is None else "%s:%d-%d" % region)
    print("=" * 84)
    print("BEARING robustness summary  |  ref = %s  |  scope: %s" % (ref_label, rtxt))
    print("-" * 84)

    all_rows = []   # (regime, per-cmp dict)
    per_regime = []  # aggregate
    for reg_label, reg_dir in regimes:
        if not os.path.isdir(reg_dir):
            print("  %-12s (dir not found: %s)" % (reg_label, reg_dir))
            continue
        recs = []
        for cmp in cmps:
            r = concordance_one(ref_dir, reg_dir, cmp, region)
            if r is not None:
                recs.append(r)
                all_rows.append((reg_label, r))
        if not recs:
            print("  %-12s (no matching comparisons)" % reg_label)
            continue
        jaccs = np.array([r["jacc"] for r in recs])
        jdirs = np.array([r["jacc_dir"] for r in recs])
        rhos = np.array([r["rho"] for r in recs], dtype=float)
        per_regime.append({
            "regime": reg_label, "n_cmp": len(recs),
            "sum_ref_sig": sum(r["ref_sig"] for r in recs),
            "sum_reg_sig": sum(r["reg_sig"] for r in recs),
            "mean_jacc": float(np.mean(jaccs)), "min_jacc": float(np.min(jaccs)),
            "mean_jacc_dir": float(np.mean(jdirs)),
            "mean_rho": float(np.nanmean(rhos)) if np.isfinite(rhos).any() else float("nan"),
        })

    hdr = ("%-12s %5s %11s %11s %9s %8s %12s %8s"
           % ("regime", "n_cmp", "sum_ref_sig", "sum_reg_sig",
              "mean_jacc", "min_jacc", "mean_jaccDir", "mean_rho"))
    print(hdr)
    for a in per_regime:
        print("%-12s %5d %11d %11d %9.3f %8.3f %12.3f %8.3f"
              % (a["regime"], a["n_cmp"], a["sum_ref_sig"], a["sum_reg_sig"],
                 a["mean_jacc"], a["min_jacc"], a["mean_jacc_dir"], a["mean_rho"]))
    print("-" * 84)
    print("jacc = FDR-significant bin-set overlap vs ref; jaccDir also requires same")
    print("direction; rho = Spearman of per-bin scores over shared bins. Read the")
    print("metric/normalization/perm-depth rows against the matched seed regime,")
    print("whose Jaccard is the permutation-noise floor (not a real effect).")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("regime\tcomparison\tref_sig\tregime_sig\tjaccard\t"
                     "jaccard_dir\trho_score\tn_shared\n")
            for reg_label, r in all_rows:
                fh.write("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\n"
                         % (reg_label, r["cmp"], r["ref_sig"], r["reg_sig"],
                            r["jacc"], r["jacc_dir"], r["rho"], r["n_shared"]))
        print("\nwrote per-comparison detail -> %s" % args.out)


if __name__ == "__main__":
    main()
