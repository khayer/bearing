#!/usr/bin/env python3
"""
regional_null_calibration.py
============================
Empirical calibration of the regional-enrichment test on matched random regions.

WHY
The regional-enrichment test (regional_enrichment.py) combines a spatial
concentration binomial with a directional-concordance binomial via Fisher's
method, and assumes bins are exchangeable within a region. Bins are spatially
autocorrelated (the differential-score independence length is ~967 bp against a
200 bp grid), so exchangeability is violated. The manuscript ASSERTS the test is
conservative under this violation; this script measures it.

WHAT
Applies the SAME test, through the SAME production code path
(regional_enrichment.compute_regional_enrichment), to random regions drawn from
the same analysis locus and matched to each real region on the number of scored
bins. Under the null the combined p-values should be approximately Uniform(0,1).

Reported per real region (i.e. per size class):
  - empirical type-I error at alpha = 0.05 and 0.01 for p_spatial,
    p_directional and p_combined
  - KS distance from Uniform(0,1) (descriptive; see CAVEAT)
  - the real region's p_combined and its rank among the matched random regions

READING THE OUTPUT
  FPR@0.05 ~ 0.05     -> calibrated; drop the conservativeness claim.
  FPR@0.05 << 0.05    -> conservative, as claimed; now demonstrated.
  FPR@0.05 >> 0.05    -> anti-conservative; the exchangeability assumption fails
                         and the regional q-values must be revised.

CAVEAT
Random windows drawn from a single ~900 kb locus overlap one another, so they
are not independent draws. The rejection rate remains an unbiased estimate of
the type-I error for a region of that size in that locus, but the KS test's
p-value would not be valid. We therefore report the KS DISTANCE as a descriptive
statistic and attach no p-value to it.

USAGE (production settings; see Snakefile REGIONAL_LOCI / REGIONAL_PTHRESH)
  python regional_null_calibration.py \
      --stats   workflow/results/pvalue/diff_DN_vs_DP.stats.tsv \
      --regions workflow/annotations/tcrb_regions_v5.bed \
      --locus   chr6:40790000-41690000 \
      --p-thresh 0.05 --n-random 1000 --seed 42 \
      --out     regional_null_calibration_tcrb_DN_vs_DP.tsv

  # Igh:
  #   --regions workflow/annotations/igh_regions_v4.bed
  #   --locus   chr12:113000000-116100000

NOTE the .qcat.bgz cannot be used: regional_enrichment.parse_diff_qcat() opens
files with plain open(). The production diff_<cmp>.stats.tsv carries the
bearing_score column that parser expects, and batch mode is invoked on it, so a
single --stats file serves as both inputs.

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import bisect
import logging
import os
import sys

import numpy as np

logger = logging.getLogger(__name__)


def import_production(repo_root):
    sys.path.insert(0, repo_root)
    try:
        from regional_enrichment import (parse_diff_qcat, parse_diff_pvals,
                                         parse_regions_bed, parse_ucsc_region,
                                         compute_regional_enrichment)
    except Exception as e:
        sys.exit("could not import regional_enrichment from %s: %s\n"
                 "run from the repo root or pass --repo" % (repo_root, e))
    return (parse_diff_qcat, parse_diff_pvals, parse_regions_bed,
            parse_ucsc_region, compute_regional_enrichment)


def scored_starts(pvals, chrom, lo, hi):
    """Sorted start coordinates of every SCORED bin inside [lo,hi) on chrom."""
    return sorted(s for (c, s, _e) in pvals if c == chrom and lo <= s < hi)


def bins_in(starts, start, end):
    """Number of scored bins in [start,end) -- O(log n) via bisect."""
    return bisect.bisect_left(starts, end) - bisect.bisect_left(starts, start)


def sample_matched(starts, chrom, lo, hi, target_bins, n_want, rng,
                   avoid, max_tries=None):
    """Random windows inside [lo,hi) containing EXACTLY target_bins scored bins.

    Taking target_bins CONSECUTIVE scored starts guarantees the count exactly,
    so no rejection sampling on bin count is needed. (The previous version
    recounted with an O(N) scan per candidate; with 13.6M bins and up to 2e5
    tries that was ~1e12 operations.)
    """
    if target_bins == 0 or len(starts) <= target_bins:
        return []
    if max_tries is None:
        max_tries = n_want * 20          # repo convention (tcrb_contact_isolation.py)
    out, tries, seen = [], 0, set()
    while len(out) < n_want and tries < max_tries:
        tries += 1
        i = int(rng.integers(0, len(starts) - target_bins))
        if i in seen:
            continue
        seen.add(i)
        start = starts[i]
        end = starts[i + target_bins - 1] + 200
        if end > hi:
            continue
        if any(not (end <= a0 or start >= a1) for a0, a1 in avoid):
            continue                       # overlaps a real region
        out.append((chrom, start, end, "rand_%d" % len(out)))
    return out


def ks_distance_uniform(ps):
    """One-sample KS distance from Uniform(0,1). Descriptive only (see CAVEAT)."""
    if len(ps) == 0:
        return float("nan")
    xs = np.sort(np.asarray(ps, dtype=float))
    n = xs.size
    upper = np.arange(1, n + 1) / n - xs
    lower = xs - np.arange(0, n) / n
    return float(max(upper.max(), lower.max()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", default=None,
                    help="production diff_<cmp>.stats.tsv; used for BOTH the qcat "
                         "and p-value inputs (parse_diff_qcat reads stats files). "
                         "Preferred: the .qcat.bgz cannot be read by the production "
                         "parser, which uses plain open().")
    ap.add_argument("--diff-qcat", default=None)
    ap.add_argument("--diff-pvals", default=None)
    ap.add_argument("--regions", required=True, help="BED of the real regions")
    ap.add_argument("--locus", required=True, help="chrom:start-end, as in production")
    ap.add_argument("--p-thresh", type=float, default=0.05)
    ap.add_argument("--n-random", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    (parse_qcat, parse_pvals, parse_bed, parse_region,
     compute) = import_production(os.path.abspath(a.repo))

    if a.stats:
        a.diff_qcat = a.diff_pvals = a.stats
    if not (a.diff_qcat and a.diff_pvals):
        sys.exit("give --stats, or both --diff-qcat and --diff-pvals")
    if a.diff_qcat.endswith(".bgz") or a.diff_qcat.endswith(".gz"):
        sys.exit("regional_enrichment.parse_diff_qcat() uses plain open() and cannot "
                 "read a bgz file. Pass --stats <diff_*.stats.tsv> instead.")

    print("loading differential scores + p-values ...")
    qcat = parse_qcat(a.diff_qcat)
    pvals = parse_pvals(a.diff_pvals)
    real = parse_bed(a.regions)
    lchrom, lstart, lend = parse_region(a.locus)
    print("  %d qcat bins, %d p-value bins, %d real regions"
          % (len(qcat), len(pvals), len(real)))

    avoid = [(s, e) for (c, s, e, _n) in real if c == lchrom]
    rng = np.random.default_rng(a.seed)
    starts = scored_starts(pvals, lchrom, lstart, lend)     # built once
    print("  %d scored bins inside the analysis locus" % len(starts))

    # real regions first, through the production code path
    real_res = {r["region"] if "region" in r else r.get("name", "?"): r
                for r in compute(qcat, pvals, real, lchrom, lstart, lend, a.p_thresh)}

    rows = []
    print("\n%-28s %7s %10s %10s %10s %9s %8s"
          % ("region (size class)", "n_bins", "FPR@0.05", "FPR@0.01", "KS dist",
             "real p", "rank"))
    for (c, s, e, name) in real:
        if c != lchrom:
            continue
        nb = bins_in(starts, s, e)
        rand = sample_matched(starts, lchrom, lstart, lend, nb, a.n_random, rng, avoid)
        if len(rand) < 50:
            print("%-28s %7d   (only %d matched random regions; skipped)"
                  % (name[:28], nb, len(rand)))
            continue
        res = compute(qcat, pvals, rand, lchrom, lstart, lend, a.p_thresh)
        pc = [r["p_combined"] for r in res]
        ps = [r["p_spatial"] for r in res]
        pd = [r["p_directional"] for r in res]

        fpr05 = sum(1 for p in pc if p < 0.05) / len(pc)
        fpr01 = sum(1 for p in pc if p < 0.01) / len(pc)
        ksd = ks_distance_uniform(pc)

        rp = None
        for k, r in real_res.items():
            if str(k).startswith(name[:12]):
                rp = r["p_combined"]
                break
        rank = (sum(1 for p in pc if p <= rp) / len(pc)) if rp is not None else float("nan")

        print("%-28s %7d %10.4f %10.4f %10.3f %9s %8s"
              % (name[:28], nb, fpr05, fpr01, ksd,
                 ("%.3g" % rp) if rp is not None else "-",
                 ("%.4f" % rank) if rp is not None else "-"))
        for r, p_c, p_s, p_d in zip(rand, pc, ps, pd):
            rows.append((name, nb, r[1], r[2], p_s, p_d, p_c))

    with open(a.out, "w") as fh:
        fh.write("real_region\tn_bins\trand_start\trand_end\tp_spatial\tp_directional\tp_combined\n")
        for r in rows:
            fh.write("%s\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\n" % r)
    print("\nwrote %d random-region tests -> %s" % (len(rows), a.out))
    print("\nHow to read this:")
    print("  FPR@0.05 close to 0.05  -> the test is calibrated.")
    print("  FPR@0.05 well BELOW 0.05 -> conservative, as the manuscript claims.")
    print("  FPR@0.05 ABOVE 0.05      -> anti-conservative; the exchangeability")
    print("                              assumption fails and the claim must be revised.")
    print("  'rank' is the fraction of matched random regions with p <= the real")
    print("  region's p. A real finding should sit far in the tail (rank near 0).")


if __name__ == "__main__":
    main()
