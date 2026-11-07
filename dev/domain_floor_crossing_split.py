#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : dev/domain_floor_crossing_split.py
# Version  : 1.0.0
# Date     : 2026-11-07
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Reviewer back-pocket check: in a given results folder, how much of the Tcrb
# domain differential sits in FLOOR-CROSSING bins, and what is the per-track
# composition of the differential in floor-crossing vs other bins?
#
# Question it answers: "Before the composition prior, cohesin + CTCF carried
# 55-71% of the domain differential (v89). Was that share carried by
# floor-crossing (noise-composition) bins?"
#
# Definitions (same as claude/floor_crossing_and_swapnull_2026-09-28.md):
#   floor crossing = one condition has per-sample total score exactly 0 in BOTH
#   replicates (raw signal below min_signal -> zeroed) while the other condition
#   has signal in at least one replicate.
#   domain bins = bins in --region with |bearing_score| >= --min-abs (v89
#   decomposition rule: |BES| >= 0.5). Use --col to filter on another column
#   (e.g. bearing_score_tested for the l1 run).
#   share = sum over bins of |kl_i| / sum over bins and tracks of |kl_j|.
#
# Inputs (read-only): <results>/pvalue/diff_<A>_vs_<B>.stats.tsv and the per-
# sample qcats <results>/<sample>.qcat.bgz (+ .tbi). Nothing is written except
# --out (optional TSV).
#
# Usage:
#   python dev/domain_floor_crossing_split.py --results workflow/results_seed42 \
#       --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3
#   (repeat with workflow/results for the v89 run and with
#    workflow/results_prior1_n100_l1 --col bearing_score_tested as a check)
# ASCII only.
import argparse
import csv
import os
import re
import sys

import pysam

TRACKS = ["ATAC", "RNA+", "RNA-", "CTCF", "Cohesin", "H3K27ac"]  # categories/mm10_6track_panel.yaml
PAIR = re.compile(r"\[\s*(-?[0-9.eE+-]+)\s*,\s*(\d+)\s*\]")


def qcat_totals(path, chrom, start, end):
    """bin start -> total |score| for one sample over the region."""
    if not os.path.exists(path):
        sys.exit("missing qcat: %s" % path)
    out = {}
    with pysam.TabixFile(path) as tb:
        for line in tb.fetch(chrom, start, end):
            f = line.split("\t")
            vals = PAIR.findall(f[3]) if len(f) > 3 else []
            out[int(f[1])] = sum(abs(float(v)) for v, _ in vals)
    return out


def parse_region(s):
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, int(a), int(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--comparisons", nargs="+",
                    default=["DN_vs_DP", "DN_vs_EbKO", "DN_vs_ProB", "DN_vs_S3T3"])
    ap.add_argument("--region", default="chr6:40400000-42400000",
                    help="Tcrb recombination domain (default)")
    ap.add_argument("--col", default="bearing_score",
                    help="column used for the bin filter (default bearing_score = BES)")
    ap.add_argument("--min-abs", type=float, default=0.5)
    ap.add_argument("--reps", nargs=2, default=["rep1", "rep2"])
    ap.add_argument("--out", default=None, help="optional TSV")
    a = ap.parse_args()
    chrom, start, end = parse_region(a.region)

    rows = []
    print("results: %s | region %s | filter |%s| >= %g" % (a.results, a.region, a.col, a.min_abs))
    for cmp_ in a.comparisons:
        ca, cb = cmp_.split("_vs_")
        stats = os.path.join(a.results, "pvalue", "diff_%s.stats.tsv" % cmp_)
        if not os.path.exists(stats):
            print("%s: no stats table (%s)" % (cmp_, stats))
            continue
        tot = {}
        for cond in (ca, cb):
            for r in a.reps:
                s = "%s_%s" % (cond, r)
                tot[s] = qcat_totals(os.path.join(a.results, s + ".qcat.bgz"), chrom, start, end)
        acc = {"fc": [0.0] * 6, "other": [0.0] * 6}
        n = {"fc": 0, "other": 0}
        missing = 0
        with open(stats) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            kcols = ["kl_%d" % i for i in range(1, 7)]
            for k in kcols + [a.col]:
                if k not in rd.fieldnames:
                    sys.exit("%s: column %s not in %s" % (cmp_, k, stats))
            ccol = rd.fieldnames[0]
            scol = rd.fieldnames[1]
            for r in rd:
                if r[ccol] != chrom:
                    continue
                b0 = int(float(r[scol]))
                if b0 < start or b0 >= end:
                    continue
                try:
                    v = float(r[a.col])
                except ValueError:
                    continue
                if abs(v) < a.min_abs:
                    continue
                t = [tot[f"{c}_{rp}"].get(b0) for c in (ca, cb) for rp in a.reps]
                if any(x is None for x in t):
                    missing += 1
                    t = [x or 0.0 for x in t]
                a_zero = t[0] == 0 and t[1] == 0
                b_zero = t[2] == 0 and t[3] == 0
                fc = (a_zero and (t[2] > 0 or t[3] > 0)) or (b_zero and (t[0] > 0 or t[1] > 0))
                key = "fc" if fc else "other"
                n[key] += 1
                for i, k in enumerate(kcols):
                    acc[key][i] += abs(float(r[k] or 0))
        all_ = [x + y for x, y in zip(acc["fc"], acc["other"])]
        tot_all = sum(all_) or 1.0
        print("\n## %s: %d bins (%d floor-crossing, %d other)%s"
              % (cmp_, n["fc"] + n["other"], n["fc"], n["other"],
                 ("; %d bins absent from a qcat (treated as 0)" % missing) if missing else ""))
        print("  share of the summed |kl| carried by floor-crossing bins: %.1f%%"
              % (100 * sum(acc["fc"]) / tot_all))
        print("  %-8s %8s %8s %8s" % ("track", "all", "fc", "other"))
        for i, tr in enumerate(TRACKS):
            sh = [100 * v[i] / (sum(v) or 1.0) for v in (all_, acc["fc"], acc["other"])]
            print("  %-8s %7.1f%% %7.1f%% %7.1f%%" % (tr, *sh))
        arch = [100 * (v[3] + v[4]) / (sum(v) or 1.0) for v in (all_, acc["fc"], acc["other"])]
        print("  %-8s %7.1f%% %7.1f%% %7.1f%%" % ("CTCF+Coh", *arch))
        rows.append([cmp_, n["fc"] + n["other"], n["fc"], round(100 * sum(acc["fc"]) / tot_all, 2)]
                    + [round(x, 2) for x in arch])
    if a.out and rows:
        with open(a.out, "w") as fh:
            fh.write("comparison\tbins\tfc_bins\tfc_share_of_differential_pct\t"
                     "ctcf_cohesin_all_pct\tctcf_cohesin_fc_pct\tctcf_cohesin_other_pct\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
        print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
