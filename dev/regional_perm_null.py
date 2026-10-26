#!/usr/bin/env python3
# regional_perm_null.py
# Authors: Katharina E. Hayer (Nina) with Claude (Opus 5.5), co-created
# Version: 1.0   Date: 2026-10-25
# Purpose: B10 option (c) - regional empirical test against the existing permutation rounds.
#   For each region in a BED (e.g. annotations/tcrb_regions_v5.bed):
#     tested bins  = bins of the observed comparison that appear in pvalue/diff_<cmp>.stats.tsv
#     region score = sum over the tested bins of the per-bin l1 diff score (sum_i |d_i| from the diff qcat)
#     the same score is computed on the SAME bins in every perm round (perm*/diff_comparison/diff_<cmp>.qcat.bgz)
#     emp p = (1 + #{perm score >= observed score}) / (1 + n_perm); BH over the regions in the BED.
#   This replaces the contrast-matched background null for the V1P 7-track run, whose 327 draws collapse
#   into 31 overlapping clusters on chr6 (A18 / B10 in claude/phd_open_tasks.md).
# Usage (from workflow/):
#   python $(git rev-parse --show-toplevel)/dev/regional_perm_null.py \
#     --results results_v1p7_prior1_l1 --cmp DN_vs_V1P \
#     --regions $(git rev-parse --show-toplevel)/annotations/tcrb_regions_v5.bed \
#     --window chr6:40790000-41690000 \
#     --out results_v1p7_prior1_l1/table_sources/regional_perm_null_tcrb_DN_vs_V1P.tsv \
#     --perm-out results_v1p7_prior1_l1/table_sources/regional_perm_null_tcrb_DN_vs_V1P_perms.tsv
# Needs: tabix-indexed .qcat.bgz files (pysam if installed, else the tabix binary on PATH).
# ASCII-only source.

import argparse
import glob
import gzip
import json
import os
import re
import subprocess
import sys

QCAT_RE = re.compile(r"qcat:(\[.*\])")


def parse_window(s):
    chrom, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return chrom, int(a), int(b)


def read_regions(path):
    regs = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            name = f[3] if len(f) > 3 else "%s:%s-%s" % (f[0], f[1], f[2])
            regs.append((name, f[0], int(f[1]), int(f[2])))
    return regs


def read_tested_bins(stats_path, chrom, lo, hi):
    """Set of bin starts in the window that appear in the per-bin stats table (= tested bins)."""
    op = gzip.open if stats_path.endswith(".gz") else open
    starts = set()
    with op(stats_path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if "chrom" not in head or "start" not in head:
            sys.exit("stats table has no 'chrom'/'start' header columns: %s" % head[:8])
        ci, si = head.index("chrom"), head.index("start")
        for line in fh:
            f = line.split("\t", max(ci, si) + 1)
            if f[ci] != chrom:
                continue
            s = int(f[si])
            if lo <= s < hi:
                starts.add(s)
    return starts


def fetch_lines(path, chrom, lo, hi):
    try:
        import pysam
        tb = pysam.TabixFile(path)
        try:
            return list(tb.fetch(chrom, lo, hi))
        finally:
            tb.close()
    except ImportError:
        out = subprocess.run(["tabix", path, "%s:%d-%d" % (chrom, lo + 1, hi)],
                             check=True, capture_output=True, text=True).stdout
        return [x for x in out.split("\n") if x]


def l1_by_bin(path, chrom, lo, hi):
    """bin start -> l1 score (sum of |d| over categories) for one diff qcat."""
    res = {}
    for line in fetch_lines(path, chrom, lo, hi):
        f = line.split("\t")
        m = QCAT_RE.search(f[3])
        if not m:
            continue
        pairs = json.loads(m.group(1))
        res[int(f[1])] = sum(abs(float(v)) for v, _c in pairs)
    return res


def region_scores(l1, regs, tested):
    out = []
    for name, c, s, e in regs:
        bins = [b for b in tested[name]]
        out.append(sum(l1.get(b, 0.0) for b in bins))
    return out


def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True, help="results folder, e.g. results_v1p7_prior1_l1")
    ap.add_argument("--cmp", required=True, help="comparison, e.g. DN_vs_V1P")
    ap.add_argument("--regions", required=True, help="regions BED (name in column 4)")
    ap.add_argument("--window", required=True, help="locus window chrom:start-end")
    ap.add_argument("--out", required=True)
    ap.add_argument("--perm-out", default=None, help="optional per-perm region scores")
    a = ap.parse_args()

    chrom, lo, hi = parse_window(a.window)
    regs = [r for r in read_regions(a.regions) if r[1] == chrom and r[2] >= lo and r[3] <= hi]
    if not regs:
        sys.exit("no regions of %s inside %s" % (a.regions, a.window))

    stats = os.path.join(a.results, "pvalue", "diff_%s.stats.tsv" % a.cmp)
    if not os.path.exists(stats) and os.path.exists(stats + ".gz"):
        stats += ".gz"
    obs_q = os.path.join(a.results, "compare", "diff_%s.qcat.bgz" % a.cmp)
    perms = sorted(glob.glob(os.path.join(a.results, "perm", "perm*", "diff_comparison",
                                          "diff_%s.qcat.bgz" % a.cmp)),
                   key=lambda p: int(re.search(r"perm(\d+)", p).group(1)))
    if not perms:
        sys.exit("no perm diff qcats found under %s/perm" % a.results)

    window_tested = read_tested_bins(stats, chrom, lo, hi)
    tested = {name: sorted(b for b in window_tested if s <= b < e) for name, c, s, e in regs}

    obs_l1 = l1_by_bin(obs_q, chrom, lo, hi)
    obs = region_scores(obs_l1, regs, tested)
    missing = sum(1 for n in tested for b in tested[n] if b not in obs_l1)
    if missing:
        print("WARNING: %d tested bins have no observed qcat line (scored 0)" % missing)

    null = []  # list of per-perm lists
    for i, p in enumerate(perms, 1):
        null.append(region_scores(l1_by_bin(p, chrom, lo, hi), regs, tested))
        if i % 20 == 0:
            print("  %d / %d perm rounds read" % (i, len(perms)))
    n = len(null)

    emp = []
    rows = []
    for j, (name, c, s, e) in enumerate(regs):
        col = [r[j] for r in null]
        n_ge = sum(1 for x in col if x >= obs[j])
        p = (1.0 + n_ge) / (1.0 + n)
        emp.append(p)
        mean = sum(col) / n
        sd = (sum((x - mean) ** 2 for x in col) / max(n - 1, 1)) ** 0.5
        rows.append([name, c, s, e, len(tested[name]), obs[j], mean, sd, max(col), n, n_ge, p])
    q = bh(emp)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("region\tchrom\tstart\tend\ttested_bins\tobs_score_l1sum\tnull_mean\tnull_sd\t"
                 "null_max\tn_perm\tn_ge\temp_p\tq_bh\n")
        for r, qq in zip(rows, q):
            fh.write("%s\t%s\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%d\t%.6g\t%.6g\n"
                     % tuple(r + [qq]))
    if a.perm_out:
        with open(a.perm_out, "w") as fh:
            fh.write("perm\t" + "\t".join(r[0] for r in regs) + "\n")
            for p, r in zip(perms, null):
                tag = re.search(r"(perm\d+)", p).group(1)
                fh.write(tag + "\t" + "\t".join("%.6g" % x for x in r) + "\n")

    print("\n%d perm rounds; emp p floor = 1/%d = %.4g" % (n, n + 1, 1.0 / (n + 1)))
    print("region\ttested\tobs\tnull_mean\tnull_max\tn_ge\temp_p\tq_bh")
    for r, qq in zip(rows, q):
        print("%s\t%d\t%.4g\t%.4g\t%.4g\t%d\t%.4g\t%.4g" % (r[0], r[4], r[5], r[6], r[8], r[10], r[11], qq))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
