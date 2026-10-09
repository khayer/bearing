#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : dev/regional_empirical_table.py
# Version  : 1.0.0
# Date     : 2026-10-09
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""
Regional-enrichment table with contrast-matched empirical p-values
(Table S2 / S11 at prior_strength 1).

Inputs
  --nominal  consolidated_enrichment_<locus>.tsv (regional_consolidate rule):
             one row per comparison x region with p_spatial, p_directional,
             p_combined and L_region_bins.
  --nulls    one regional_null_calibration.py output per comparison
             (real_region, n_bins, null_chrom, null_start, p_combined), named
             cmnull_<locus>_<comparison>.tsv.

Per comparison x region
  emp_count  number of null regions with p_combined <= the real p_combined
  emp_n      number of null regions drawn for that region
  emp_p      (emp_count + 1) / (emp_n + 1)        (add-one; floor 1/(n+1))
  testable   L_region_bins >= --min-bins (default 5, the same cut
             regional_null_calibration.py uses); untestable regions get no
             emp_p.
  q_testable BH over the testable regions of the locus (primary; agreed with
             the advisors 2026-10-09)
  q_all      BH over every comparison x region of the locus (the original
             20-test family; untestable regions enter with p = 1)

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import csv
import os
import re
import sys


def bh(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        q[i] = min(prev, 1.0)
    return q


def cmp_name(raw):
    c = raw
    if c.startswith("diff_"):
        c = c[len("diff_"):]
    if c.endswith(".stats"):
        c = c[:-len(".stats")]
    return c


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nominal", required=True)
    ap.add_argument("--nulls", nargs="+", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--min-bins", type=int, default=5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # comparison -> region -> list of null p_combined
    nulls = {}
    pat = re.compile(r"cmnull_%s_(.+)\.tsv$" % re.escape(a.locus))
    for path in a.nulls:
        m = pat.search(os.path.basename(path))
        if not m:
            sys.exit("cannot parse comparison from %s (expected cmnull_%s_<cmp>.tsv)"
                     % (path, a.locus))
        c = m.group(1)
        per = nulls.setdefault(c, {})
        for r in read_tsv(path):
            per.setdefault(r["real_region"], []).append(float(r["p_combined"]))

    rows = []
    for r in read_tsv(a.nominal):
        c = cmp_name(r["comparison"])
        reg = r["region_name"]
        nb = int(float(r["L_region_bins"]))
        p = float(r["p_combined"])
        testable = nb >= a.min_bins
        draws = nulls.get(c, {}).get(reg, [])
        if testable and not draws:
            sys.exit("no null draws for testable region %s / %s" % (c, reg))
        if testable:
            k = sum(1 for x in draws if x <= p)
            emp = (k + 1.0) / (len(draws) + 1.0)
        else:
            k, emp = None, None
        rows.append(dict(locus=a.locus, comparison=c, region=reg,
                         chrom=r["chrom"], start=r["start"], end=r["end"],
                         tested_bins=nb, k=r["k"], k_pos=r["k_pos"], k_neg=r["k_neg"],
                         n_locus=r["n_locus"], p_spatial=r["p_spatial"],
                         p_directional=r["p_directional"], p_combined=r["p_combined"],
                         g_dir=r["g_dir"], testable="yes" if testable else "no",
                         emp_n=len(draws) if testable else "",
                         emp_count=k if testable else "",
                         emp_p=emp))

    t_idx = [i for i, r in enumerate(rows) if r["testable"] == "yes"]
    q_t = bh([rows[i]["emp_p"] for i in t_idx])
    for i, q in zip(t_idx, q_t):
        rows[i]["q_testable"] = q
    q_a = bh([r["emp_p"] if r["emp_p"] is not None else 1.0 for r in rows])
    for r, q in zip(rows, q_a):
        r["q_all"] = q

    cols = ["locus", "comparison", "region", "chrom", "start", "end", "tested_bins",
            "k", "k_pos", "k_neg", "n_locus", "p_spatial", "p_directional",
            "p_combined", "g_dir", "testable", "emp_n", "emp_count", "emp_p",
            "q_testable", "q_all"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            out = []
            for c in cols:
                v = r.get(c, "")
                if v is None:
                    v = ""
                elif isinstance(v, float):
                    v = "%.6g" % v
                out.append(str(v))
            fh.write("\t".join(out) + "\n")
    n_sig = sum(1 for r in rows if isinstance(r.get("q_testable"), float)
                and r["q_testable"] < 0.05)
    print("%s: %d rows, %d testable, %d with q_testable < 0.05 -> %s"
          % (a.locus, len(rows), len(t_idx), n_sig, a.out))


if __name__ == "__main__":
    main()
