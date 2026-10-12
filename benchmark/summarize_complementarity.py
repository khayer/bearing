#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : benchmark/summarize_complementarity.py
# Version  : 1.0.0
# Date     : 2026-10-11
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""
Pool complementarity_benchmark.py runs (one directory per
prior<a>_depth<d>_seed<s>) into:

  --out             one row per prior x depth x block type: mean and range
                    across seeds of the call rate of production BEARING
                    (B_sum), BEARING with --diff-stat l1 (B_L1) and the
                    per-assay test, the per-assay discordant (switch-like)
                    sign-pattern share, and the empirical FDR of each test
                    (pooled false calls / pooled calls over seeds)
  --per-config-out  one row per prior x depth x configuration, means across
                    seeds

ASCII only. Reads the run outputs; fabricates nothing.
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def stats(v):
    v = [x for x in v if x == x]
    if not v:
        return "NA", "NA"
    return "%.3f" % (sum(v) / len(v)), "%.3f-%.3f" % (min(v), max(v))


def parse_fdr(summary_path):
    """{test: (false, total)} from the 'Empirical FDR' line of summary.txt."""
    out = {}
    with open(summary_path) as fh:
        for line in fh:
            if line.startswith("Empirical FDR"):
                for m in re.finditer(r"(\w+) (\d+)/(\d+)", line.split(":", 1)[1]):
                    out[m.group(1)] = (int(m.group(2)), int(m.group(3)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-config-out", required=True)
    a = ap.parse_args()

    runs = sorted(glob.glob(os.path.join(a.root, "prior*_depth*_seed*")))
    runs = [r for r in runs if os.path.isdir(r)]
    if not runs:
        sys.exit("no runs under %s" % a.root)
    pat = re.compile(r"prior([^_]+)_depth([^_]+)_seed(\d+)$")
    per_type = defaultdict(lambda: defaultdict(list))
    per_cfg = defaultdict(lambda: defaultdict(list))
    fdr = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    order_t, order_c = [], []
    for r in runs:
        m = pat.search(r)
        if not m or not os.path.exists(os.path.join(r, "per_type.tsv")):
            continue
        key = (m.group(1), m.group(2))
        for row in read(os.path.join(r, "per_type.tsv")):
            if row["type"] not in order_t:
                order_t.append(row["type"])
            for c in ("bearing_sum_called", "bearing_l1_called", "per_assay_called",
                      "per_assay_discordant"):
                per_type[key + (row["type"],)][c].append(fnum(row[c]))
        for row in read(os.path.join(r, "per_config.tsv")):
            if row["config"] not in order_c:
                order_c.append((row["config"]))
            k = key + (row["config"],)
            per_cfg[k]["type"] = row["type"]
            for c in ("E_stat_sum", "E_stat_l1", "bearing_sum_tested",
                      "bearing_sum_called", "bearing_l1_called", "per_assay_called",
                      "per_assay_discordant"):
                per_cfg[k].setdefault(c, []).append(fnum(row[c]))
        for t, (f, n) in parse_fdr(os.path.join(r, "summary.txt")).items():
            fdr[key][t][0] += f
            fdr[key][t][1] += n

    keys = sorted({k[:2] for k in per_type}, key=lambda k: (float(k[0]), float(k[1])))
    cols = ["prior_strength", "depth", "type", "n_seeds",
            "bearing_sum_called", "bearing_sum_range",
            "bearing_l1_called", "bearing_l1_range",
            "per_assay_called", "per_assay_range",
            "per_assay_discordant", "fdr_bearing_sum", "fdr_bearing_l1",
            "fdr_per_assay"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k in keys:
            fk = fdr[k]

            def fd(t):
                f, n = fk.get(t, [0, 0])
                return "%d/%d" % (f, n) + ((" (%.3f)" % (f / n)) if n else "")
            for ty in order_t:
                d = per_type.get(k + (ty,))
                if not d:
                    continue
                bs, bsr = stats(d["bearing_sum_called"])
                bl, blr = stats(d["bearing_l1_called"])
                pa, par = stats(d["per_assay_called"])
                pd, _ = stats(d["per_assay_discordant"])
                fh.write("\t".join([k[0], k[1], ty, str(len(d["per_assay_called"])),
                                    bs, bsr, bl, blr, pa, par, pd,
                                    fd("bearing_sum"), fd("bearing_l1"),
                                    fd("per_assay")]) + "\n")
    ccols = ["prior_strength", "depth", "config", "type", "E_stat_sum", "E_stat_l1",
             "bearing_sum_tested", "bearing_sum_called", "bearing_l1_called",
             "per_assay_called", "per_assay_discordant"]
    with open(a.per_config_out, "w") as fh:
        fh.write("\t".join(ccols) + "\n")
        for k in keys:
            for c in order_c:
                d = per_cfg.get(k + (c,))
                if not d:
                    continue
                vals = [stats(d[x])[0] for x in ccols[4:]]
                fh.write("\t".join([k[0], k[1], c, d["type"]] + vals) + "\n")
    print("wrote %s and %s (%d prior x depth cells, %d runs)"
          % (a.out, a.per_config_out, len(keys), len(runs)))


if __name__ == "__main__":
    main()
