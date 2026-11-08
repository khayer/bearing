#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : benchmark/summarize_transform_sweep.py
# Version  : 1.0.0
# Date     : 2026-11-08
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""Collect run_transform_sweep.sh cells into one table.

One row per transform x rna_tail x depth (seeds pooled: mean, min, max) with
the l1 call rate of block groups, the empirical FDR of l1 and of the per-assay
test, the l1 tested-bin count and the RNA-led share of the l1 calls.
Groups (fraction of block-interior bins called by BEARING l1):
  magnitude, switch, rebalance, both, onoff   per_type.tsv
  arch_onoff_s0.6 / s0.9   onoff_{gain,loss}_arch_s*
  txn_onoff_s0.9           onoff_{gain,loss}_txn_s0.9
  rna_solo_onoff_s0.9      onoff_{gain,loss}_solo_rnaplus_s0.9
  ctcf_rad21_s0.6 / s0.9   switch_ctcf_to_rad21_s*
  arch_txn_s0.6 / s0.9     switch_arch_to_txn_s*
SIMULATED data only. ASCII only.
"""

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

CFG_GROUPS = {
    "arch_onoff_s0.6": ("onoff_gain_arch_s0.6", "onoff_loss_arch_s0.6"),
    "arch_onoff_s0.9": ("onoff_gain_arch_s0.9", "onoff_loss_arch_s0.9"),
    "txn_onoff_s0.9": ("onoff_gain_txn_s0.9", "onoff_loss_txn_s0.9"),
    "rna_solo_onoff_s0.9": ("onoff_gain_solo_rnaplus_s0.9",
                            "onoff_loss_solo_rnaplus_s0.9"),
    "ctcf_rad21_s0.6": ("switch_ctcf_to_rad21_s0.6",),
    "ctcf_rad21_s0.9": ("switch_ctcf_to_rad21_s0.9",),
    "arch_txn_s0.6": ("switch_arch_to_txn_s0.6",),
    "arch_txn_s0.9": ("switch_arch_to_txn_s0.9",),
}
TYPE_GROUPS = ("magnitude", "switch", "rebalance", "both", "onoff")
CELL = re.compile(r"tf_(?P<tf>[a-z0-9_]+)_tail(?P<tail>[0-9.]+)_d(?P<d>[0-9.]+)_s(?P<s>\d+)$")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def one_cell(od):
    out = {}
    for r in read_tsv(os.path.join(od, "per_type.tsv")):
        if r["type"] in TYPE_GROUPS:
            out[r["type"]] = float(r["bearing_l1_called"])
    cfg = {r["config"]: r for r in read_tsv(os.path.join(od, "per_config.tsv"))}
    for g, names in CFG_GROUPS.items():
        n = sum(int(cfg[c]["n_bins"]) for c in names)
        k = sum(float(cfg[c]["bearing_l1_called"]) * int(cfg[c]["n_bins"]) for c in names)
        out[g] = k / n if n else float("nan")
    txt = open(os.path.join(od, "summary.txt")).read()
    m = re.search(r"bearing_l1 (\d+)/(\d+)", txt)
    out["l1_false"], out["l1_calls"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"per_assay (\d+)/(\d+)", txt)
    out["pa_fdr"] = int(m.group(1)) / int(m.group(2)) if int(m.group(2)) else float("nan")
    m = re.search(r"--diff-stat l1 sum\|d\|\s+(\d+) bins tested", txt)
    out["l1_tested"] = int(m.group(1))
    m = re.search(r"\(RNA ([0-9.]+)%\)", txt)
    out["rna_lead_pct"] = float(m.group(1)) if m else float("nan")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    cells = {}
    for od in sorted(glob.glob(os.path.join(a.root, "tf_*"))):
        m = CELL.search(os.path.basename(od))
        if not m or not os.path.exists(os.path.join(od, "per_type.tsv")):
            continue
        key = (m.group("tf"), float(m.group("tail")), float(m.group("d")))
        cells.setdefault(key, []).append(one_cell(od))
    if not cells:
        sys.exit("ERROR: no finished cells under %s" % a.root)
    cols = list(TYPE_GROUPS) + list(CFG_GROUPS) + ["rna_lead_pct", "l1_tested", "pa_fdr"]
    order = {"none": 0, "nzq": 1, "log1p_rna": 2, "asinh_auto": 3}
    rows = []
    for key in sorted(cells, key=lambda k: (k[1], order.get(k[0], 9), k[2])):
        cs = cells[key]
        r = dict(transform=key[0], rna_tail_sigma=key[1], depth=key[2], seeds=len(cs),
                 l1_false_calls="%d/%d" % (sum(c["l1_false"] for c in cs),
                                           sum(c["l1_calls"] for c in cs)))
        for c in cols:
            v = np.array([x[c] for x in cs], dtype=float)
            r[c] = "%.3g (%.3g-%.3g)" % (np.nanmean(v), np.nanmin(v), np.nanmax(v))
        rows.append(r)
    head = ["transform", "rna_tail_sigma", "depth", "seeds", "l1_false_calls"] + cols
    with open(a.out, "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[h]) for h in head) + "\n")
    sys.stdout.write(open(a.out).read())


if __name__ == "__main__":
    main()
