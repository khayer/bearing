#!/usr/bin/env python3
"""Compare two binwise dumps on a shared grid, per bin.

Tests whether two perturbations change the SAME bins in the SAME direction.
If dump A (e.g. V1TxS) phenocopies dump B (e.g. V1P), the per-bin BES and
delta_contact should be positively correlated between them.

Both dumps must be the SAME grid (same chrom/start/end), e.g. both 4-track
native dumps. Usage:

  python compare_dumps.py --a mcc_results_v1txs/binwise_DN_vs_V1TxS.tsv \
                          --b mcc_results_v1_4track/binwise_DN_vs_V1P.tsv \
                          --region chr6:40800000-41650000

Reports Spearman(bes_a, bes_b) and Spearman(delta_a, delta_b) over the
overlapping bins, dropping non-finite values. ASCII only.
"""
import argparse
import sys
import numpy as np


def load(path):
    rows = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[ci["chrom"]], int(f[ci["start"]]), int(f[ci["end"]]))
            try:
                bes = float(f[ci["bes"]])
                dc = float(f[ci["delta_contact"]])
            except (ValueError, KeyError):
                continue
            rows[key] = (bes, dc)
    return rows


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        return float("nan"), 0
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if denom == 0:
        return float("nan"), x.size
    return float((rx * ry).sum() / denom), x.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="dump A (e.g. V1TxS)")
    ap.add_argument("--b", required=True, help="dump B (e.g. V1P)")
    ap.add_argument("--region", default=None, help="chrom:start-end to restrict")
    args = ap.parse_args()

    A = load(args.a)
    B = load(args.b)

    lo = hi = None
    chrom = None
    if args.region:
        chrom, rng = args.region.split(":")
        lo, hi = (int(x.replace(",", "")) for x in rng.split("-"))

    keys = [k for k in A if k in B]
    if chrom is not None:
        keys = [k for k in keys if k[0] == chrom and k[1] >= lo and k[2] <= hi]

    if not keys:
        sys.stderr.write("ERROR: no overlapping bins\n")
        sys.exit(1)

    bes_a = [A[k][0] for k in keys]
    bes_b = [B[k][0] for k in keys]
    dc_a = [A[k][1] for k in keys]
    dc_b = [B[k][1] for k in keys]

    r_bes, n1 = spearman(bes_a, bes_b)
    r_dc, n2 = spearman(dc_a, dc_b)

    print("region        : %s" % (args.region or "all"))
    print("overlap bins  : %d" % len(keys))
    print("rho(BES_a, BES_b)            : %+.4f  (n=%d)" % (r_bes, n1))
    print("rho(dContact_a, dContact_b)  : %+.4f  (n=%d)" % (r_dc, n2))
    print()
    print("Interpretation: high positive on BOTH -> A changes the same bins")
    print("the same way as B (phenocopy). Near zero -> A's signal is in")
    print("different bins (not a phenocopy).")


if __name__ == "__main__":
    main()
