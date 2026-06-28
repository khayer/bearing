#!/usr/bin/env python3
"""Excess-over-null tail diagnostic for binwise dumps.

Uses one dump as the empirical NULL distribution of |delta_contact| and asks how
far each TEST dump exceeds the null's 95th/99th percentiles, plus how the test's
largest changes concentrate in the RC window. ASCII only.

  python tail_diag_generic.py --null NULL.tsv --test A.tsv B.tsv \
        --region chr6:40800000-41650000 [--rc 41450000-41650000]
"""
import argparse
import numpy as np


def load(path, chrom, lo, hi):
    starts, vals = [], []
    with open(path) as f:
        next(f)
        for ln in f:
            c = ln.split("\t")
            if c[0] == chrom and int(c[1]) >= lo and int(c[2]) <= hi:
                starts.append(int(c[1]))
                vals.append(abs(float(c[4])))
    return np.array(starts), np.array(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", required=True)
    ap.add_argument("--test", nargs="+", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--rc", default="41450000-41650000")
    args = ap.parse_args()

    chrom, rng = args.region.split(":")
    lo, hi = (int(x) for x in rng.split("-"))
    rc_lo, rc_hi = (int(x) for x in args.rc.split("-"))

    _, nv = load(args.null, chrom, lo, hi)
    p95, p99 = np.percentile(nv, 95), np.percentile(nv, 99)
    import os
    print("NULL = %s" % os.path.basename(args.null))
    print("  p50=%.4f  p95=%.4f  p99=%.4f  (n=%d)" % (
        np.percentile(nv, 50), p95, p99, len(nv)))
    print("  %-22s %%>p95  %%>p99    max     top20_in_RC" % "test")
    rows = [("(null itself)", args.null)] + [("", t) for t in args.test]
    seen = set()
    for _, path in rows:
        if path in seen:
            continue
        seen.add(path)
        s, v = load(path, chrom, lo, hi)
        order = np.argsort(-v)[:20]
        in_rc = int(np.sum((s[order] >= rc_lo) & (s[order] <= rc_hi)))
        print("  %-22s %4.1f%%  %4.1f%%  %.4f   %d/20" % (
            os.path.basename(path), 100 * np.mean(v > p95),
            100 * np.mean(v > p99), v.max(), in_rc))
    print("  (chance is 5.0%/1.0%; the null row is the by-definition reference)")


if __name__ == "__main__":
    main()
