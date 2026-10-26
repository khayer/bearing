#!/usr/bin/env python3
# regional_4c_cis_null.py
# Authors: Katharina E. Hayer (Nina) with Claude (Opus 5.5), co-created
# Version: 1.0   Date: 2026-10-25
# Purpose: B10 split design, part 2 - regional test for a VIEWPOINT track (4C / virtual 4C O/E) against a
#   cis-reach null. A viewpoint track only has signal in cis around its viewpoint, so genome-wide background
#   loci cannot serve as its null. Instead, each test region is compared with every non-overlapping window of
#   the same width inside the 4C reach on the same chromosome (test regions themselves excluded).
#   Statistic per window = sum over its bins of the signed differential d of the chosen category
#   (from the diff qcat; (+) = first-named condition higher, A8). Upper tail tests "first condition higher",
#   lower tail tests "second condition higher"; both are reported with BH over the regions per tail.
#   O/E tracks already remove most distance decay; a distance-band sensitivity (windows whose distance to the
#   viewpoint is within +/- --band of the region's) is reported alongside.
# Usage (from workflow/):
#   python $(git rev-parse --show-toplevel)/dev/regional_4c_cis_null.py \
#     --diff-qcat results_v1p7_prior1_l1/compare/diff_DN_vs_V1P.qcat.bgz --cat 7 \
#     --regions $(git rev-parse --show-toplevel)/annotations/tcrb_regions_v5.bed \
#     --reach chr6:35000000-48000000 --viewpoint 41505538 \
#     --out results_v1p7_prior1_l1/table_sources/regional_4c_cis_null_tcrb_DN_vs_V1P.tsv
# Needs tabix-indexed qcat (pysam or the tabix binary). ASCII-only source.

import argparse
import json
import re
import subprocess
import sys

QCAT_RE = re.compile(r"qcat:(\[.*\])")


def parse_window(s):
    c, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return c, int(a), int(b)


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


def read_regions(path, chrom):
    regs = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom:
                regs.append((f[3] if len(f) > 3 else "%s:%s-%s" % tuple(f[:3]), int(f[1]), int(f[2])))
    return regs


def bh(p):
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i])
    q, prev = [0.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * m / rank)
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diff-qcat", required=True)
    ap.add_argument("--cat", type=int, required=True, help="category number of the viewpoint track")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--reach", required=True, help="chrom:start-end of the 4C reach used as the null")
    ap.add_argument("--viewpoint", type=int, required=True, help="viewpoint coordinate (bp)")
    ap.add_argument("--band", type=float, default=0.5,
                    help="distance-band sensitivity: keep windows with |dist - region dist| <= band * region dist")
    ap.add_argument("--bin", type=int, default=200)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    chrom, rlo, rhi = parse_window(a.reach)
    d = {}
    for line in fetch_lines(a.diff_qcat, chrom, rlo, rhi):
        f = line.split("\t")
        m = QCAT_RE.search(f[3])
        if not m:
            continue
        for v, c in json.loads(m.group(1)):
            if int(c) == a.cat:
                d[int(f[1])] = float(v)
                break
    if not d:
        sys.exit("no bins of category %d in %s" % (a.cat, a.reach))
    nz = sorted(b for b, v in d.items() if v != 0.0)
    print("category %d: %d bins read in %s; nonzero d from %s to %s (%d bins)"
          % (a.cat, len(d), a.reach, nz[0] if nz else "-", nz[-1] if nz else "-", len(nz)))

    regs = read_regions(a.regions, chrom)
    regs = [r for r in regs if r[1] >= rlo and r[2] <= rhi]
    if not regs:
        sys.exit("no regions inside the reach")

    def wsum(lo, hi):
        return sum(d.get(b, 0.0) for b in range(lo - lo % a.bin, hi, a.bin) if lo <= b < hi)

    def mid_dist(lo, hi):
        return abs((lo + hi) // 2 - a.viewpoint)

    rows = []
    for name, lo, hi in regs:
        w = hi - lo
        obs = wsum(lo, hi)
        dist = mid_dist(lo, hi)
        null, nullband = [], []
        start = rlo + (lo - rlo) % w  # tile aligned to the region's own phase
        for s in range(start, rhi - w + 1, w):
            e = s + w
            if any(not (e <= rlo2 or s >= rhi2) for _n, rlo2, rhi2 in regs):
                continue
            x = wsum(s, e)
            null.append(x)
            if abs(mid_dist(s, e) - dist) <= a.band * max(dist, 1):
                nullband.append(x)
        n, nb = len(null), len(nullband)
        up = (1.0 + sum(1 for x in null if x >= obs)) / (1.0 + n)
        dn = (1.0 + sum(1 for x in null if x <= obs)) / (1.0 + n)
        upb = (1.0 + sum(1 for x in nullband if x >= obs)) / (1.0 + nb) if nb else float("nan")
        dnb = (1.0 + sum(1 for x in nullband if x <= obs)) / (1.0 + nb) if nb else float("nan")
        mean = sum(null) / n if n else float("nan")
        rows.append([name, lo, hi, w // a.bin, dist, obs, n, mean, max(null) if n else float("nan"),
                     min(null) if n else float("nan"), up, dn, nb, upb, dnb])
    qu = bh([r[10] for r in rows])
    qd = bh([r[11] for r in rows])

    with open(a.out, "w") as fh:
        fh.write("region\tstart\tend\tbins\tdist_to_viewpoint\tobs_sum_d\tn_null\tnull_mean\tnull_max\t"
                 "null_min\temp_p_first_higher\temp_p_second_higher\tq_first_higher\tq_second_higher\t"
                 "n_null_band\temp_p_first_higher_band\temp_p_second_higher_band\n")
        for r, a1, a2 in zip(rows, qu, qd):
            fh.write("%s\t%d\t%d\t%d\t%d\t%.6g\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%.6g\t%.6g\n"
                     % tuple(r[:12] + [a1, a2] + r[12:]))
    print("region\tobs\tn_null\tnull_mean\tnull_max\tp_first_higher\tq\tp_second_higher\tq\tn_band\tp_first_band")
    for r, a1, a2 in zip(rows, qu, qd):
        print("%s\t%.4g\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%d\t%.4g"
              % (r[0], r[5], r[6], r[7], r[8], r[10], a1, r[11], a2, r[12], r[13]))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
