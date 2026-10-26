#!/usr/bin/env python3
# regional_null_matched_v2.py
# Authors: Katharina E. Hayer (Nina) with Claude (Opus 5.5), co-created
# Version: 1.0   Date: 2026-10-25
# Purpose: B10 option (a) - contrast-matched background null without overlapping draws.
#   Problem (A18): dev/regional_null_calibration.py sample_background_loci() draws random positions with
#   replacement and no overlap check; in the V1P 7-track run (RC 4C viewpoint in the panel) only windows near
#   Tcrb match the scored-bin count, so 327 draws collapsed into 31 overlapping clusters on chr6.
#   Fix: tile every autosome with NON-OVERLAPPING windows of the locus span (random phase per chromosome),
#   keep tiles whose scored-bin count is within +/- tol of the real locus, and sample up to n of them
#   without replacement. All statistics are computed by the ORIGINAL script (unchanged file on disk);
#   this wrapper only swaps in the new sampler at run time.
#
# Two modes:
#   diagnose  standalone, read-only: how many non-overlapping matched loci exist at each tolerance, where
#             they are, and where the real locus sits in the genome-wide distribution of scored-bin counts.
#   run       runs dev/regional_null_calibration.py with the tiled sampler (pass its normal arguments
#             after --). Write to NEW output paths; production files stay untouched.
#
# Usage (from workflow/; TOP=$(git rev-parse --show-toplevel)):
#   python $TOP/dev/regional_null_matched_v2.py diagnose \
#       --stats results_v1p7_prior1_l1/pvalue/diff_DN_vs_V1P.stats.tsv \
#       --chrom-sizes resources/mm10.chrom.sizes --locus chr6:40790000-41690000
#   python $TOP/dev/regional_null_matched_v2.py run --base $TOP/dev/regional_null_calibration.py --tol 0.25 -- \
#       --stats ... --regions ... --locus ... --chrom-sizes ... --mode background --n-random 2000 \
#       --p-thresh 0.05 --seed 42 --verify 5 --repo $TOP --out NEW.tsv --summary-out NEW_summary.tsv
# ASCII-only source.

import argparse
import gzip
import re
import sys

import numpy as np

AUTOSOME_SKIP = ("chrM", "chrX", "chrY")


# ---------------------------------------------------------------------------------------------- diagnose
def read_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    return sizes


def read_tested(stats_path):
    op = gzip.open if stats_path.endswith(".gz") else open
    per = {}
    with op(stats_path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci, si = head.index("chrom"), head.index("start")
        for line in fh:
            f = line.split("\t", max(ci, si) + 1)
            per.setdefault(f[ci], []).append(int(f[si]))
    return {c: np.sort(np.array(v, dtype=np.int64)) for c, v in per.items()}


def diagnose(a):
    lc, rest = a.locus.split(":")
    llo, lhi = [int(x) for x in rest.replace(",", "").split("-")]
    span = lhi - llo
    sizes = read_sizes(a.chrom_sizes)
    tested = read_tested(a.stats)
    s = tested.get(lc, np.array([], dtype=np.int64))
    ref = int(np.searchsorted(s, lhi) - np.searchsorted(s, llo))
    rng = np.random.default_rng(a.seed)
    tiles = []  # (chrom, lo, n_tested)
    for c in sizes:
        if c not in tested or "_" in c or c in AUTOSOME_SKIP or sizes[c] <= span:
            continue
        phase = int(rng.integers(0, span))
        st = tested[c]
        for lo in range(phase, sizes[c] - span + 1, span):
            hi = lo + span
            if c == lc and not (hi <= llo or lo >= lhi):
                continue  # never use the real locus itself
            tiles.append((c, lo, int(np.searchsorted(st, hi) - np.searchsorted(st, lo))))
    counts = np.array([t[2] for t in tiles])
    print("locus %s  span %d bp  scored bins in locus (ref) = %d" % (a.locus, span, ref))
    print("non-overlapping autosomal tiles: %d; scored bins per tile: median %d, 90th %d, 99th %d, max %d"
          % (len(tiles), np.median(counts), np.percentile(counts, 90), np.percentile(counts, 99), counts.max()))
    print("real locus is at the %.2f percentile of tiles" % (100.0 * (counts < ref).mean()))
    for tol in [float(x) for x in a.tol_grid.split(",")]:
        keep = [t for t in tiles if abs(t[2] - ref) <= tol * ref]
        chroms = {}
        for t in keep:
            chroms[t[0]] = chroms.get(t[0], 0) + 1
        top = ", ".join("%s:%d" % kv for kv in sorted(chroms.items(), key=lambda kv: -kv[1])[:6])
        print("  tol %.2f: %4d matched non-overlapping loci on %2d chromosomes  (%s)"
              % (tol, len(keep), len(chroms), top))
    if a.at_least:
        keep = [t for t in tiles if t[2] >= ref]
        print("  at least ref scored bins: %d loci" % len(keep))


# ---------------------------------------------------------------------------------------------- run
NEW_SAMPLER = '''
def sample_background_loci(idx, sizes, span, n_want, ref_bins, rng, exclude, tol=0.25):
    """v2 (regional_null_matched_v2.py): non-overlapping tiles of the locus span, random phase per
    chromosome, matched on scored-bin count within +/- _V2_TOL; up to n_want sampled without replacement."""
    tol = _V2_TOL
    chroms = [c for c in sizes if c in idx and "_" not in c and c not in ("chrM", "chrX", "chrY")]
    cand = []
    for c in chroms:
        if sizes[c] <= span:
            continue
        phase = int(rng.integers(0, span))
        for lo in range(phase, sizes[c] - span + 1, span):
            hi = lo + span
            if overlaps(c, lo, hi, exclude):
                continue
            sl = slice_locus(idx, c, lo, hi)
            if sl is None:
                continue
            nb = sl[0].size
            if abs(nb - ref_bins) > tol * ref_bins:
                continue
            cand.append((c, lo, hi, sl))
    per = {}
    for t in cand:
        per[t[0]] = per.get(t[0], 0) + 1
    print("  [v2] %d non-overlapping matched loci available (tol %.2f) on %d chromosomes: %s"
          % (len(cand), tol, len(per), ", ".join("%s:%d" % kv for kv in sorted(per.items(), key=lambda kv: -kv[1])[:8])))
    if len(cand) > n_want:
        pick = sorted(rng.choice(len(cand), size=n_want, replace=False).tolist())
        cand = [cand[i] for i in pick]
    return cand

'''


def run(a, passthrough):
    src = open(a.base).read()
    pat = re.compile(r"\ndef sample_background_loci\(.*?(?=\ndef )", re.S)
    if len(pat.findall(src)) != 1:
        sys.exit("could not find exactly one sample_background_loci() in %s" % a.base)
    if "def overlaps(" not in src or "def slice_locus(" not in src:
        sys.exit("base script lacks overlaps()/slice_locus(); refusing to patch")
    src = pat.sub(lambda m: "\n" + NEW_SAMPLER, src)
    sys.argv = [a.base] + passthrough
    g = {"__name__": "__main__", "__file__": a.base, "_V2_TOL": a.tol}
    print("[v2] running %s with the tiled non-overlapping sampler, tol %.2f" % (a.base, a.tol))
    exec(compile(src, a.base + " [v2 sampler]", "exec"), g)


def main():
    argv = sys.argv[1:]
    passthrough = []
    if "--" in argv:
        i = argv.index("--")
        argv, passthrough = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser(description="v2 contrast-matched null (non-overlapping tiles)")
    sub = ap.add_subparsers(dest="mode", required=True)
    d = sub.add_parser("diagnose")
    d.add_argument("--stats", required=True)
    d.add_argument("--chrom-sizes", required=True)
    d.add_argument("--locus", required=True)
    d.add_argument("--tol-grid", default="0.25,0.5,1.0")
    d.add_argument("--at-least", action="store_true", help="also count tiles with >= ref scored bins")
    d.add_argument("--seed", type=int, default=42)
    r = sub.add_parser("run")
    r.add_argument("--base", required=True, help="path to dev/regional_null_calibration.py")
    r.add_argument("--tol", type=float, default=0.25)
    a = ap.parse_args(argv)
    if a.mode == "diagnose":
        diagnose(a)
    else:
        run(a, passthrough)


if __name__ == "__main__":
    main()
