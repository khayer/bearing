#!/usr/bin/env python3
"""
build_superbin_grid.py -- coarsen the BEARING adaptive segmentation into a Hi-C
super-bin grid by merging WHOLE adaptive bins until each super-bin clears a
minimum-width floor.

Because every super-bin is a union of whole adaptive bins, the BEARING grid nests
exactly inside it: each BEARING bin maps to exactly one super-bin, and each
super-bin is an exact set of BEARING bins. This is what makes the 1D-to-3D
overlay lossless. Gap-aware: never merges across a gap in the segmentation
(dropped / dead / blacklisted space).

Within each contiguous island the bins are partitioned into
  n = max(1, floor(island_width / min_width))
balanced groups by cumulative width, so groups are ~equal width and >= min_width,
with no runt remainder. An island narrower than min_width becomes a single
sub-floor super-bin (unavoidable; these are isolated, contact-starved islands
that Hi-C low-coverage masking drops anyway).

Inputs:
  --seg        adaptive_segmentation.bed (chrom start end [width cov]); MUST be
               sorted: chromosomes in contiguous blocks, starts increasing.
  --min-width  width floor in bp (default 5000)

Outputs:
  --out-bed    super-bin BED: chrom start end superbin_id width n_bearing_bins
  --out-map    (optional) bearing_to_superbin.tsv:
               chrom start end bearing_idx superbin_id

ASCII only.
"""
import argparse
import sys
import numpy as np


def read_segmentation(path):
    bins = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            bins.append((f[0], int(f[1]), int(f[2])))
    return bins


def check_sorted(bins):
    """Verify chromosomes are in contiguous blocks and starts increase."""
    seen = set()
    prev_chrom = None
    prev_start = -1
    for i, (c, s, e) in enumerate(bins):
        if c != prev_chrom:
            if c in seen:
                sys.exit("[ERROR] chromosome %s is not in one contiguous block "
                         "(seen again at line %d). Sort the BED first." % (c, i + 1))
            seen.add(c)
            prev_chrom = c
            prev_start = -1
        if s < prev_start:
            sys.exit("[ERROR] starts not increasing within %s at line %d "
                     "(%d after %d). Sort the BED first." % (c, i + 1, s, prev_start))
        prev_start = s


def islands(bins):
    """Group consecutive bins with no gap (end[i]==start[i+1], same chrom).

    Yields (global_start_index, list_of_bins).
    """
    out = []
    cur = []
    cur_g0 = 0
    for i, b in enumerate(bins):
        if not cur:
            cur = [b]
            cur_g0 = i
            continue
        pc, ps, pe = cur[-1]
        if b[0] == pc and b[1] == pe:
            cur.append(b)
        else:
            out.append((cur_g0, cur))
            cur = [b]
            cur_g0 = i
    if cur:
        out.append((cur_g0, cur))
    return out


def partition_by_width(island, n):
    """Split contiguous (chrom,start,end) bins into n groups of WHOLE bins,
    balanced by cumulative width. Returns list of (a, b) half-open index pairs."""
    if n <= 1:
        return [(0, len(island))]
    w = np.array([e - s for _, s, e in island], dtype=np.float64)
    cum = np.concatenate([[0.0], np.cumsum(w)])
    targets = np.linspace(0.0, cum[-1], n + 1)
    bnd = np.searchsorted(cum, targets, side="left")
    bnd[0] = 0
    bnd[-1] = len(island)
    for k in range(1, len(bnd)):
        if bnd[k] <= bnd[k - 1]:
            bnd[k] = min(bnd[k - 1] + 1, len(island))
    groups = []
    for k in range(len(bnd) - 1):
        a, b = int(bnd[k]), int(bnd[k + 1])
        if b > a:
            groups.append((a, b))
    return groups


def enforce_floor(island, groups, floor):
    """Merge any sub-floor group into its smaller-width neighbor (within the
    island only), so every group clears the floor. The balanced partition can
    orphan a narrow remainder next to a very wide adaptive bin; this repairs it.
    If the whole island is narrower than the floor it stays a single group."""
    def gw(g):
        a, b = g
        return island[b - 1][2] - island[a][1]
    changed = True
    while changed and len(groups) > 1:
        changed = False
        for i, g in enumerate(groups):
            if gw(g) < floor:
                if i == 0:
                    j = 1
                elif i == len(groups) - 1:
                    j = i - 1
                else:
                    j = i - 1 if gw(groups[i - 1]) <= gw(groups[i + 1]) else i + 1
                lo, hi = min(i, j), max(i, j)
                groups = groups[:lo] + [(groups[lo][0], groups[hi][1])] + groups[hi + 1:]
                changed = True
                break
    return groups


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seg", required=True)
    ap.add_argument("--min-width", type=int, default=5000)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--out-map", default=None)
    ap.add_argument("--report-region", default=None,
                    help="CHR:START-END to summarize super-bin widths in, e.g. "
                         "chr6:40400000-42400000")
    args = ap.parse_args()

    bins = read_segmentation(args.seg)
    if not bins:
        sys.exit("[ERROR] empty segmentation")
    check_sorted(bins)

    superbins = []   # (chrom, start, end, width, n_bearing)
    mapping = []     # (chrom, start, end, bearing_idx, superbin_id)
    sb_id = 0
    sub_floor = 0

    for g0, island in islands(bins):
        W = sum(e - s for _, s, e in island)
        n = max(1, int(W // args.min_width))
        groups = enforce_floor(island, partition_by_width(island, n), args.min_width)
        for (a, b) in groups:
            chrom = island[a][0]
            start = island[a][1]
            end = island[b - 1][2]
            width = end - start
            superbins.append((chrom, start, end, width, b - a))
            if width < args.min_width:
                sub_floor += 1
            if args.out_map:
                for j in range(a, b):
                    c, s, e = island[j]
                    mapping.append((c, s, e, g0 + j, sb_id))
            sb_id += 1

    with open(args.out_bed, "w") as out:
        for sid, (c, s, e, w, nb) in enumerate(superbins):
            out.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (c, s, e, sid, w, nb))
    if args.out_map:
        with open(args.out_map, "w") as out:
            out.write("chrom\tstart\tend\tbearing_idx\tsuperbin_id\n")
            for c, s, e, bi, sid in mapping:
                out.write("%s\t%d\t%d\t%d\t%d\n" % (c, s, e, bi, sid))

    widths = np.array([w for _, _, _, w, _ in superbins])
    n_bear = np.array([nb for _, _, _, _, nb in superbins])
    sys.stderr.write(
        "[done] %d super-bins from %d BEARING bins -> %s\n"
        "       width median=%d mean=%d min=%d max=%d ; sub-floor (isolated)=%d\n"
        "       BEARING bins per super-bin: median=%d max=%d\n"
        % (len(superbins), len(bins), args.out_bed,
           int(np.median(widths)), int(widths.mean()), int(widths.min()),
           int(widths.max()), sub_floor,
           int(np.median(n_bear)), int(n_bear.max())))

    if args.report_region:
        chrom, rng = args.report_region.split(":")
        a, b = (int(x.replace(",", "")) for x in rng.split("-"))
        sel = [w for (c, s, e, w, _) in superbins if c == chrom and s >= a and e <= b]
        if sel:
            sel = np.array(sel)
            sys.stderr.write(
                "[region %s] %d super-bins; width median=%d min=%d max=%d ; "
                "under floor=%d\n"
                % (args.report_region, len(sel), int(np.median(sel)),
                   int(sel.min()), int(sel.max()), int((sel < args.min_width).sum())))
        else:
            sys.stderr.write("[region %s] no super-bins\n" % args.report_region)


if __name__ == "__main__":
    main()
