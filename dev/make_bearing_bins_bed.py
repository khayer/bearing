#!/usr/bin/env python3
"""
make_bearing_bins_bed.py -- emit the exact bin grid BEARING scores, as a BED, for
bin-to-bin concordance analyses (e.g. the CTCF edgeR check).

The grid is generated with the SAME function BEARING uses (bins_for_chrom), so
the bins match the qcat coordinates exactly -- no off-by-one, no re-binning, no
peak mapping. Coordinates are 0-based half-open (standard BED): bin 0 is
(0, 200), bin 1 is (200, 400), ... and the last bin of each chromosome may be
shorter than bin-size, exactly as BEARING produces it.

Blacklist filtering (recommended): BEARING zeroes blacklisted bins during
scoring, so they are never significant. Excluding them here makes the analysis
universe match the bins BEARING can actually call. Pass the SAME blacklist used
in the BEARING run. Low-signal bins are NOT dropped here -- let edgeR's
filterByExpr handle near-empty bins, so the universe stays a clean function of
(chrom grid + blacklist) only.

Examples:
  # genome-wide, blacklist-filtered  -> bearing_bins_200bp.bed
  python make_bearing_bins_bed.py \
    --chrom-sizes workflow/resources/mm10.chrom.sizes \
    --blacklist   workflow/resources/mm10-blacklist.v2.bed \
    --out bearing_bins_200bp.bed

  # just the wide Tcrb locus -> bearing_bins_tcrb_wide_200bp.bed
  python make_bearing_bins_bed.py \
    --chrom-sizes workflow/resources/mm10.chrom.sizes \
    --blacklist   workflow/resources/mm10-blacklist.v2.bed \
    --region chr6:40400000-42400000 \
    --out bearing_bins_tcrb_wide_200bp.bed

ASCII only.
"""
import argparse
import sys

from bigwig_to_qcat import (
    bins_for_chrom,
    load_blacklist,
    bins_overlapping_blacklist,
    BIN_SIZE,
)


def parse_region(s):
    """'chr6:40400000-42400000' -> ('chr6', 40400000, 42400000)."""
    if ":" not in s or "-" not in s:
        sys.exit("[ERROR] --region must look like chr6:40400000-42400000")
    chrom, rng = s.split(":", 1)
    a, b = rng.replace(",", "").split("-", 1)
    start, end = int(a), int(b)
    if end <= start:
        sys.exit("[ERROR] --region end must be greater than start")
    return chrom, start, end


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom-sizes", required=True,
                    help="two-column chrom<TAB>size file (same one BEARING used)")
    ap.add_argument("--bin-size", type=int, default=BIN_SIZE,
                    help="bin width in bp (default %d, the BEARING grid)" % BIN_SIZE)
    ap.add_argument("--blacklist", default=None,
                    help="BED of blacklist regions; bins overlapping these are "
                         "excluded (use the SAME blacklist as the BEARING run)")
    ap.add_argument("--region", default=None,
                    help="restrict to CHR:START-END (e.g. chr6:40400000-42400000); "
                         "bins are still on the genome-aligned grid")
    ap.add_argument("--chroms", default=None,
                    help="comma-separated chromosome subset (ignored if --region)")
    ap.add_argument("--out", required=True, help="output BED path")
    args = ap.parse_args()

    sizes = {}
    with open(args.chrom_sizes) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    if not sizes:
        sys.exit("[ERROR] no chromosomes parsed from --chrom-sizes")

    blacklist = load_blacklist(args.blacklist) if args.blacklist else None

    region = parse_region(args.region) if args.region else None
    if region:
        chroms = [region[0]]
    elif args.chroms:
        chroms = args.chroms.split(",")
    else:
        chroms = list(sizes.keys())
    missing = [c for c in chroms if c not in sizes]
    if missing:
        sys.exit("[ERROR] chromosome(s) not in --chrom-sizes: %s" % ",".join(missing))

    n_total = 0
    n_bl = 0
    with open(args.out, "w") as out:
        for chrom in chroms:
            bins = bins_for_chrom(sizes[chrom], args.bin_size)
            if region is not None:
                _, rs, re_ = region
                # keep bins whose start falls inside the region (genome-aligned)
                bins = [(s, e) for (s, e) in bins if rs <= s < re_]
            if blacklist is not None and bins:
                mask = bins_overlapping_blacklist(bins, blacklist, chrom)
                n_bl += int(mask.sum())
                bins = [b for b, m in zip(bins, mask) if not m]
            for (s, e) in bins:
                out.write("%s\t%d\t%d\n" % (chrom, s, e))
            n_total += len(bins)

    msg = "[done] %d bins (%d bp) -> %s" % (n_total, args.bin_size, args.out)
    if blacklist is not None:
        msg += " (excluded %d blacklisted)" % n_bl
    if region is not None:
        msg += " [region %s:%d-%d]" % region
    sys.stderr.write(msg + "\n")


if __name__ == "__main__":
    main()
