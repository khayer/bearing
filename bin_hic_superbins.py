#!/usr/bin/env python3
"""
bin_hic_superbins.py -- bin Hi-C valid-pair BAMs onto the BEARING super-bin grid
at one or more locus windows (the pilot), producing one raw-count cooler per
condition for the integration overlay.

What it does, per condition:
  - downsamples every library to the lowest pooled condition (the floor), with a
    reproducible per-pair Bernoulli decision keyed on read name (rep-independent,
    deterministic across machines via CRC32);
  - for each --region window, fetches only that window from the indexed BAM,
    keeps within-window CIS read pairs (both ends on the same chromosome and
    inside the window), assigns each end to a super-bin by searchsorted with a
    gap-drop, and pools the two replicates;
  - writes results/<cond>.superbin5kb.pilot.cool holding RAW counts plus a per-bin
    'width' column. Width-density normalization (count / (w_i * w_j)) and the A/B
    ratio are deliberately left to plot time so the stored matrix stays raw and
    inspectable. Use raw_*_valid_pairs.bam, NOT raw_corrected_KR_* (already
    balanced on the old fixed grid).

Inputs:
  --superbins  superbins_5kb.bed (chrom start end superbin_id width n_bearing)
  --counts     pair_counts_by_lib.tsv: cond <tab> bam <tab> pair_count
               (the file you already built; one row per library)
  --bam-dir    directory holding the bam basenames in --counts (ignored for any
               row whose bam path is already absolute)
  --regions    comma-separated CHR:START-END windows, e.g.
               chr6:40400000-42400000,chr12:113200000-116000000
  --conditions optional comma-separated subset (default: all in --counts)
  --out-dir    output directory for the .cool files and the QC table
  --seed       integer folded into the downsample hash (default 1)

Requires: pysam, cooler, numpy, pandas. BAMs must be coordinate-sorted+indexed
(.bai). ASCII only.
"""
import argparse
import os
import sys
import zlib
import numpy as np
import pandas as pd


def parse_regions(s):
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        chrom, rng = tok.split(":")
        a, b = (int(x.replace(",", "")) for x in rng.split("-"))
        if b <= a:
            sys.exit("[ERROR] region end <= start: %s" % tok)
        out.append((chrom, a, b))
    return out


def load_superbins(path):
    """chrom -> dict(starts, ends, gid, width) as numpy arrays, sorted by start."""
    by = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            c = f[0]
            by.setdefault(c, []).append((int(f[1]), int(f[2]), int(f[3]), int(f[4])))
    out = {}
    for c, rows in by.items():
        rows.sort()
        out[c] = {
            "starts": np.array([r[0] for r in rows], dtype=np.int64),
            "ends":   np.array([r[1] for r in rows], dtype=np.int64),
            "gid":    np.array([r[2] for r in rows], dtype=np.int64),
            "width":  np.array([r[3] for r in rows], dtype=np.int64),
        }
    return out


def region_slice(sb_chrom, wstart, wend):
    """Sub-arrays of a chromosome's super-bins that OVERLAP [wstart, wend)."""
    s, e = sb_chrom["starts"], sb_chrom["ends"]
    m = (e > wstart) & (s < wend)
    return {k: v[m] for k, v in sb_chrom.items()}


def assign(pos, starts, ends):
    """Local index of the super-bin containing pos, or -1 (gap / out of range)."""
    i = int(np.searchsorted(starts, pos, side="right")) - 1
    if i < 0 or pos >= ends[i]:
        return -1
    return i


def try_pair(r, c, a, b, starts, ends, off, keep_p, seedstr, counts):
    """Process a primary read1 record; returns 1 if a contact was tallied.
    Keeps only within-window CIS pairs (both ends on c, inside [a,b)), applies
    the reproducible per-pair downsample, gap-drops ends, and accumulates."""
    if r.next_reference_id < 0 or r.next_reference_name != c:
        return 0
    p1 = r.reference_start
    if p1 < a or p1 >= b:
        return 0
    p2 = r.next_reference_start
    if p2 < a or p2 >= b:
        return 0
    if zlib.crc32((seedstr + r.query_name).encode()) / 4294967296.0 >= keep_p:
        return 0
    b1 = assign(p1, starts, ends)
    if b1 < 0:
        return 0
    b2 = assign(p2, starts, ends)
    if b2 < 0:
        return 0
    g1, g2 = off + b1, off + b2
    if g1 > g2:
        g1, g2 = g2, g1
    counts[(g1, g2)] = counts.get((g1, g2), 0) + 1
    return 1


def load_counts(path):
    """Return list of (cond, bam_basename_or_path, count)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                continue
            rows.append((f[0], f[1], int(f[2])))
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--superbins", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--bam-dir", default=".")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mode", choices=["auto","fetch","scan"], default="auto",
                    help="auto: fetch if BAM is coordinate-sorted+indexed, else "
                         "scan (one streaming pass, no index needed). HiC-Pro "
                         "valid_pairs are usually name-sorted -> scan.")
    args = ap.parse_args()

    import pysam
    import cooler

    regions = parse_regions(args.regions)
    sb = load_superbins(args.superbins)
    for c, a, b in regions:
        if c not in sb:
            sys.exit("[ERROR] no super-bins on %s (check chrom naming in %s)"
                     % (c, args.superbins))

    # global pilot bins table: region super-bins concatenated in region order,
    # each region's bins sorted by start. local id = row index. We also keep, per
    # region, the local-id offset and the sliced start/end arrays for assignment.
    bin_chrom, bin_start, bin_end, bin_w, bin_gid = [], [], [], [], []
    region_arrays = []  # (chrom, wstart, wend, starts, ends, local_offset, n)
    for (c, a, b) in regions:
        sl = region_slice(sb[c], a, b)
        n = len(sl["starts"])
        if n == 0:
            sys.exit("[ERROR] region %s:%d-%d covers no super-bins" % (c, a, b))
        off = len(bin_chrom)
        bin_chrom += [c] * n
        bin_start += sl["starts"].tolist()
        bin_end += sl["ends"].tolist()
        bin_w += sl["width"].tolist()
        bin_gid += sl["gid"].tolist()
        region_arrays.append((c, a, b, sl["starts"], sl["ends"], off, n))
    bins_df = pd.DataFrame({
        "chrom": bin_chrom, "start": bin_start, "end": bin_end,
        "width": bin_w, "superbin_id": bin_gid,
    })
    sys.stderr.write("[grid] %d pilot super-bins across %d region(s)\n"
                     % (len(bins_df), len(regions)))

    # condition pooled totals -> floor -> per-condition keep fraction
    rows = load_counts(args.counts)
    pooled = {}
    libs = {}
    for cond, bam, cnt in rows:
        pooled[cond] = pooled.get(cond, 0) + cnt
        libs.setdefault(cond, []).append((bam, cnt))
    if args.conditions:
        want = set(x.strip() for x in args.conditions.split(","))
        pooled = {k: v for k, v in pooled.items() if k in want}
        libs = {k: v for k, v in libs.items() if k in want}
    if not pooled:
        sys.exit("[ERROR] no conditions selected")
    floor = min(pooled.values())
    floor_cond = min(pooled, key=pooled.get)
    sys.stderr.write("[floor] %s = %d pooled pairs (downsample target)\n"
                     % (floor_cond, floor))

    os.makedirs(args.out_dir, exist_ok=True)
    qc_path = os.path.join(args.out_dir, "pilot_qc.tsv")
    qc = open(qc_path, "w")
    qc.write("condition\tregion\tkept_pairs\tnnz_pixels\ttotal_contacts\t"
             "n_bins\tzero_marginal_bins\tmedian_marginal\n")

    seedstr = "%d:" % args.seed
    thr = {}  # per-condition CRC threshold (keep if crc/2^32 < frac)
    for cond in pooled:
        thr[cond] = floor / pooled[cond]

    for cond in sorted(pooled):
        keep_p = thr[cond]
        counts = {}                  # (b1,b2) local -> count (pooled over reps)
        kept_per_region = {i: 0 for i in range(len(regions))}
        for bam_name, _cnt in libs[cond]:
            path = bam_name if os.path.isabs(bam_name) else os.path.join(args.bam_dir, bam_name)
            if not os.path.exists(path):
                sys.exit("[ERROR] missing BAM: %s" % path)
            bam = pysam.AlignmentFile(path, "rb")
            refset = set(bam.references)
            for (c, a, b, starts, ends, off, n) in region_arrays:
                if c not in refset:
                    sys.exit("[ERROR] %s not a reference in %s (chrom naming?)" % (c, path))
            so = bam.header.to_dict().get("HD", {}).get("SO", "")
            use_fetch = (args.mode == "fetch") or (
                args.mode == "auto" and so == "coordinate" and bam.has_index())
            if args.mode == "fetch" and not bam.has_index():
                sys.exit("[ERROR] --mode fetch needs an index for %s" % path)
            sys.stderr.write("  %s : %s mode (SO=%s, indexed=%s)\n"
                             % (os.path.basename(path),
                                "fetch" if use_fetch else "scan", so or "?",
                                bam.has_index()))
            if use_fetch:
                for ri, (c, a, b, starts, ends, off, n) in enumerate(region_arrays):
                    for r in bam.fetch(c, a, b):
                        if (r.flag & 0x900) or not (r.flag & 0x40):
                            continue
                        kept_per_region[ri] += try_pair(
                            r, c, a, b, starts, ends, off, keep_p, seedstr, counts)
            else:
                by_chrom = {}
                for ri, (c, a, b, starts, ends, off, n) in enumerate(region_arrays):
                    by_chrom.setdefault(c, []).append((ri, c, a, b, starts, ends, off))
                for r in bam.fetch(until_eof=True):
                    if (r.flag & 0x900) or not (r.flag & 0x40):
                        continue
                    lst = by_chrom.get(r.reference_name)
                    if not lst:
                        continue
                    for (ri, c, a, b, starts, ends, off) in lst:
                        kept_per_region[ri] += try_pair(
                            r, c, a, b, starts, ends, off, keep_p, seedstr, counts)
            bam.close()

        # per-region QC + write cooler
        marg = np.zeros(len(bins_df), dtype=np.int64)
        for (g1, g2), c0 in counts.items():
            marg[g1] += c0
            if g2 != g1:
                marg[g2] += c0
        for ri, (c, a, b, starts, ends, off, n) in enumerate(region_arrays):
            seg = marg[off:off + n]
            nnz = sum(1 for (g1, g2) in counts
                      if off <= g1 < off + n and off <= g2 < off + n)
            tot = sum(c0 for (g1, g2), c0 in counts.items()
                      if off <= g1 < off + n and off <= g2 < off + n)
            zero = int((seg == 0).sum())
            med = int(np.median(seg)) if n else 0
            qc.write("%s\t%s:%d-%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
                     % (cond, c, a, b, kept_per_region[ri], nnz, tot, n, zero, med))

        if not counts:
            sys.stderr.write("[WARN] %s: zero kept pairs in all regions; "
                             "no cooler written\n" % cond)
            continue
        px = sorted(counts.items())
        pixels = pd.DataFrame({
            "bin1_id": [k[0] for k, _ in px],
            "bin2_id": [k[1] for k, _ in px],
            "count":   [v for _, v in px],
        })
        cool_path = os.path.join(args.out_dir, "%s.superbin5kb.pilot.cool" % cond)
        cooler.create_cooler(cool_path, bins_df, pixels,
                             dtypes={"count": "int32"},
                             symmetric_upper=True, ordered=True,
                             assembly="mm10")
        sys.stderr.write("[cool] %s -> %s (%d pixels, %d contacts)\n"
                         % (cond, cool_path, len(pixels), int(pixels["count"].sum())))

    qc.close()
    sys.stderr.write("[qc] %s\n" % qc_path)


if __name__ == "__main__":
    main()
