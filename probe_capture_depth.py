#!/usr/bin/env python3
"""
probe_capture_depth.py -- how fine can the capture-Hi-C be binned at a locus
before contacts go sparse? Answers the gating question for feature-anchored
binning: it is only worth re-testing adaptive vs uniform binning at a resolution
where bins are actually populated.

Scans each valid-pairs BAM ONCE, keeps cis pairs with both ends inside the
window, pools replicates per condition, then for each requested resolution
reports per-bin occupancy. The finest resolution where bins stay populated
(low zero-marginal fraction, non-trivial median pixel) is the floor that
feature-binning could exploit.

INPUT
  --manifest TSV (header) columns: cond  bam   (one row per replicate file)
  --bam-dir  dir holding the bam basenames (default .)
  --region   chrom:start-end (the capture/locus window)
  --resolutions comma list in bp (default 250,500,1000,2000,3000,5000,10000)
  --out      output TSV

OUTPUT per (cond, resolution): n_pairs_in_window, n_bins, median_marginal,
  frac_zero_marginal_bins, median_nonzero_pixel, max_pixel, mean_pairs_per_bin

ASCII only. Requires pysam, numpy.
"""
import argparse
import os
import sys

import numpy as np


def scan_positions(bam_paths, chrom, a, b):
    """Return (p1, p2) int arrays for cis pairs with both ends in [a,b)."""
    import pysam
    P1, P2 = [], []
    for path in bam_paths:
        bam = pysam.AlignmentFile(path, "rb")
        so = bam.header.to_dict().get("HD", {}).get("SO", "")
        use_fetch = (so == "coordinate") and bam.has_index()
        itr = bam.fetch(chrom, a, b) if use_fetch else bam.fetch(until_eof=True)
        sys.stderr.write("  scan %s (%s)\n" % (os.path.basename(path),
                                               "fetch" if use_fetch else "scan"))
        for r in itr:
            if (r.flag & 0xD00) or not (r.flag & 0x40):
                continue
            if not use_fetch and r.reference_name != chrom:
                continue
            if r.next_reference_id < 0 or r.next_reference_name != chrom:
                continue
            p1, p2 = r.reference_start, r.next_reference_start
            if p1 < a or p1 >= b or p2 < a or p2 >= b:
                continue
            P1.append(p1); P2.append(p2)
        bam.close()
    return np.asarray(P1, dtype=np.int64), np.asarray(P2, dtype=np.int64)


def occupancy(p1, p2, a, b, res):
    n = int(np.ceil((b - a) / res))
    i = ((p1 - a) // res).astype(np.int64)
    j = ((p2 - a) // res).astype(np.int64)
    lo = np.minimum(i, j); hi = np.maximum(i, j)
    # pixel counts on upper triangle via flattened index
    flat = lo * n + hi
    uniq, cnt = np.unique(flat, return_counts=True)
    # per-bin marginal (each pixel contributes to both bins)
    marg = np.zeros(n, dtype=np.int64)
    np.add.at(marg, lo, 1)
    np.add.at(marg, hi, 1)
    nz_pix = cnt  # counts of occupied pixels
    return {
        "n_bins": n,
        "median_marginal": float(np.median(marg)),
        "frac_zero_marginal": float((marg == 0).mean()),
        "median_nonzero_pixel": float(np.median(nz_pix)) if len(nz_pix) else 0.0,
        "max_pixel": int(nz_pix.max()) if len(nz_pix) else 0,
        "mean_pairs_per_bin": float(len(p1) * 2.0 / n),
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--bam-dir", default=".")
    ap.add_argument("--region", required=True)
    ap.add_argument("--resolutions", default="250,500,1000,2000,3000,5000,10000")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    chrom, rng = args.region.split(":")
    a, b = (int(x.replace(",", "")) for x in rng.split("-"))
    res_list = [int(x) for x in args.resolutions.split(",")]

    libs = {}
    with open(args.manifest) as fh:
        hdr = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if hdr is None:
                hdr = [c.lower() for c in f]; continue
            rec = dict(zip(hdr, f))
            libs.setdefault(rec["cond"], []).append(rec["bam"])

    rows = []
    for cond, bams in libs.items():
        paths = [x if os.path.isabs(x) else os.path.join(args.bam_dir, x) for x in bams]
        sys.stderr.write("[cond] %s (%d libs)\n" % (cond, len(paths)))
        p1, p2 = scan_positions(paths, chrom, a, b)
        npairs = len(p1)
        sys.stderr.write("  %d cis pairs in window\n" % npairs)
        for res in res_list:
            if npairs == 0:
                rows.append([cond, res, 0, "NA", "NA", "NA", "NA", "NA", "NA"]); continue
            o = occupancy(p1, p2, a, b, res)
            rows.append([cond, res, npairs, o["n_bins"],
                         "%.0f" % o["median_marginal"],
                         "%.3f" % o["frac_zero_marginal"],
                         "%.0f" % o["median_nonzero_pixel"],
                         o["max_pixel"], "%.1f" % o["mean_pairs_per_bin"]])
            sys.stderr.write("  res=%-6d n_bins=%-5d med_marg=%-6.0f "
                             "zero=%.3f med_pix=%.0f\n" % (
                res, o["n_bins"], o["median_marginal"],
                o["frac_zero_marginal"], o["median_nonzero_pixel"]))

    hdr = ["cond", "resolution_bp", "n_pairs_window", "n_bins",
           "median_marginal", "frac_zero_marginal_bins",
           "median_nonzero_pixel", "max_pixel", "mean_pairs_per_bin"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    sys.stderr.write("[out] %s (%d rows)\n" % (args.out, len(rows)))


if __name__ == "__main__":
    main()
