#!/usr/bin/env python3
"""
crosslocus_fdr_counts.py
========================
Count FDR-significant DN-vs-EbKO bins per survey locus (ALL bins, not just the
top-N reported by the decomposition figure). This is the number that settles
whether the cross-loci are null after multiple-testing correction.

Reads a single diff stats TSV (e.g. diff_DN_vs_EbKO.stats.tsv, which carries a
per-bin `significant_fdr<alpha>` column written by bearing_pvalue.py) plus the
survey region BED, and reports per locus:
  - total bins in the locus
  - bins above a descriptive |BES| threshold (context for the figure)
  - FDR-significant bins (the decision-relevant count)
  - which track dominates the FDR-significant bins (by summed |kl|)

ASCII-only. Transparently handles gzip-compressed input.
"""

import argparse
import csv
import gzip
import sys


def _open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def read_bed(path):
    regions = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            name = p[3] if len(p) > 3 else "%s:%s-%s" % (p[0], p[1], p[2])
            # strip any trailing inline comment on the name field
            name = name.split("#")[0].strip()
            regions.append((p[0], int(p[1]), int(p[2]), name))
    return regions


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff", required=True,
                    help="diff_<A>_vs_<B>.stats.tsv (has significant_fdr* column)")
    ap.add_argument("--regions", required=True, help="survey region BED")
    ap.add_argument("--bes-threshold", type=float, default=0.5,
                    help="descriptive |BES| threshold for context (default 0.5)")
    ap.add_argument("--out", required=True, help="output per-locus TSV")
    args = ap.parse_args()

    regions = read_bed(args.regions)

    # Read the diff stats once.
    fh = _open(args.diff)
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    fdr_col = next((c for c in fields if c.startswith("significant_fdr")), None)
    if fdr_col is None:
        sys.exit("ERROR: no significant_fdr* column in %s" % args.diff)
    kl_cols = [c for c in fields if c.startswith("kl_")]
    rows = list(reader)
    fh.close()

    # Pre-parse rows to typed tuples for speed.
    parsed = []
    for r in rows:
        try:
            c = r["chrom"]; s = int(r["start"]); e = int(r["end"])
            bes = float(r.get("bearing_score", "nan"))
            sig = int(float(r.get(fdr_col, "0") or 0))
        except (ValueError, KeyError):
            continue
        parsed.append((c, s, e, bes, sig, r))

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus", "chrom", "start", "end", "n_bins",
                    "n_above_bes_thresh", "n_fdr_significant",
                    "frac_fdr_significant", "dominant_track_in_sig_bins"])
        for (chrom, start, end, name) in regions:
            n_bins = n_above = n_sig = 0
            kl_abs_sum = {k: 0.0 for k in kl_cols}
            for (c, s, e, bes, sig, r) in parsed:
                if c != chrom or e <= start or s >= end:
                    continue
                n_bins += 1
                if abs(bes) >= args.bes_threshold:
                    n_above += 1
                if sig == 1:
                    n_sig += 1
                    for k in kl_cols:
                        try:
                            kl_abs_sum[k] += abs(float(r.get(k, "0") or 0))
                        except ValueError:
                            pass
            dom = ""
            if n_sig > 0 and kl_cols:
                dom = max(kl_abs_sum, key=kl_abs_sum.get).replace("kl_", "")
            frac = (n_sig / n_bins) if n_bins else 0.0
            w.writerow([name, chrom, start, end, n_bins, n_above, n_sig,
                        "%.4f" % frac, dom])

    # Echo a compact summary to stderr.
    sys.stderr.write("Wrote per-locus FDR counts: %s\n" % args.out)
    with open(args.out) as fh:
        for ln in fh:
            sys.stderr.write("  " + ln)


if __name__ == "__main__":
    main()
