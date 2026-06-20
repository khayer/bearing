#!/usr/bin/env python3
"""Pre-flight validator for the BEARING permutation null set.

Run this BEFORE launching the p-value rebuild. It confirms every expected
permutation qcat exists -- per-sample (perm{N}/{sample}/{sample}_perm{N}.qcat.bgz)
and differential (perm{N}/diff_comparison/diff_{comp}.qcat.bgz) -- and that none
is a size OUTLIER relative to its peers.

DESIGN: a truncated/failed perm write produces a small file. Rather than a fixed
size floor (which a 1 MB-but-broken file can sneak past while its 150 MB siblings
pass), this compares each file to the MEDIAN size of its class (per-sample vs
diff) and flags anything below --min-frac of that median. This is stat-only --
no decompression -- so the whole 100-perm set checks in milliseconds, and it
adapts automatically to KL vs JSD vs whatever the healthy size happens to be.

The truncation failure mode (small files) is what actually occurs from crashed/
interrupted jobs, and the outlier check catches it cheaply. For the rarer
"full-sized but semantically empty" case (all bins below the score floor),
add --content-check to additionally parse each file and count above-floor bins
(slower: decompresses, but early-exits at --min-bins for healthy files).

It tells you which perm{N}/{target} to rebuild AND prints the cleanup commands,
including the perm{N}.obs.done / perm{N}.diff.done sentinels (which must be
removed or Snakemake treats the broken perm as complete and skips it).

Usage:
  python check_perm_set.py \
      --perm-dir workflow/results/perm --n-perms 100 \
      --samples DN_rep1 DN_rep2 DP_rep1 DP_rep2 EbKO_rep1 EbKO_rep2 \
                ProB_rep1 ProB_rep2 S3T3_rep1 S3T3_rep2 \
      --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 DP_vs_EbKO \
                    DP_vs_ProB DP_vs_S3T3 EbKO_vs_ProB EbKO_vs_S3T3 ProB_vs_S3T3
      [--content-check --min-signal 0.5]   # optional deep pass

Exit 0 if complete and intact; 1 otherwise (gate a launch:
`python check_perm_set.py ... && snakemake ...`). ASCII only.
"""
import argparse
import gzip
import json
import os
import sys
from statistics import median


def count_above_floor(path, min_signal, min_bins, diff_mode):
    """Count qcat bins with (|score| in diff mode) >= min_signal, early-exiting
    at min_bins. Mirrors bearing_pvalue.parse_qcat extraction."""
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            col = parts[3]
            if col.startswith("{"):
                pairs = json.loads(col).get("qcat", [])
            else:
                i = col.find("qcat:")
                if i < 0:
                    continue
                j = col.find(",raw:", i)
                payload = col[i + 5:j] if j >= 0 else col[i + 5:]
                pairs = json.loads(payload)
            score = sum(s for s, _ in pairs)
            s = abs(score) if diff_mode else score
            if s >= min_signal:
                n += 1
                if n >= min_bins:
                    return n
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perm-dir", required=True,
                    help="the perm/ directory (e.g. workflow/results/perm)")
    ap.add_argument("--n-perms", type=int, required=True)
    ap.add_argument("--samples", nargs="*", default=[],
                    help="per-sample perm qcats to check (omit to skip).")
    ap.add_argument("--comparisons", nargs="*", default=[],
                    help="differential perm qcats to check (omit to skip).")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="flag a file smaller than this fraction of its class's "
                         "median size (default 0.5 = half the median).")
    ap.add_argument("--abs-min-bytes", type=int, default=1024,
                    help="also flag any file below this absolute size, to catch "
                         "the case where MANY files are bad and skew the median "
                         "(default 1024).")
    ap.add_argument("--content-check", action="store_true",
                    help="additionally parse each file and require --min-bins "
                         "above-floor bins (catches full-size-but-empty files).")
    ap.add_argument("--min-signal", type=float, default=0.5,
                    help="score floor for --content-check (0.5 KL, 0.05 JSD).")
    ap.add_argument("--min-bins", type=int, default=1000,
                    help="min above-floor bins for --content-check (default 1000).")
    args = ap.parse_args()

    # Gather (path, diff_mode) for every expected file; record existence + size.
    expected = []  # (path, diff_mode)
    for p in range(1, args.n_perms + 1):
        for s in args.samples:
            expected.append((os.path.join(args.perm_dir, "perm%d" % p, s,
                             "%s_perm%d.qcat.bgz" % (s, p)), False))
        for c in args.comparisons:
            expected.append((os.path.join(args.perm_dir, "perm%d" % p,
                             "diff_comparison", "diff_%s.qcat.bgz" % c), True))

    missing = [p for p, _ in expected if not os.path.exists(p)]
    present = [(p, d) for p, d in expected if os.path.exists(p)]

    # Class medians (per-sample vs diff differ in size), so compare like with like.
    sizes_sample = [os.path.getsize(p) for p, d in present if not d]
    sizes_diff = [os.path.getsize(p) for p, d in present if d]
    med_sample = median(sizes_sample) if sizes_sample else 0
    med_diff = median(sizes_diff) if sizes_diff else 0

    outliers = []   # (path, size, why)
    for p, d in present:
        sz = os.path.getsize(p)
        med = med_diff if d else med_sample
        thresh = max(args.abs_min_bytes, args.min_frac * med)
        if sz < thresh:
            cls = "diff" if d else "sample"
            outliers.append((p, sz, "%d B < %.0f (%.0f%% of %s-median %d)"
                             % (sz, thresh, 100 * args.min_frac, cls, med)))

    lowcontent = []
    if args.content_check:
        bad_set = set(missing) | {p for p, _, _ in outliers}
        for p, d in present:
            if p in bad_set:
                continue
            try:
                c = count_above_floor(p, args.min_signal, args.min_bins, d)
            except Exception as e:
                lowcontent.append((p, "unreadable: %s" % e)); continue
            if c < args.min_bins:
                lowcontent.append((p, "%d bins above floor" % c))

    n_ok = len(present) - len(outliers) - len(lowcontent)
    mode = "size-outlier (<%.0f%% of class median)" % (100 * args.min_frac)
    if args.content_check:
        mode += " + content (floor=%g)" % args.min_signal
    print("Checked %d expected perm qcats across %d perms [%s]." %
          (len(expected), args.n_perms, mode))
    print("  median size: per-sample=%d B, diff=%d B" % (med_sample, med_diff))
    print("  %d OK, %d missing, %d size-outliers, %d low-content."
          % (n_ok, len(missing), len(outliers), len(lowcontent)))

    if missing:
        print("\nMISSING (%d):" % len(missing))
        for p in sorted(missing):
            print("  " + p)
    if outliers:
        print("\nSIZE OUTLIERS / TRUNCATED (%d):" % len(outliers))
        for p, sz, why in sorted(outliers):
            print("  %s  -- %s" % (p, why))
    if lowcontent:
        print("\nFULL-SIZE BUT EMPTY (%d):" % len(lowcontent))
        for p, why in sorted(lowcontent):
            print("  %s  -- %s" % (p, why))

    bad_paths = list(missing) + [x[0] for x in outliers] + [x[0] for x in lowcontent]
    if bad_paths:
        bad_perms = set()
        for p in bad_paths:
            for part in p.split(os.sep):
                if part.startswith("perm") and part[4:].isdigit():
                    bad_perms.add(int(part[4:])); break
        print("\nPerm indices to rebuild: %s"
              % ", ".join(str(i) for i in sorted(bad_perms)))
        print("\nCleanup commands (remove the bad outputs AND the .done sentinels,")
        print("else Snakemake treats the broken perm as complete and skips it):")
        for i in sorted(bad_perms):
            print("  rm -f  %s/perm%d.obs.done %s/perm%d.diff.done"
                  % (args.perm_dir, i, args.perm_dir, i))
        for p in sorted(bad_paths):
            print("  rm -f  %s" % p)
        print("\nRefusing to proceed -- rebuild the above before computing p-values.")
        sys.exit(1)

    print("\nPerm set is complete and intact. Safe to compute p-values.")
    sys.exit(0)


if __name__ == "__main__":
    main()
