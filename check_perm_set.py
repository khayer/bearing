#!/usr/bin/env python3
"""Pre-flight validator for the BEARING permutation null set.

Run this BEFORE launching the p-value rebuild. It confirms every expected
permutation qcat exists -- per-sample (perm{N}/{sample}/{sample}_perm{N}.qcat.bgz)
and differential (perm{N}/diff_comparison/diff_{comp}.qcat.bgz) -- and that each
actually contains bins above the min-signal floor.

WHY CONTENT, NOT SIZE (learned the hard way): a broken perm qcat can be FULL
SIZE on disk yet have ZERO bins above the floor. Observed on real data:
perm34/DP_rep1 was 84 MB (about half the ~171 MB median) but contributed 0
usable bins. A size-outlier rule cannot separate these from legitimately small
healthy perms (perm qcat sizes vary with how much signal each random shift lands
in mappable regions), so it both misses real failures and false-positives on
fine files. The only reliable signal is the bin count, so the content check is
the default. It early-exits at --min-bins, so healthy files (millions of bins)
confirm in a fraction of a second each.

Use --size-prescreen for a fast stat-only first pass (flags gross truncation,
e.g. a few-KB file) before the content check -- useful as a quick sanity glance,
but it is NOT sufficient on its own and does not replace the content check.

It tells you which perm{N}/{target} to rebuild AND prints the cleanup commands,
including the perm{N}.obs.done / perm{N}.diff.done sentinels (which must be
removed or Snakemake treats the broken perm as complete and skips it).

Usage:
  python check_perm_set.py \
      --perm-dir workflow/results/perm --n-perms 100 \
      --samples DN_rep1 DN_rep2 DP_rep1 DP_rep2 EbKO_rep1 EbKO_rep2 \
                ProB_rep1 ProB_rep2 S3T3_rep1 S3T3_rep2 \
      --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 DP_vs_EbKO \
                    DP_vs_ProB DP_vs_S3T3 EbKO_vs_ProB EbKO_vs_S3T3 ProB_vs_S3T3 \
      --min-signal 0.5    # 0.5 for KL, 0.05 for JSD

Exit 0 if complete and intact; 1 otherwise (gate a launch:
`python check_perm_set.py ... && snakemake ...`). ASCII only.
"""
import argparse
import gzip
import json
import os
import sys


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
    ap.add_argument("--min-signal", type=float, default=0.5,
                    help="score floor for the content check (0.5 KL, 0.05 JSD).")
    ap.add_argument("--min-bins", type=int, default=1000,
                    help="minimum above-floor bins each file must contain "
                         "(default 1000; healthy files have millions).")
    ap.add_argument("--size-prescreen", action="store_true",
                    help="fast stat-only first pass (flags a few-KB truncation) "
                         "before the content check. A glance, not a substitute.")
    ap.add_argument("--prescreen-min-bytes", type=int, default=1000000,
                    help="prescreen floor in bytes (default 1e6).")
    args = ap.parse_args()

    expected = []  # (path, diff_mode)
    for p in range(1, args.n_perms + 1):
        for s in args.samples:
            expected.append((os.path.join(args.perm_dir, "perm%d" % p, s,
                             "%s_perm%d.qcat.bgz" % (s, p)), False))
        for c in args.comparisons:
            expected.append((os.path.join(args.perm_dir, "perm%d" % p,
                             "diff_comparison", "diff_%s.qcat.bgz" % c), True))

    missing, prescreen_fail, lowcontent = [], [], []
    ok = 0
    for path, diff_mode in expected:
        if not os.path.exists(path):
            missing.append(path); continue
        if args.size_prescreen and os.path.getsize(path) < args.prescreen_min_bytes:
            prescreen_fail.append((path, os.path.getsize(path))); continue
        try:
            c = count_above_floor(path, args.min_signal, args.min_bins, diff_mode)
        except Exception as e:
            lowcontent.append((path, "unreadable: %s" % e)); continue
        if c < args.min_bins:
            lowcontent.append((path, "%d bins above floor" % c)); continue
        ok += 1

    mode = "content (floor=%g, min_bins=%d)" % (args.min_signal, args.min_bins)
    if args.size_prescreen:
        mode = "size-prescreen + " + mode
    print("Checked %d expected perm qcats across %d perms [%s]: "
          "%d OK, %d missing, %d prescreen-fail, %d low-content."
          % (len(expected), args.n_perms, mode, ok, len(missing),
             len(prescreen_fail), len(lowcontent)))

    if missing:
        print("\nMISSING (%d):" % len(missing))
        for p in sorted(missing):
            print("  " + p)
    if prescreen_fail:
        print("\nTRUNCATED (prescreen, %d):" % len(prescreen_fail))
        for p, sz in sorted(prescreen_fail):
            print("  %s  -- %d bytes" % (p, sz))
    if lowcontent:
        print("\nEMPTY / LOW-CONTENT -- size OK but no usable bins (%d):"
              % len(lowcontent))
        for p, why in sorted(lowcontent):
            print("  %s  -- %s" % (p, why))

    bad_paths = list(missing) + [x[0] for x in prescreen_fail] + [x[0] for x in lowcontent]
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
