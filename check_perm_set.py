#!/usr/bin/env python3
"""Pre-flight validator for the BEARING permutation null set.

Run this BEFORE launching the p-value rebuild. It verifies that every expected
permutation qcat -- per-sample (perm{N}/{sample}/{sample}_perm{N}.qcat.bgz) and
differential (perm{N}/diff_comparison/diff_{comp}.qcat.bgz) -- exists, is a
plausible size, AND actually contains bins above the min-signal floor.

WHY THE CONTENT CHECK MATTERS: a qcat can be >1 MB on disk yet contribute ZERO
bins above the floor (the perm34/DP_rep1 failure). A size-only check passes such
a file; the scorer then aborts on it after loading 99 good files. This validator
counts parseable above-floor bins (early-exiting at --min-bins for healthy files,
so it stays fast) and flags any file that comes up short.

It tells you exactly which perm{N}/{target} to rebuild, AND prints the cleanup
commands -- including the perm{N}.obs.done / perm{N}.diff.done sentinels, which
must be removed or Snakemake will treat the broken perm as complete and skip it.

Usage:
  python check_perm_set.py \
      --perm-dir workflow/results/perm --n-perms 100 \
      --samples DN_rep1 DN_rep2 DP_rep1 DP_rep2 EbKO_rep1 EbKO_rep2 \
                ProB_rep1 ProB_rep2 S3T3_rep1 S3T3_rep2 \
      --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 DP_vs_EbKO \
                    DP_vs_ProB DP_vs_S3T3 EbKO_vs_ProB EbKO_vs_S3T3 ProB_vs_S3T3 \
      --min-signal 0.5   # 0.5 for KL, 0.05 for JSD

Exit 0 if complete and intact; 1 otherwise (so it can gate a launch:
`python check_perm_set.py ... && snakemake ...`). ASCII only.
"""
import argparse
import gzip
import json
import os
import sys


def count_above_floor(path, min_signal, min_bins, diff_mode):
    """Count qcat bins with (|score| in diff mode) >= min_signal, early-exiting
    once min_bins is reached. Mirrors bearing_pvalue.parse_qcat extraction."""
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
    ap.add_argument("--min-bytes", type=int, default=1000000,
                    help="minimum plausible qcat size in bytes (default 1e6).")
    ap.add_argument("--min-signal", type=float, default=0.5,
                    help="score floor for the content check (0.5 KL, 0.05 JSD).")
    ap.add_argument("--min-bins", type=int, default=1000,
                    help="minimum above-floor bins each file must contain "
                         "(default 1000; healthy files have millions).")
    ap.add_argument("--size-only", action="store_true",
                    help="skip the (slower) content check; stat sizes only.")
    args = ap.parse_args()

    missing, toosmall, lowcontent = [], [], []
    ok = [0]

    def check(path, diff_mode):
        if not os.path.exists(path):
            missing.append(path); return
        if os.path.getsize(path) < args.min_bytes:
            toosmall.append((path, os.path.getsize(path))); return
        if not args.size_only:
            try:
                c = count_above_floor(path, args.min_signal, args.min_bins, diff_mode)
            except Exception as e:
                lowcontent.append((path, "unreadable: %s" % e)); return
            if c < args.min_bins:
                lowcontent.append((path, "%d bins above floor" % c)); return
        ok[0] += 1

    for p in range(1, args.n_perms + 1):
        for s in args.samples:
            check(os.path.join(args.perm_dir, "perm%d" % p, s,
                               "%s_perm%d.qcat.bgz" % (s, p)), diff_mode=False)
        for c in args.comparisons:
            check(os.path.join(args.perm_dir, "perm%d" % p, "diff_comparison",
                               "diff_%s.qcat.bgz" % c), diff_mode=True)

    total = ok[0] + len(missing) + len(toosmall) + len(lowcontent)
    mode = "size-only" if args.size_only else ("content, floor=%g" % args.min_signal)
    print("Checked %d expected perm qcats across %d perms [%s]: "
          "%d OK, %d missing, %d too-small, %d low-content."
          % (total, args.n_perms, mode, ok[0], len(missing),
             len(toosmall), len(lowcontent)))

    for label, items in (("MISSING", [(p, "") for p in missing]),
                         ("TOO SMALL / TRUNCATED", toosmall),
                         ("EMPTY / LOW-CONTENT (size OK, no usable bins)", lowcontent)):
        if items:
            print("\n%s (%d):" % (label, len(items)))
            for p, why in sorted(items):
                print("  %s%s" % (p, ("  -- %s" % why) if why != "" else ""))

    bad_paths = list(missing) + [x[0] for x in toosmall] + [x[0] for x in lowcontent]
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
