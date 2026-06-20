#!/usr/bin/env python3
"""Pre-flight validator for the BEARING permutation null set.

Run this BEFORE launching the p-value rebuild. It verifies that every expected
permutation qcat -- per-sample (perm{N}/{sample}/{sample}_perm{N}.qcat.bgz) and
differential (perm{N}/diff_comparison/diff_{comp}.qcat.bgz) -- exists and is a
plausible, non-truncated size. It does NOT load the bins (that is what costs
hours); it stats the files, so it runs in seconds and tells you exactly which
perm{N}/{target} to rebuild.

This catches the failure mode where a single truncated perm file (e.g. an empty
perm34/DP_rep1) aborts a multi-hour p-value job at the very end. Run it, fix the
handful it reports, then launch -- the job then can't die on a bad perm.

Usage:
  python check_perm_set.py \
      --perm-dir workflow/results/perm \
      --n-perms 100 \
      --samples DN_rep1 DN_rep2 DP_rep1 DP_rep2 EbKO_rep1 EbKO_rep2 \
                ProB_rep1 ProB_rep2 S3T3_rep1 S3T3_rep2 \
      --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 DP_vs_EbKO \
                    DP_vs_ProB DP_vs_S3T3 EbKO_vs_ProB EbKO_vs_S3T3 ProB_vs_S3T3

Exit code 0 if the set is complete; 1 if anything is missing/too small (so it
can gate a launch script: `python check_perm_set.py ... && snakemake ...`).
ASCII only.
"""
import argparse
import os
import sys


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
    ap.add_argument("--min-bytes", type=int, default=1_000_000,
                    help="minimum plausible qcat size in bytes (default 1e6; a "
                         "healthy per-sample qcat is ~150M, a diff qcat ~300M, "
                         "so anything under 1M is truncated/empty).")
    args = ap.parse_args()

    missing = []   # does not exist
    toosmall = []  # exists but below min-bytes
    ok = 0

    def check(path):
        nonlocal ok
        if not os.path.exists(path):
            missing.append(path)
        elif os.path.getsize(path) < args.min_bytes:
            toosmall.append((path, os.path.getsize(path)))
        else:
            ok += 1

    for p in range(1, args.n_perms + 1):
        for s in args.samples:
            check(os.path.join(args.perm_dir, "perm%d" % p, s,
                               "%s_perm%d.qcat.bgz" % (s, p)))
        for c in args.comparisons:
            check(os.path.join(args.perm_dir, "perm%d" % p, "diff_comparison",
                               "diff_%s.qcat.bgz" % c))

    total = ok + len(missing) + len(toosmall)
    print("Checked %d expected perm qcats across %d perms: %d OK, %d missing, "
          "%d too-small." % (total, args.n_perms, ok, len(missing), len(toosmall)))

    if missing:
        print("\nMISSING (%d):" % len(missing))
        for p in sorted(missing):
            print("  " + p)
    if toosmall:
        print("\nTOO SMALL / TRUNCATED (%d):" % len(toosmall))
        for p, sz in sorted(toosmall):
            print("  %s  (%d bytes)" % (p, sz))

    if missing or toosmall:
        # Summarize which perm indices are implicated, for targeted rebuild.
        bad_perms = set()
        for p in missing + [x[0] for x in toosmall]:
            for part in p.split(os.sep):
                if part.startswith("perm") and part[4:].isdigit():
                    bad_perms.add(int(part[4:]))
                    break
        print("\nPerm indices to rebuild: %s"
              % ", ".join(str(i) for i in sorted(bad_perms)))
        print("Refusing to proceed -- rebuild the above before computing p-values.")
        sys.exit(1)

    print("\nPerm set is complete and intact. Safe to compute p-values.")
    sys.exit(0)


if __name__ == "__main__":
    main()
