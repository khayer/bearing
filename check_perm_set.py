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


def gzip_intact(path, chunk=1 << 20):
    """Verify a gzip stream is complete by draining it to the end (the CRC/length
    trailer is at the END, so this MUST read the whole file -- an early-exiting
    content check cannot see a truncation that occurs after its first chunk, e.g.
    perm85 loaded 13.6M bins fine then truncated). No JSON parsing, just byte
    decompression, so it is faster than the content count. Returns (ok, err)."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(chunk):
                pass
        return True, None
    except Exception as e:
        return False, "%s: %s" % (type(e).__name__, e)


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
    ap.add_argument("--integrity", action="store_true",
                    help="ALSO verify each file's gzip stream is complete (full "
                         "decompress, no parse). Catches truncation anywhere in "
                         "the file -- including after a healthy start, which the "
                         "early-exit content check misses. Slower (reads every "
                         "byte); RECOMMENDED before a long run.")
    ap.add_argument("--calibration-dir", default=None,
                    help="also check the calibration tree (e.g. "
                         "workflow/results/calibration). Requires --conditions.")
    ap.add_argument("--conditions", nargs="*", default=[],
                    help="conditions for --calibration-dir (e.g. DN DP EbKO "
                         "ProB S3T3). Calibration diff qcats truncate (the "
                         "perm85 EOFError), so pair this with --integrity.")
    args = ap.parse_args()

    expected = []  # (path, diff_mode)
    for p in range(1, args.n_perms + 1):
        for s in args.samples:
            expected.append((os.path.join(args.perm_dir, "perm%d" % p, s,
                             "%s_perm%d.qcat.bgz" % (s, p)), False))
        for c in args.comparisons:
            expected.append((os.path.join(args.perm_dir, "perm%d" % p,
                             "diff_comparison", "diff_%s.qcat.bgz" % c), True))

    # Calibration tree: per condition, an observed self-vs-self replicate diff
    # qcat plus one per permutation, at
    #   {cal}/{cond}/{cond}/observed/diff_{cond}_rep1_vs_{cond}_rep2.qcat.bgz
    #   {cal}/{cond}/{cond}/perm/perm{N}/diff_{cond}_rep1_vs_{cond}_rep2.qcat.bgz
    # These are differential qcats (diff_mode=True) and their failure mode is
    # truncation, so --integrity matters most here.
    if args.calibration_dir:
        for cond in args.conditions:
            base = os.path.join(args.calibration_dir, cond, cond)
            fn = "diff_%s_rep1_vs_%s_rep2.qcat.bgz" % (cond, cond)
            expected.append((os.path.join(base, "observed", fn), True))
            for p in range(1, args.n_perms + 1):
                expected.append((os.path.join(base, "perm", "perm%d" % p, fn), True))

    missing, prescreen_fail, lowcontent, corrupt = [], [], [], []
    ok = 0
    for path, diff_mode in expected:
        if not os.path.exists(path):
            missing.append(path); continue
        if args.size_prescreen and os.path.getsize(path) < args.prescreen_min_bytes:
            prescreen_fail.append((path, os.path.getsize(path))); continue
        if args.integrity:
            good, err = gzip_intact(path)
            if not good:
                corrupt.append((path, err)); continue
        try:
            c = count_above_floor(path, args.min_signal, args.min_bins, diff_mode)
        except Exception as e:
            corrupt.append((path, "unreadable: %s" % e)); continue
        if c < args.min_bins:
            lowcontent.append((path, "%d bins above floor" % c)); continue
        ok += 1

    mode = "content (floor=%g, min_bins=%d)" % (args.min_signal, args.min_bins)
    if args.integrity:
        mode = "integrity + " + mode
    if args.size_prescreen:
        mode = "size-prescreen + " + mode
    if args.calibration_dir:
        mode += " [+calibration]"
    print("Checked %d expected qcats across %d perms [%s]: "
          "%d OK, %d missing, %d prescreen-fail, %d corrupt, %d low-content."
          % (len(expected), args.n_perms, mode, ok, len(missing),
             len(prescreen_fail), len(corrupt), len(lowcontent)))

    if missing:
        print("\nMISSING (%d):" % len(missing))
        for p in sorted(missing):
            print("  " + p)
    if prescreen_fail:
        print("\nTRUNCATED (prescreen, %d):" % len(prescreen_fail))
        for p, sz in sorted(prescreen_fail):
            print("  %s  -- %d bytes" % (p, sz))
    if corrupt:
        print("\nCORRUPT / TRUNCATED GZIP -- incomplete stream (%d):" % len(corrupt))
        for p, why in sorted(corrupt):
            print("  %s  -- %s" % (p, why))
    if lowcontent:
        print("\nEMPTY / LOW-CONTENT -- size OK but no usable bins (%d):"
              % len(lowcontent))
        for p, why in sorted(lowcontent):
            print("  %s  -- %s" % (p, why))

    bad_paths = (list(missing) + [x[0] for x in prescreen_fail]
                 + [x[0] for x in corrupt] + [x[0] for x in lowcontent])
    if bad_paths:
        # Split perm-tree files (need .done sentinel cleanup so Snakemake rebuilds)
        # from calibration-tree files (rebuilt by the calibration rule; no perm
        # sentinels apply). Both live under .../perm{N}/..., so distinguish by
        # whether the path is under --calibration-dir.
        cal_dir = os.path.normpath(args.calibration_dir) if args.calibration_dir else None
        perm_bad = [p for p in bad_paths
                    if not (cal_dir and os.path.normpath(p).startswith(cal_dir))]
        cal_bad = [p for p in bad_paths if p not in perm_bad]

        if perm_bad:
            bad_perms = set()
            for p in perm_bad:
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
            for p in sorted(perm_bad):
                print("  rm -f  %s" % p)

        if cal_bad:
            print("\nCalibration files to rebuild (delete, then re-run the")
            print("calibration rule for the affected condition(s)):")
            for p in sorted(cal_bad):
                print("  rm -f  %s" % p)

        print("\nRefusing to proceed -- rebuild the above before computing p-values.")
        sys.exit(1)

    print("\nAll checked qcats are complete and intact. Safe to proceed.")
    sys.exit(0)


if __name__ == "__main__":
    main()
