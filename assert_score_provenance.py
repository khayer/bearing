#!/usr/bin/env python3
"""Verify (or stamp) per-sample scoring-provenance signatures against the active
config, and HALT the pipeline (exit 1) on any mismatch.

This is the data-layer backstop for the stale-qcat-reuse failure: a regime whose
normalize / score-method / adaptive-binning differs from a previous run can
silently reuse the previous run's qcats (because those settings live in
Snakemake params, not file inputs, and rerun-triggers=mtime ignores param
changes). Each qcat carries a sidecar <qcat>.sig stamped by bigwig_to_qcat.py
recording the settings it ACTUALLY scored with; this script recomputes the
EXPECTED signature from the settings the active config resolves to and compares.

The expected settings are passed in as the SAME resolved flag string the score
rule builds (--norm-flags), plus the other scoring axes, so there is no second,
divergent re-resolution of the config here.

Modes:
  verify (default) -- compare each <outdir>/<sample>.qcat.bgz.sig to expected;
                      mismatch -> exit 1. Missing .sig -> exit 1 unless
                      --allow-missing (transition aid for pre-guard trees).
  --stamp          -- write the expected signature to each .sig. Use ONLY on a
                      tree you have INDEPENDENTLY verified is correctly scored
                      (one-time backfill for trees built before this guard).

Exit codes: 0 = all good (or stamped); 1 = mismatch / missing / error.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_provenance import score_provenance_signature


def parse_norm_flags(s):
    """Derive (normalize_tracks, normalize_method, cohort_reference) from the
    resolved normalize-flag string the score rule passes to bigwig_to_qcat.py,
    interpreting it exactly as that script's argparse would."""
    toks = (s or "").split()
    normalize_tracks = "--normalize-tracks" in toks
    normalize_method = "nonzero-quantile"
    if "--normalize-method" in toks:
        i = toks.index("--normalize-method")
        if i + 1 < len(toks):
            normalize_method = toks[i + 1]
    cohort_reference = None
    if "--cohort-reference" in toks:
        i = toks.index("--cohort-reference")
        if i + 1 < len(toks):
            cohort_reference = toks[i + 1]
    return normalize_tracks, normalize_method, cohort_reference


def sig_path_for(outdir, sample):
    return os.path.join(outdir, "%s.qcat.bgz.sig" % sample)


def read_sig_digest(path):
    """First non-comment, non-blank line of a .sig file (the digest)."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                return line
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True,
                    help="Results dir holding <sample>.qcat.bgz[.sig] files.")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--norm-flags", default="",
                    help="Resolved normalize-flag string the score rule used "
                         "(e.g. '--normalize-tracks --normalize-method nonzero-quantile').")
    ap.add_argument("--score-method", default="kl")
    ap.add_argument("--min-signal", type=float, required=True)
    ap.add_argument("--categories", required=True)
    ap.add_argument("--bins-bed", default="",
                    help="Adaptive-binning bed path, or empty for fixed binning.")
    ap.add_argument("--stamp", action="store_true",
                    help="Write expected sig to each .sig (trusted backfill only).")
    ap.add_argument("--allow-missing", action="store_true",
                    help="Treat a missing .sig as a warning instead of a failure.")
    args = ap.parse_args()

    ntr, meth, cref = parse_norm_flags(args.norm_flags)
    expected, payload = score_provenance_signature(
        ntr, meth, args.score_method, args.min_signal,
        args.categories, (args.bins_bed or None), cref)

    if args.stamp:
        for s in args.samples:
            p = sig_path_for(args.outdir, s)
            with open(p, "w") as fh:
                fh.write(expected + "\n")
                for ln in payload.splitlines():
                    fh.write("# " + ln + "\n")
            print("stamped %s" % p)
        print("Stamped %d signatures with expected = %s" % (len(args.samples), expected))
        return 0

    mismatches = []
    missing = []
    for s in args.samples:
        p = sig_path_for(args.outdir, s)
        if not os.path.exists(p):
            missing.append(s)
            continue
        got = read_sig_digest(p)
        if got != expected:
            mismatches.append((s, got))

    if missing:
        msg = "missing scoring-provenance signatures for: %s" % ", ".join(missing)
        if args.allow_missing:
            print("WARNING: %s (--allow-missing set, not failing)" % msg)
        else:
            sys.stderr.write(
                "SCORE-PROVENANCE FAIL: %s\n"
                "  These qcats have no .sig and cannot be verified. Re-score them, "
                "or backfill with --stamp ONLY if you have verified they are "
                "correctly scored.\n" % msg)
            return 1

    if mismatches:
        sys.stderr.write(
            "SCORE-PROVENANCE FAIL: %d qcat(s) were scored with settings that do "
            "NOT match the active config.\n" % len(mismatches))
        sys.stderr.write("  expected signature: %s\n" % expected)
        for ln in payload.splitlines():
            sys.stderr.write("    %s\n" % ln)
        for s, got in mismatches:
            sys.stderr.write("  %-14s on-disk sig: %s\n" % (s, got))
        sys.stderr.write(
            "  This is the stale-qcat-reuse trap: the qcats were NOT re-scored for "
            "this regime. Delete the regime's score layer (rm the .qcat.bgz files "
            "or the whole outdir) and re-run so scoring is forced.\n")
        return 1

    print("score-provenance OK: %d/%d qcats match expected signature %s"
          % (len(args.samples) - len(missing), len(args.samples), expected))
    return 0


if __name__ == "__main__":
    sys.exit(main())
