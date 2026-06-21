#!/usr/bin/env python3
"""Assert a freshly-written qcat has real content; exit non-zero if not.

Used as a post-write gate in the Snakefile perm rule so a permutation qcat that
comes out FULL SIZE but with zero bins above the score floor (the perm34/DP_rep1
failure: 84 MB, 0 usable bins) fails AT CREATION instead of being caught later
by the p-value step. On failure Snakemake deletes the bad output, so the empty
perm can never be recorded as complete.

Counts bins with (|score| in --diff mode) >= --min-signal, early-exiting at
--min-bins, so healthy files (millions of bins) confirm in a fraction of a
second. Mirrors bearing_pvalue.parse_qcat extraction. ASCII only.

  python assert_qcat_nonempty.py FILE.qcat.bgz --min-signal 0.5 [--min-bins 1000] [--diff]
"""
import argparse
import gzip
import json
import os
import sys


def gzip_intact(path, chunk=1 << 20):
    """Verify the gzip stream is complete by draining to the end (the CRC trailer
    is at the END, so a truncation after a healthy start is only visible on a full
    read). No parsing, just byte decompression. Returns (ok, err)."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(chunk):
                pass
        return True, None
    except Exception as e:
        return False, "%s: %s" % (type(e).__name__, e)


def count_above_floor(path, min_signal, min_bins, diff_mode):
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
    ap.add_argument("qcat", help="qcat.bgz to validate")
    ap.add_argument("--min-signal", type=float, default=0.5,
                    help="score floor (0.5 KL, 0.05 JSD).")
    ap.add_argument("--min-bins", type=int, default=1000,
                    help="minimum above-floor bins required (default 1000).")
    ap.add_argument("--diff", action="store_true",
                    help="use |score| (differential qcat).")
    args = ap.parse_args()

    if not os.path.exists(args.qcat):
        sys.exit("ASSERT FAIL: qcat does not exist: %s" % args.qcat)
    # Integrity first: a truncated gzip (incomplete stream) would otherwise pass
    # an early-exit content check that stops before reaching the cut point.
    good, err = gzip_intact(args.qcat)
    if not good:
        sys.exit(
            "ASSERT FAIL: %s is a truncated/corrupt gzip stream (%s). The write "
            "did not complete. Failing so this output is not recorded as complete."
            % (args.qcat, err)
        )
    try:
        n = count_above_floor(args.qcat, args.min_signal, args.min_bins, args.diff)
    except Exception as e:
        sys.exit("ASSERT FAIL: could not read %s (%s)" % (args.qcat, e))
    if n < args.min_bins:
        sys.exit(
            "ASSERT FAIL: %s has only %d bins above the min-signal floor "
            "(%g); expected >= %d. The permutation produced no usable signal "
            "(full-size-but-empty). Failing so this output is not recorded as "
            "complete." % (args.qcat, n, args.min_signal, args.min_bins)
        )
    print("OK: %s intact, >= %d bins above floor %g"
          % (args.qcat, args.min_bins, args.min_signal))


if __name__ == "__main__":
    main()
