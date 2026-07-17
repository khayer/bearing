#!/usr/bin/env python3
"""
pvminsig_sweep.py -- differential-floor sensitivity, parallel and low-memory.

WHAT CHANGED FROM v1 (and why v1 was killed)
--------------------------------------------
v1 pooled every null score into one sorted array. At floor 0.1 that is ~11M bins
x 100 perms = ~1.1 BILLION float32 = 4.4 GB, doubled to ~8.8 GB by
np.concatenate. It OOM-killed after 8265 s of single-threaded parsing.

The array was never needed. Every floor only ever asks for the null SURVIVAL
FUNCTION at the observed scores and at the floors:

    p(s | f) = (count(null >= s) + 1) / (count(null >= f) + 1)

So use the sorted unique observed scores (plus the floors) as histogram edges
and accumulate null counts straight into them:

    hist[i]   = count of null scores in [edges[i], edges[i+1])
    suffix[i] = sum(hist[i:]) = count(null >= edges[i])        <- exact

This is EXACT, not a binned approximation, because the edges ARE the query
points: every value we ever need a count for is its own edge. Memory drops from
~8.8 GB to ~80 MB (one int64 counter per distinct observed score), and per-perm
work is independent, so it parallelises across files.

~5-10 min per comparison on 12 cores, in well under 1 GB.

CORRECTNESS
-----------
--verify compares the 0.5 arm against the production stats.tsv. Production runs
bearing_pvalue.py with no --min-signal, which resolves to 0.5 for KL, so the two
MUST agree on bins tested and bins significant. If they do not, trust
bearing_pvalue.py and discard this output.

USAGE -- on a compute node, NOT the login node
----------------------------------------------
    srun -p dbhiq,defq --mem=32G --cpus-per-task=12 --time=2:00:00 --pty \
      python pvminsig_sweep.py \
        --diff-dir workflow/results/compare \
        --perm-glob 'workflow/results/perm/perm*/diff_comparison' \
        --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 \
        --floors 0.1 0.25 0.5 1.0 1.5 \
        --threads 12 \
        --out sens/floor_sweep.tsv \
        --verify workflow/results/pvalue
"""

import argparse
import glob
import gzip
import multiprocessing as mp
import os
import re
import sys
import time

import numpy as np

SCORE_RE = re.compile(rb"\[(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),")

# Set before the Pool forks; workers inherit these copy-on-write on Linux.
_EDGES = None
_MIN_FLOOR = None
_CMP = None


def perm_label(d):
    """perm7 from '.../perm7/diff_comparison'.

    v1 printed os.path.basename(d), which was the literal string
    'diff_comparison' for all 100 dirs -- so 100 distinct files scrolled past
    under identical labels and looked like an infinite loop. It was not.
    """
    b = os.path.basename(d.rstrip("/"))
    if b in ("diff_comparison", ""):
        return os.path.basename(os.path.dirname(d.rstrip("/")))
    return b


def floor_tag(f):
    """Format a floor for filenames: 0.25 -> '0p25', 1.0 -> '1p0' (config_minsig
    style). Uses str() so the trailing zero in 1.0 is kept."""
    return str(f).replace(".", "p")


def read_observed_bins(path, min_abs):
    """Like iter_abs_scores but keeps per-bin coordinates and the SIGNED score,
    for --per-bin-out. Returns (abs_scores, chroms, starts, ends, signed), all
    aligned. Only used when per-bin output is requested (higher memory than the
    streaming abs-only path)."""
    absv = []
    chroms = []
    starts = []
    ends = []
    signed = []
    intern = sys.intern
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if not line or line[0:1] == b"#":
                continue
            parts = line.rstrip(b"\n").split(b"\t")
            if len(parts) < 4:
                continue
            col = parts[3]
            qs = col.find(b"qcat:")
            if qs == -1:
                continue
            rs = col.find(b",raw:", qs)
            payload = col[qs + 5:rs] if rs >= 0 else col[qs + 5:]
            total = sum(float(m.group(1)) for m in SCORE_RE.finditer(payload))
            a = abs(total)
            if a < min_abs:
                continue
            chroms.append(intern(parts[0].decode()))
            starts.append(int(parts[1]))
            ends.append(int(parts[2]))
            signed.append(total)
            absv.append(a)
    return (np.array(absv, dtype=np.float64), chroms,
            np.array(starts, dtype=np.int64), np.array(ends, dtype=np.int64),
            np.array(signed, dtype=np.float64))


def iter_abs_scores(path, min_abs, block=1000000):
    """Yield float64 arrays of abs(signed summed score) for bins >= min_abs.

    Mirrors bearing_pvalue.py parse_qcat in diff_mode: the score is the SIGNED
    sum of the per-track values; filtering uses its absolute value.
    """
    buf = []
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if not line or line[0:1] == b"#":
                continue
            parts = line.rstrip(b"\n").split(b"\t")
            if len(parts) < 4:
                continue
            col = parts[3]
            qs = col.find(b"qcat:")
            if qs == -1:
                continue
            rs = col.find(b",raw:", qs)
            payload = col[qs + 5:rs] if rs >= 0 else col[qs + 5:]
            # Use builtin sum(), NOT a naive `total += x` loop. bearing_pvalue.py
            # computes score_total = sum(per_track.values()), and CPython 3.12+
            # gives sum() Neumaier compensated summation. Naive accumulation of
            # the same six values in the same order can land 0.2 ULP lower --
            # enough to flip a bin sitting exactly on the floor. Observed once in
            # 5,013,001 bins (chr12:45070600, sums to exactly 0.5 under sum() and
            # to 0.49999999999999994 naively).
            total = sum(float(m.group(1)) for m in SCORE_RE.finditer(payload))
            a = abs(total)
            if a < min_abs:
                continue
            buf.append(a)
            if len(buf) >= block:
                yield np.array(buf, dtype=np.float64)
                buf = []
    if buf:
        yield np.array(buf, dtype=np.float64)


def _null_hist_worker(perm_dir):
    """One perm file -> histogram of its null scores over the shared edges."""
    path = os.path.join(perm_dir, "diff_%s.qcat.bgz" % _CMP)
    n_edges = len(_EDGES)
    hist = np.zeros(n_edges, dtype=np.int64)
    n = 0
    if not os.path.exists(path):
        return perm_label(perm_dir), hist, -1
    for chunk in iter_abs_scores(path, _MIN_FLOOR):
        idx = np.searchsorted(_EDGES, chunk, side="right") - 1
        np.clip(idx, 0, n_edges - 1, out=idx)
        hist += np.bincount(idx, minlength=n_edges)
        n += len(chunk)
    return perm_label(perm_dir), hist, n


def bh_fdr(pvals):
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ps = pvals[order]
    adj = np.minimum(1.0, ps * n / np.arange(1, n + 1))
    np.minimum.accumulate(adj[::-1], out=adj[::-1])
    out = np.empty(n, dtype=np.float64)
    out[order] = adj
    return out


def main():
    global _EDGES, _MIN_FLOOR, _CMP
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff-dir", required=True)
    ap.add_argument("--perm-glob", required=True)
    ap.add_argument("--comparisons", nargs="+", required=True)
    ap.add_argument("--floors", nargs="+", type=float,
                    default=[0.1, 0.25, 0.5, 1.0, 1.5])
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--max-perms", type=int, default=None,
                    help="Use only the first N perm dirs (smoke test). "
                         "--verify needs all of them.")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    ap.add_argument("--verify", default=None, metavar="PVALUE_DIR")
    ap.add_argument("--per-bin-out", default=None, metavar="DIR",
                    help="Write DIR/<CMP>_floor<TAG>.tsv.gz per comparison and "
                         "floor: chrom start end bearing_score pval pval_adj_bh "
                         "significant. Enables Panel C to cite a REAL flipping "
                         "bin. Written incrementally + gzipped.")
    args = ap.parse_args()
    if args.per_bin_out:
        os.makedirs(args.per_bin_out, exist_ok=True)

    floors = sorted(args.floors)
    min_floor = floors[0]
    perm_dirs = sorted(glob.glob(args.perm_glob))
    if args.max_perms:
        perm_dirs = perm_dirs[:args.max_perms]
    if not perm_dirs:
        print("ERROR: --perm-glob matched nothing", file=sys.stderr)
        return 1

    print("Perm dirs : %d" % len(perm_dirs), file=sys.stderr)
    print("Threads   : %d" % args.threads, file=sys.stderr)
    print("Floors    : %s" % ", ".join("%g" % f for f in floors), file=sys.stderr)
    if args.max_perms:
        print("WARNING: --max-perms %d -- p-values NOT comparable to production; "
              "--verify disabled." % args.max_perms, file=sys.stderr)

    rows = []
    for cmp_name in args.comparisons:
        t0 = time.time()
        print("\n=== %s ===" % cmp_name, file=sys.stderr)
        obs_path = os.path.join(args.diff_dir, "diff_%s.qcat.bgz" % cmp_name)
        if not os.path.exists(obs_path):
            print("  ERROR: missing %s" % obs_path, file=sys.stderr)
            continue

        print("  Reading observed...", file=sys.stderr)
        obs_chrom = obs_start = obs_end = obs_signed = None
        if args.per_bin_out:
            (obs_abs, obs_chrom, obs_start,
             obs_end, obs_signed) = read_observed_bins(obs_path, min_floor)
        else:
            obs_chunks = list(iter_abs_scores(obs_path, min_floor))
            obs_abs = np.concatenate(obs_chunks) if obs_chunks else np.array([])
            del obs_chunks
        if len(obs_abs) == 0:
            print("  no observed bins >= %g" % min_floor, file=sys.stderr)
            continue
        print("    %s observed bins >= %g" % (format(len(obs_abs), ","), min_floor),
              file=sys.stderr)

        edges = np.unique(np.concatenate(
            [obs_abs, np.array(floors, dtype=np.float64)]))
        print("    %s distinct edges (%.0f MB of counters)"
              % (format(len(edges), ","), len(edges) * 8 / 1e6), file=sys.stderr)

        _EDGES = edges
        _MIN_FLOOR = min_floor
        _CMP = cmp_name

        print("  Null histogram: %d perms on %d cores..."
              % (len(perm_dirs), args.threads), file=sys.stderr)
        total_hist = np.zeros(len(edges), dtype=np.int64)
        n_null_rows = 0
        done = 0
        with mp.Pool(processes=args.threads) as pool:
            for label, hist, n in pool.imap_unordered(_null_hist_worker, perm_dirs):
                done += 1
                if n < 0:
                    print("    [%3d/%3d] %-10s MISSING FILE"
                          % (done, len(perm_dirs), label), file=sys.stderr)
                    continue
                total_hist += hist
                n_null_rows += n
                if done % 10 == 0 or done == len(perm_dirs):
                    print("    [%3d/%3d] %-10s %s null bins; %.0f s elapsed"
                          % (done, len(perm_dirs), label, format(n, ","),
                             time.time() - t0), file=sys.stderr)
        print("    pooled null: %s scores in %.0f s"
              % (format(n_null_rows, ","), time.time() - t0), file=sys.stderr)

        suffix = np.cumsum(total_hist[::-1])[::-1]   # suffix[i] = count(null >= edges[i])

        for f in floors:
            keep = obs_abs >= f
            s_kept = obs_abs[keep]
            if len(s_kept) == 0:
                rows.append((cmp_name, f, 0, 0, 0, float("nan")))
                continue
            fi = int(np.searchsorted(edges, f, side="left"))
            n_null_f = int(suffix[fi])
            ji = np.searchsorted(edges, s_kept, side="left")
            cnt = suffix[ji]
            p = (cnt + 1.0) / (n_null_f + 1.0)
            p = np.clip(p, 1.0 / (n_null_f + 1.0), 1.0)
            q = bh_fdr(p)
            n_sig = int((q <= args.fdr).sum())
            rows.append((cmp_name, f, len(s_kept), n_null_f, n_sig, float(p.min())))

            if args.per_bin_out:
                # p/q are in obs_abs[keep] order == idx_kept order.
                idx_kept = np.nonzero(keep)[0]
                sig = (q <= args.fdr)
                pb = os.path.join(args.per_bin_out,
                                  "%s_floor%s.tsv.gz" % (cmp_name, floor_tag(f)))
                with gzip.open(pb, "wt") as pf:   # incremental, gzipped
                    pf.write("chrom\tstart\tend\tbearing_score\tpval\t"
                             "pval_adj_bh\tsignificant\n")
                    for k in range(len(idx_kept)):
                        oi = int(idx_kept[k])
                        pf.write("%s\t%d\t%d\t%.6g\t%.6g\t%.6g\t%d\n" % (
                            obs_chrom[oi], obs_start[oi], obs_end[oi],
                            obs_signed[oi], p[k], q[k], 1 if sig[k] else 0))
            print("    floor %-5g tested=%-12s null=%-15s sig=%-8s min_p=%.3g"
                  % (f, format(len(s_kept), ","), format(n_null_f, ","),
                     format(n_sig, ","), p.min()), file=sys.stderr)

            if args.verify and abs(f - 0.5) < 1e-9 and not args.max_perms:
                prod = os.path.join(args.verify, "diff_%s.stats.tsv" % cmp_name)
                if os.path.exists(prod):
                    n_prod = n_prod_sig = 0
                    with open(prod) as fh:
                        next(fh)
                        for line in fh:
                            n_prod += 1
                            if line.rstrip("\n").split("\t")[9] == "1":
                                n_prod_sig += 1
                    ok_n = n_prod == len(s_kept)
                    ok_s = n_prod_sig == n_sig
                    print("    VERIFY: tested %s vs prod %s [%s]; sig %s vs prod %s [%s]"
                          % (format(len(s_kept), ","), format(n_prod, ","),
                             "OK" if ok_n else "MISMATCH",
                             format(n_sig, ","), format(n_prod_sig, ","),
                             "OK" if ok_s else "MISMATCH"), file=sys.stderr)
                    if not (ok_n and ok_s):
                        print("    -> the 0.5 arm MUST reproduce production. Do not "
                              "trust the other floors until it does.", file=sys.stderr)

        del obs_abs, edges, total_hist, suffix
        print("  %s done in %.0f s" % (cmp_name, time.time() - t0), file=sys.stderr)

    outdir = os.path.dirname(args.out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("comparison\tfloor\tbins_tested\tnull_scores\t"
                 "bins_bh_significant\tmin_pval\n")
        for r in rows:
            fh.write("%s\t%g\t%d\t%d\t%d\t%.6g\n" % r)
    print("\nWrote %s" % args.out, file=sys.stderr)

    print("\n=== bins BH-significant by floor ===", file=sys.stderr)
    print("%-14s" % "comparison" + "".join("%14s" % ("floor %g" % f) for f in floors),
          file=sys.stderr)
    for c in sorted(set(r[0] for r in rows)):
        line = "%-14s" % c
        for f in floors:
            m = [r for r in rows if r[0] == c and abs(r[1] - f) < 1e-9]
            line += "%14s" % (format(m[0][4], ",") if m else "-")
        print(line, file=sys.stderr)
    print("", file=sys.stderr)
    print("Stable across floors -> the inherited 0.5 is not driving the calls, and",
          file=sys.stderr)
    print("one Methods sentence closes it. Unstable -> the floor is a real analytic",
          file=sys.stderr)
    print("choice that needs justifying before submission.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
