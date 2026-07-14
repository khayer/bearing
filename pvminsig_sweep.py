#!/usr/bin/env python3
"""
pvminsig_sweep.py -- differential-floor sensitivity in ONE pass over the null.

WHY THIS EXISTS
---------------
Running bearing_pvalue.py once per floor re-reads the whole permutation null
each time. With 100 perms x ~5M bins that is ~500M JSON parses PER FLOOR, and
the null is bit-identical across floors -- only the threshold applied to it
changes. Five floors x four comparisons = 20 full null reads for what is really
four.

This script reads the observed diff qcat once and the null once per comparison,
keeps the pooled null as a single sorted array, and derives every floor from it
by binary search. The results are EXACT, not approximate:

    null_f            = N[N >= f]                       (a suffix of sorted N)
    n_null_f          = len(N) - searchsorted(N, f)
    count(null_f >= s) = len(N) - searchsorted(N, s)     for s >= f
    p(s | f)          = (count + 1) / (n_null_f + 1)     (Davison-Hinkley)

which is exactly what bearing_pvalue.py computes when handed a null collected
at floor f (bearing_pvalue.py: empirical_pvals, lines ~258-275).

CORRECTNESS
-----------
--verify compares this script's floor-0.5 output against the production
stats.tsv. Production uses the bearing_pvalue.py default floor for KL, which is
0.5, so the two must agree. If they do not, trust bearing_pvalue.py and discard
this script's output.

NOTE ON SUBSAMPLING
-------------------
Do not subsample the null to save time. With ~5M tests, BH significance needs
p ~1e-8, and the smallest attainable p is 1/(n_null+1). A null of 500M supports
p ~2e-9; a subsampled null of 5M caps p at 2e-7, and nothing can ever be called
significant. The null must stay large; that is why the fix is fewer passes, not
a smaller null.

USAGE
-----
    python pvminsig_sweep.py \
        --diff-dir workflow/results/compare \
        --perm-glob 'workflow/results/perm/perm*/diff_comparison' \
        --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 \
        --floors 0.1 0.25 0.5 1.0 1.5 \
        --out sens/floor_sweep.tsv \
        --verify workflow/results/pvalue
"""

import argparse
import glob
import gzip
import os
import sys
import time

import numpy as np

# Matches the leading score of each [score, track] pair: [-0.158967,2]
import re
SCORE_RE = re.compile(rb"\[(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),")


def parse_scores(path, min_abs):
    """
    Yield (chrom, start, signed_score) for bins with abs(score) >= min_abs.

    Mirrors bearing_pvalue.py parse_qcat in diff_mode: score_total is the SIGNED
    sum of the per-track values, and filtering uses its absolute value. Parsing
    is done on bytes with a regex rather than json.loads, which is materially
    faster over hundreds of millions of rows and yields identical sums (the
    payload is a flat [[score,track],...] list).
    """
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
            total = 0.0
            for m in SCORE_RE.finditer(payload):
                total += float(m.group(1))
            if abs(total) < min_abs:
                continue
            yield parts[0], int(parts[1]), total


def load_null(perm_dirs, cmp_name, min_floor):
    """Pool |score| >= min_floor across all perm dirs into one sorted array."""
    chunks = []
    total_rows = 0
    for d in perm_dirs:
        p = os.path.join(d, "diff_%s.qcat.bgz" % cmp_name)
        if not os.path.exists(p):
            print("    WARNING: missing %s" % p, file=sys.stderr)
            continue
        vals = np.fromiter(
            (abs(s) for _, _, s in parse_scores(p, min_floor)),
            dtype=np.float32,
        )
        chunks.append(vals)
        total_rows += len(vals)
        print("    %s: %s bins >= %g" % (os.path.basename(d),
                                         format(len(vals), ","), min_floor),
              file=sys.stderr)
    if not chunks:
        return None
    N = np.concatenate(chunks)
    del chunks
    N.sort()
    return N


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values (same procedure as the pipeline)."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ps = pvals[order]
    ranks = np.arange(1, n + 1)
    adj = np.minimum(1.0, ps * n / ranks)
    # enforce monotonicity from the top down
    for i in range(n - 2, -1, -1):
        if adj[i] > adj[i + 1]:
            adj[i] = adj[i + 1]
    out = np.empty(n, dtype=np.float64)
    out[order] = adj
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff-dir", required=True,
                    help="Directory holding the observed diff_<CMP>.qcat.bgz.")
    ap.add_argument("--perm-glob", required=True,
                    help="Glob for the per-perm diff directories, quoted.")
    ap.add_argument("--comparisons", nargs="+", required=True)
    ap.add_argument("--floors", nargs="+", type=float,
                    default=[0.1, 0.25, 0.5, 1.0, 1.5])
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out", required=True, help="Summary TSV to write.")
    ap.add_argument("--per-bin-dir", default=None,
                    help="Optional: write per-bin results per floor here.")
    ap.add_argument("--verify", default=None, metavar="PVALUE_DIR",
                    help="Compare the 0.5 floor against production stats.tsv "
                         "in this directory. Production uses the KL default "
                         "floor of 0.5, so they must agree.")
    args = ap.parse_args()

    floors = sorted(args.floors)
    min_floor = floors[0]
    perm_dirs = sorted(glob.glob(args.perm_glob))
    if not perm_dirs:
        print("ERROR: --perm-glob matched nothing: %s" % args.perm_glob,
              file=sys.stderr)
        return 1
    print("Perm dirs: %d" % len(perm_dirs), file=sys.stderr)
    print("Floors: %s (null collected once at %g)"
          % (", ".join("%g" % f for f in floors), min_floor), file=sys.stderr)

    rows = []
    for cmp_name in args.comparisons:
        t0 = time.time()
        print("\n=== %s ===" % cmp_name, file=sys.stderr)

        obs_path = os.path.join(args.diff_dir, "diff_%s.qcat.bgz" % cmp_name)
        if not os.path.exists(obs_path):
            print("  ERROR: missing %s" % obs_path, file=sys.stderr)
            continue
        print("  Reading observed...", file=sys.stderr)
        chrom_l, start_l, score_l = [], [], []
        for c, s, v in parse_scores(obs_path, min_floor):
            chrom_l.append(c)
            start_l.append(s)
            score_l.append(v)
        obs_score = np.array(score_l, dtype=np.float64)
        obs_abs = np.abs(obs_score)
        obs_start = np.array(start_l, dtype=np.int64)
        obs_chrom = np.array([c.decode() for c in chrom_l])
        del chrom_l, start_l, score_l
        print("    %s observed bins >= %g" % (format(len(obs_abs), ","), min_floor),
              file=sys.stderr)

        print("  Reading null (ONCE, reused for all floors)...", file=sys.stderr)
        N = load_null(perm_dirs, cmp_name, min_floor)
        if N is None:
            print("  ERROR: no null loaded for %s" % cmp_name, file=sys.stderr)
            continue
        n_tot = len(N)
        print("    pooled null: %s scores; %.1f s so far"
              % (format(n_tot, ","), time.time() - t0), file=sys.stderr)

        for f in floors:
            keep = obs_abs >= f
            s_kept = obs_abs[keep]
            if len(s_kept) == 0:
                rows.append((cmp_name, f, 0, 0, 0, float("nan")))
                continue
            k_f = int(np.searchsorted(N, f, side="left"))
            n_null_f = n_tot - k_f
            cnt = n_tot - np.searchsorted(N, s_kept, side="left")
            p = (cnt + 1.0) / (n_null_f + 1.0)
            p = np.clip(p, 1.0 / (n_null_f + 1.0), 1.0)
            q = bh_fdr(p)
            n_sig = int((q <= args.fdr).sum())
            rows.append((cmp_name, f, len(s_kept), n_null_f, n_sig,
                         float(p.min())))
            print("    floor %-5g tested=%-12s null=%-13s sig=%-8s min_p=%.3g"
                  % (f, format(len(s_kept), ","), format(n_null_f, ","),
                     format(n_sig, ","), p.min()), file=sys.stderr)

            if args.per_bin_dir:
                os.makedirs(args.per_bin_dir, exist_ok=True)
                tag = ("%g" % f).replace(".", "p")
                out_p = os.path.join(args.per_bin_dir,
                                     "%s_floor%s.tsv" % (cmp_name, tag))
                with open(out_p, "w") as fh:
                    fh.write("chrom\tstart\tbearing_score\tpval\tpval_adj_bh\tsignificant\n")
                    ss = obs_score[keep]
                    cc = obs_chrom[keep]
                    st = obs_start[keep]
                    for i in range(len(ss)):
                        fh.write("%s\t%d\t%.6f\t%.6e\t%.6e\t%d\n"
                                 % (cc[i], st[i], ss[i], p[i], q[i],
                                    1 if q[i] <= args.fdr else 0))

            if args.verify and abs(f - 0.5) < 1e-9:
                prod = os.path.join(args.verify, "diff_%s.stats.tsv" % cmp_name)
                if os.path.exists(prod):
                    n_prod = 0
                    n_prod_sig = 0
                    with open(prod) as fh:
                        next(fh)
                        for line in fh:
                            fl = line.rstrip("\n").split("\t")
                            n_prod += 1
                            if fl[9] == "1":
                                n_prod_sig += 1
                    ok_n = (n_prod == len(s_kept))
                    ok_s = (n_prod_sig == n_sig)
                    print("    VERIFY vs production: tested %s vs %s [%s]; "
                          "significant %s vs %s [%s]"
                          % (format(len(s_kept), ","), format(n_prod, ","),
                             "OK" if ok_n else "MISMATCH",
                             format(n_sig, ","), format(n_prod_sig, ","),
                             "OK" if ok_s else "MISMATCH"),
                          file=sys.stderr)
                    if not (ok_n and ok_s):
                        print("    -> floor 0.5 must reproduce production. "
                              "Investigate before trusting other floors.",
                              file=sys.stderr)
        del N, obs_abs, obs_score, obs_start, obs_chrom
        print("  %s done in %.1f s" % (cmp_name, time.time() - t0), file=sys.stderr)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("comparison\tfloor\tbins_tested\tnull_scores\tbins_bh_significant\tmin_pval\n")
        for r in rows:
            fh.write("%s\t%g\t%d\t%d\t%d\t%.6g\n" % r)
    print("\nWrote %s" % args.out, file=sys.stderr)

    print("\n=== SUMMARY: bins BH-significant by floor ===", file=sys.stderr)
    cmps = sorted(set(r[0] for r in rows))
    hdr = "%-14s" % "comparison" + "".join("%12s" % ("floor " + ("%g" % f)) for f in floors)
    print(hdr, file=sys.stderr)
    for c in cmps:
        line = "%-14s" % c
        for f in floors:
            m = [r for r in rows if r[0] == c and abs(r[1] - f) < 1e-9]
            line += "%12s" % (format(m[0][4], ",") if m else "-")
        print(line, file=sys.stderr)
    print("\nRead across each row: if the calls are stable across floors, the "
          "0.5 default is not driving the results. If they are not, the floor "
          "is a real analytic choice and needs justifying in Methods.",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
