#!/usr/bin/env python3
# baseline_comparison.py  (v3 -- reads bigWigs; unpacks (scores, n_masked) from the scorer)
#
# Reviewer question (NAR): "is BEARING simply rediscovering what a simpler
# statistic would show?"  For the SAME 200 bp bins and the SAME min_signal mask,
# compute the BEARING score and the simple baselines a typical analyst would
# reach for, then report rank agreement and overlap of the top-scoring bins.
#
# Per-bin scalar statistics (each summarizes the 6-track vector at that bin):
#   bearing_kl     clamped per-track KL, summed   <- the production score
#                  score_i = max(P_i * log2(P_i/Q_i), 0)
#   bearing_jsd    the pipeline's bounded JSD variant (score_method="jsd")
#   unclamped_kl   full KL(P||Q) = sum_i P_i log2(P_i/Q_i)   (no clamp)
#   jsd_full       full Jensen-Shannon divergence JSD(P||Q)
#   neg_entropy    -H(P): composition "peakedness"
#   sum_z          sum_i z_i of raw signal (z per track, genome-wide mean/sd)
#   max_z          max_i |z_i|
#
# P and Q are reconstructed from the BigWigs (200 bp bin means, abs() on the
# negative-strand track, pseudocount, row-normalize), mirroring bigwig_to_qcat.
# bearing_kl / bearing_jsd call the PRODUCTION scorer (kl_scores_per_bin), so
# the BEARING side is the pipeline's own code, not a paraphrase.
#
# CAVEAT: this reconstruction does not apply blacklist masking or per-track
# noise floors. Every statistic (BEARING included) is computed from the same
# reconstructed matrix, so the COMPARISON is internally consistent. Use
# --check-qcat to verify the reconstruction against the stored qcat scores on
# one chromosome before quoting any number.
#
#   conda activate bearing
#   python baseline_comparison.py --sheet workflow/config/samples.tsv \
#       --sample-name DN_rep1 --chrom-sizes workflow/resources/mm10.chrom.sizes \
#       --check-qcat workflow/results/DN_rep1.qcat.bgz --check-chrom chr19
#
#   # differential contrast (the one the paper makes):
#   python baseline_comparison.py --sheet workflow/config/samples.tsv \
#       --sample-name DN_rep1 --vs ProB_rep1 \
#       --chrom-sizes workflow/resources/mm10.chrom.sizes
#
# Tip: start with --chroms chr6,chr12 (Tcrb + Igh) for a fast first look.
# ASCII only. Reads real data; fabricates nothing.

import argparse
import gzip
import os
import re
import sys

import numpy as np

SCORE_RE = re.compile(r"\[\s*(-?[0-9.eE+]+)\s*,\s*\d+\s*\]")
BIN_SIZE = 200
# panel order from the sheet: ATAC,RNAseq+,RNAseq-,CTCF,Cohesin,H3K27ac
NEG_STRAND_IDX = 2


def import_production_scorer(repo_root):
    sys.path.insert(0, repo_root)
    try:
        from bigwig_to_qcat import kl_scores_per_bin, signals_to_prob
    except Exception as e:
        sys.exit("could not import bigwig_to_qcat from %s: %s\n"
                 "run from the repo root or pass --repo" % (repo_root, e))
    return kl_scores_per_bin, signals_to_prob


def read_sheet(path, name):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0].strip() in ("", "sample"):
                continue
            if f[0].strip() == name:
                return [p.strip() for p in f[4].split(",") if p.strip()]
    sys.exit("sample %s not found in %s" % (name, path))


def read_chrom_sizes(path, chroms):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 2:
                continue
            c, n = f[0], int(f[1])
            if chroms:
                if c in chroms:
                    out.append((c, n))
                continue
            if "_" in c or c == "chrM":
                continue          # skip scaffolds/mito by default
            out.append((c, n))
    if not out:
        sys.exit("no chromosomes selected from %s" % path)
    return out


def load_matrix(bw_paths, chrom_sizes, data_dir):
    """Return (raw float32 matrix [bins x tracks], list of (chrom, nbins))."""
    try:
        import pyBigWig
    except ImportError:
        sys.exit("pyBigWig is required (conda activate bearing)")

    handles = []
    for p in bw_paths:
        full = p if os.path.isabs(p) else os.path.normpath(os.path.join(data_dir, p))
        if not os.path.exists(full):
            sys.exit("bigWig not found: %s (adjust --data-dir)" % full)
        handles.append(pyBigWig.open(full))

    blocks, index = [], []
    for chrom, clen in chrom_sizes:
        nb = clen // BIN_SIZE
        if nb < 1:
            continue
        cols = []
        for bw in handles:
            if chrom not in bw.chroms():
                cols.append(np.zeros(nb, dtype=np.float32))
                continue
            v = bw.stats(chrom, 0, nb * BIN_SIZE, type="mean", nBins=nb)
            cols.append(np.array([0.0 if x is None else x for x in v],
                                 dtype=np.float32))
        blocks.append(np.vstack(cols).T)
        index.append((chrom, nb))
        print("    %-6s %9d bins" % (chrom, nb))
    for bw in handles:
        bw.close()

    R = np.vstack(blocks)
    R[~np.isfinite(R)] = 0.0
    R[:, NEG_STRAND_IDX] = np.abs(R[:, NEG_STRAND_IDX])   # negative-strand track
    R[R < 0] = 0.0
    return R, index


def statistics(R, kl_fn, prob_fn, Q=None):
    P = prob_fn(R.astype(np.float64))
    if Q is None:
        Q = P.mean(axis=0)
    logPQ = np.log2(P / Q[None, :])

    def _scores(method):
        # kl_scores_per_bin returns (scores, n_masked) -- the docstring says
        # otherwise, but line ~1755 of bigwig_to_qcat.py returns a tuple.
        res = kl_fn(P, Q, raw_signal_matrix=R, min_signal=0.0, score_method=method)
        arr = res[0] if isinstance(res, tuple) else res
        return np.asarray(arr)

    out = {}
    out["bearing_kl"] = _scores("kl").sum(axis=1)
    try:
        out["bearing_jsd"] = _scores("jsd").sum(axis=1)
    except Exception as e:
        print("  (bearing_jsd unavailable: %s)" % e, file=sys.stderr)

    out["unclamped_kl"] = (P * logPQ).sum(axis=1)

    M = 0.5 * (P + Q[None, :])
    out["jsd_full"] = 0.5 * ((P * np.log2(P / M)).sum(axis=1)
                             + (Q[None, :] * np.log2(Q[None, :] / M)).sum(axis=1))

    out["neg_entropy"] = (P * np.log2(P)).sum(axis=1)

    mu, sd = R.mean(axis=0), R.std(axis=0)
    sd = np.where(sd == 0, 1.0, sd)
    Z = (R - mu[None, :]) / sd[None, :]
    out["sum_z"] = Z.sum(axis=1)
    out["max_z"] = np.abs(Z).max(axis=1)
    return out, Q


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(np.float64)
    rb = np.argsort(np.argsort(b)).astype(np.float64)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


def top_jaccard(a, b, frac):
    k = max(1, int(round(frac * a.size)))
    ta = set(np.argpartition(-a, k - 1)[:k].tolist())
    tb = set(np.argpartition(-b, k - 1)[:k].tolist())
    return len(ta & tb) / float(len(ta | tb))


def report(stats, ref, frac):
    print("\n%-14s %10s %18s" % ("baseline", "Spearman", "top-%g%% Jaccard" % (100 * frac)))
    for k in stats:
        if k == ref:
            continue
        print("%-14s %+10.3f %18.3f"
              % (k, spearman(stats[ref], stats[k]), top_jaccard(stats[ref], stats[k], frac)))
    print("\nReading the table:")
    print("  rho ~ 1 AND Jaccard ~ 1  -> that baseline reproduces the BEARING ranking;")
    print("                              BEARING adds nothing on this axis.")
    print("  high rho, low Jaccard    -> agree globally, disagree on exactly the bins")
    print("                              that get called (which is what matters).")
    print("  None of these baselines ATTRIBUTES the change to a track. That, plus the")
    print("  permutation FDR, is BEARING's claim -- not raw per-bin detection.")


def qcat_scores_for_chrom(path, chrom):
    tot = []
    seen = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                if seen:
                    break          # qcat is coordinate-sorted; past our chrom
                continue
            seen = True
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            sc = SCORE_RE.findall(p[3])
            tot.append(sum(float(s) for s in sc) if sc else 0.0)
    return np.asarray(tot)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", default="workflow/config/samples.tsv")
    ap.add_argument("--sample-name", required=True)
    ap.add_argument("--vs", default=None, help="second sample -> compare differentials")
    ap.add_argument("--chrom-sizes", default="workflow/resources/mm10.chrom.sizes")
    ap.add_argument("--chroms", default="", help="comma list; default = primary chroms")
    ap.add_argument("--data-dir", default=".",
                    help="directory the sheet's relative bw paths resolve against")
    ap.add_argument("--repo", default=".")
    ap.add_argument("--min-signal", type=float, default=0.1)
    ap.add_argument("--stride", type=int, default=1,
                    help="use every Nth passing bin for the correlations (Q uses all)")
    ap.add_argument("--top-frac", type=float, default=0.01)
    ap.add_argument("--check-qcat", default=None,
                    help="qcat.bgz to validate the reconstruction against")
    ap.add_argument("--check-chrom", default="chr19")
    a = ap.parse_args()

    kl_fn, prob_fn = import_production_scorer(os.path.abspath(a.repo))
    chroms = [c for c in a.chroms.split(",") if c]
    cs = read_chrom_sizes(a.chrom_sizes, chroms)

    print("loading bigWigs for %s" % a.sample_name)
    bw1 = read_sheet(a.sheet, a.sample_name)
    R1, index = load_matrix(bw1, cs, a.data_dir)
    print("  matrix: %d bins x %d tracks" % R1.shape)

    keep = R1.sum(axis=1) >= a.min_signal
    print("  min_signal=%.3g -> %d bins pass (%.1f%% masked)"
          % (a.min_signal, int(keep.sum()), 100.0 * (1 - keep.mean())))
    if keep.sum() == 0:
        sys.exit("no bins pass the mask")

    st1, Q1 = statistics(R1[keep], kl_fn, prob_fn)

    if a.check_qcat:
        print("\ncross-check vs stored qcat scores on %s" % a.check_chrom)
        stored = qcat_scores_for_chrom(a.check_qcat, a.check_chrom)
        off, mine = 0, None
        for chrom, nb in index:
            if chrom == a.check_chrom:
                s, _ = statistics(R1[off:off + nb], kl_fn, prob_fn, Q=Q1)
                mine = s["bearing_kl"]
                break
            off += nb
        if mine is None or stored.size == 0:
            print("  could not align (%s absent from selection/qcat); skipped"
                  % a.check_chrom)
        else:
            n = min(mine.size, stored.size)
            rho = spearman(mine[:n], stored[:n])
            print("  n=%d  Spearman(reconstructed, stored) = %.4f" % (n, rho))
            if rho < 0.90:
                print("  WARNING: reconstruction differs from production (blacklist /")
                print("  noise floors / clip not applied here). The baseline table is")
                print("  still internally consistent, but do NOT quote it as production.")
            else:
                print("  reconstruction is faithful; baseline table is trustworthy.")

    if a.stride > 1:
        st1 = {k: v[::a.stride] for k, v in st1.items()}
        print("  correlations on every %dth passing bin (n=%d)"
              % (a.stride, st1["bearing_kl"].size))

    if a.vs is None:
        report(st1, "bearing_kl", a.top_frac)
        return

    print("\nloading bigWigs for %s" % a.vs)
    bw2 = read_sheet(a.sheet, a.vs)
    R2, _ = load_matrix(bw2, cs, a.data_dir)
    if R2.shape != R1.shape:
        sys.exit("bin grids differ (%s vs %s)" % (R1.shape, R2.shape))
    keep2 = R2.sum(axis=1) >= a.min_signal
    both = keep & keep2
    print("  bins passing in BOTH samples: %d" % int(both.sum()))
    if both.sum() == 0:
        sys.exit("no shared passing bins")

    # each sample gets its own Q, as in production
    s1, _ = statistics(R1[both], kl_fn, prob_fn)
    s2, _ = statistics(R2[both], kl_fn, prob_fn)
    diff = {k: np.abs(s1[k] - s2[k]) for k in s1 if k in s2}
    if a.stride > 1:
        diff = {k: v[::a.stride] for k, v in diff.items()}
    print("\nDIFFERENTIAL |%s - %s| (n=%d bins):"
          % (a.sample_name, a.vs, diff["bearing_kl"].size))
    report(diff, "bearing_kl", a.top_frac)


if __name__ == "__main__":
    main()
