#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : dev/replicate_swap_null.py
# Version  : 1.0.2
# Date     : 2026-09-26
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# Status   : dev/ prototype (category B until a manuscript number uses it)
# Changes  : 1.0.1 --regions also reads regions_manuscript.tsv (name, chr:start-end);
#            relative qcat paths resolve from cwd (workflow/) or sheet dir
#            1.0.2 perbin adds the four per-sample summed scores (a1 a2 b1 b2)
#            so replicate-driven bins can be inspected
# ----------------------------------------------------------------------
"""
replicate_swap_null.py
======================
Replicate-aware differential test for BEARING with n = 2 per condition,
using a replicate label-swap null instead of the circular-shift null.

WHY
---
Production (compare_qcat.py -> bearing_pvalue.py --diff) tests the summed
differential

    d(b) = mean(S_A1, S_A2) - mean(S_B1, S_B2),   S = summed per-track KL,

against a null built by circularly shifting EVERY track of EVERY sample
independently (generate_perm_nulls.py seeds per sample and file). In that
null, rep1 and rep2 of a condition are no longer aligned, so the null
differential is a difference of four unrelated samples. Real replicates are
correlated, so this null is wider than true replicate noise. That is one
candidate explanation for the CONSERVATIVE calibration (lambda 0.34-0.44,
zero BH-significant bins). Replicate disagreement never enters the test.

With two replicates per condition there are exactly two relabelings that
mix one replicate of each condition into each group:

    swap1: (A1 + B1)/2 - (A2 + B2)/2  =  [(A1 - A2) + (B1 - B2)] / 2
    swap2: (A1 + B2)/2 - (A2 + B1)/2  =  [(A1 - A2) - (B1 - B2)] / 2

The true A-vs-B effect cancels exactly in both, and each has the same
variance as the noise part of d(b) (four independent replicate errors, each
weighted 1/2). Pooled over bins, swap1 and swap2 are an empirical null for
|d| that keeps the real replicate correlation and real local chromatin
structure. Replicate disagreement enters the test through this null: a
bin is called only if its between-condition difference exceeds what
replicate-vs-replicate differences produce at comparable signal.

HETEROSCEDASTICITY (stratified null)
------------------------------------
Replicate noise grows with signal. With --strata K the null is pooled
within K quantile strata of the swap-invariant covariate
c(b) = mean(S_A1, S_A2, S_B1, S_B2). c does not change under relabeling,
so stratifying on it does not use the condition labels (same logic as
edgeR filterByExpr / IHW covariates). --strata 1 = one genome-wide pool.

BATCH DIAGNOSTIC
----------------
If rep1 of both conditions shares a batch (e.g. processed together),
a batch offset cancels in d(b) and in swap2 but doubles in swap1. The
script reports the spread of swap1 and swap2 separately; swap1 much wider
than swap2 means a replicate-batch effect, and --swap-mode swap2 then
gives the matched null.

WHAT A PER-BIN t-STATISTIC CANNOT DO HERE (documented dead end)
----------------------------------------------------------------
A moderated t = d / (SE + s0) with SE from the two within-condition
differences has no valid null at n = 2: under the swaps the numerator is
built from the same two differences as the denominator (t1^2 + t2^2 = 2
when s0 = 0), and the circular-shift null inflates the within-condition
SE because replicates are shifted independently. Stratifying the swap null
on c is the valid replicate-aware alternative.

FLOOR
-----
--min-signal applies the production differential floor to |d| for BOTH the
observed and the null values (bearing_pvalue.py --diff --min-signal, config
pvalue_min_signal = 0.5), so p = P(|d'| >= |d| given |d'| >= floor), as in
production.

FIDELITY
--------
Per-sample qcats are read with compare_qcat.parse_qcat_bgz (production
parser). Bins are the intersection across the loaded samples (production
intersects across ALL sheet samples, so a few bins can differ; use
--prod-stats to restrict the observed test set to production's exact tested
bins). p-values use bearing_pvalue.empirical_pvals and BH uses
bearing_pvalue.bh_fdr. --check-diff-qcat verifies that the reconstructed
d(b) matches the production diff qcat before any number is quoted.

INPUTS
------
--sheet       results/samples.qcat.tsv (sample, condition, replicate, qcat)
--cond-a/-b   condition names, exactly two replicates each
--chroms      optional restriction for a fast first look (e.g. chr6,chr12);
              the null is then pooled from those chromosomes only
--prod-stats  optional production diff_<A>_vs_<B>.stats.tsv (adds
              production pval / pval_adj_bh columns and restricts the
              observed tested set to production's tested bins)
--prod-null-qcats  optional production perm diff qcats (a few are enough)
              to compare the width of the circular-shift null with the
              swap null
--regions     optional BED (chrom start end name) for per-region counts

OUTPUTS (--out-prefix P)
------------------------
P.summary.tsv     key / value: bins, floor, tested, BH hits (swap and
                  production), pi0 estimates, batch diagnostic
P.quantiles.tsv   |d| quantiles: observed, swap1, swap2, production null
P.perbin.tsv.gz   tested bins: d, cov, stratum, p_swap, q_swap
                  (+ p_prod, q_prod)
P.regions.tsv     per region (if --regions)
P.pdf             survival curves of |d| + p-value histograms

USAGE
-----
  conda activate bearing
  cd /mnt/isilon/bassing_lab/integration_paper/bearing
  # 1. fidelity check + fast look on Tcrb/Igh chromosomes
  python dev/replicate_swap_null.py \\
      --sheet workflow/results/samples.qcat.tsv \\
      --cond-a DN --cond-b DP --chroms chr6,chr12 \\
      --min-signal 0.5 --strata 10 \\
      --check-diff-qcat workflow/results/compare/diff_DN_vs_DP.qcat.bgz \\
      --prod-stats workflow/results/pvalue/diff_DN_vs_DP.stats.tsv \\
      --prod-null-qcats workflow/results/perm/perm{1,2,3}/diff_comparison/diff_DN_vs_DP.qcat.bgz \\
      --regions regions_manuscript.tsv \\
      --out-prefix swapnull/DN_vs_DP_chr6_chr12
  # 2. genome-wide (drop --chroms), one job per comparison

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import csv
import gzip
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def import_production(repo):
    sys.path.insert(0, repo)
    try:
        import compare_qcat as cq
        import bearing_pvalue as bp
    except Exception as e:
        sys.exit("could not import compare_qcat / bearing_pvalue from %s: %s"
                 % (repo, e))
    for mod, fn in ((cq, "parse_qcat_bgz"), (bp, "empirical_pvals"),
                    (bp, "bh_fdr"), (bp, "collect_null_scores")):
        if not hasattr(mod, fn):
            sys.exit("production interface drifted: %s.%s missing"
                     % (mod.__name__, fn))
    return cq, bp


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


# ----------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------
def read_sheet(path):
    rows = []
    with open(path) as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    rdr = csv.DictReader(lines, delimiter="\t")
    for r in rdr:
        rows.append({k.strip(): (v or "").strip() for k, v in r.items() if k})
    need = {"sample", "condition", "qcat"}
    if not rows or not need.issubset(rows[0].keys()):
        sys.exit("sheet %s needs columns %s" % (path, sorted(need)))
    # Relative qcat paths: Snakemake writes them relative to workflow/ (the
    # cwd it runs in), not to the sheet's folder. Try cwd first, then the
    # sheet's folder; fail loudly if neither exists.
    base = os.path.dirname(os.path.abspath(path))
    for r in rows:
        q = r["qcat"]
        if os.path.isabs(q):
            cands = [q]
        else:
            cands = [os.path.abspath(q), os.path.join(base, q)]
        hit = next((c for c in cands if os.path.exists(c)), None)
        if hit is None:
            sys.exit("qcat for %s not found; tried: %s (run from workflow/)"
                     % (r["sample"], " , ".join(cands)))
        r["qcat"] = hit
    return rows


def pick_replicates(rows, cond):
    reps = [r for r in rows if r["condition"] == cond]
    if len(reps) != 2:
        sys.exit("condition %s has %d replicate(s) in the sheet; the swap "
                 "null needs exactly 2" % (cond, len(reps)))
    reps.sort(key=lambda r: (int(r.get("replicate") or 0), r["sample"]))
    return reps


def load_summed(cq, path, chroms):
    """Return dict (chrom,start,end) -> summed per-track score, and n_tracks."""
    out = cq.parse_qcat_bgz(path, chroms=chroms)
    bins, scores, n_states = out[0], out[1], out[2]
    s = scores.astype(np.float64).sum(axis=1)
    return bins, s, n_states


def align(bin_lists, value_lists):
    common = set(bin_lists[0])
    for b in bin_lists[1:]:
        common &= set(b)
    common = sorted(common, key=lambda x: (x[0], x[1]))
    mats = []
    for bins, vals in zip(bin_lists, value_lists):
        idx = {b: i for i, b in enumerate(bins)}
        mats.append(np.fromiter((vals[idx[b]] for b in common),
                                dtype=np.float64, count=len(common)))
    return common, mats


# ----------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------
def storey_pi0(p, lam=0.5):
    p = np.asarray(p)
    if len(p) == 0:
        return float("nan")
    return float(min(1.0, np.mean(p > lam) / (1.0 - lam)))


def strata_edges(cov, k):
    if k <= 1:
        return np.array([-np.inf, np.inf])
    qs = np.quantile(cov, np.linspace(0, 1, k + 1))
    qs[0], qs[-1] = -np.inf, np.inf
    return np.unique(qs)


def stratified_pvals(bp, obs_abs, obs_stratum, null_abs, null_stratum,
                     n_strata):
    """Empirical p per stratum with the production estimator."""
    p = np.ones(len(obs_abs), dtype=np.float64)
    null_n = {}
    for k in range(n_strata):
        om = obs_stratum == k
        if not om.any():
            continue
        nv = np.sort(null_abs[null_stratum == k])
        null_n[k] = len(nv)
        if len(nv) == 0:
            log("  WARNING: stratum %d has no null values above the floor; "
                "its %d bins get p = 1" % (k, int(om.sum())))
            continue
        p[om] = bp.empirical_pvals(obs_abs[om], nv)
    return p, null_n


def quantile_row(label, v, qs):
    v = np.asarray(v)
    row = {"distribution": label, "n": int(len(v))}
    for q in qs:
        row["q%s" % q] = float(np.quantile(v, q / 100.0)) if len(v) else float("nan")
    return row


# ----------------------------------------------------------------------
# Optional production joins
# ----------------------------------------------------------------------
def load_prod_stats(path, chroms):
    import pandas as pd
    use = ["chrom", "start", "end", "bearing_score", "pval", "pval_adj_bh"]
    out = []
    for ch in pd.read_csv(path, sep="\t", usecols=lambda c: c in use,
                          chunksize=2_000_000):
        if chroms:
            ch = ch[ch["chrom"].isin(chroms)]
        out.append(ch)
    df = pd.concat(out, ignore_index=True)
    missing = [c for c in use if c not in df.columns]
    if missing:
        sys.exit("prod stats %s lacks columns %s" % (path, missing))
    return df


def check_diff(cq, path, chroms, common, d):
    out = cq.parse_qcat_bgz(path, chroms=chroms)
    bins, scores = out[0], out[1]
    prod = dict(zip(bins, scores.astype(np.float64).sum(axis=1)))
    idx = [i for i, b in enumerate(common) if b in prod]
    if not idx:
        sys.exit("check-diff-qcat: no overlapping bins")
    ours = d[idx]
    theirs = np.array([prod[common[i]] for i in idx])
    err = np.abs(ours - theirs)
    tol = 1e-4 * np.maximum(1.0, np.abs(theirs))
    return {
        "check_overlap_bins": len(idx),
        "check_prod_only_bins": len(prod) - len(idx),
        "check_ours_only_bins": len(common) - len(idx),
        "check_max_abs_err": float(err.max()),
        "check_frac_within_tol": float(np.mean(err <= tol)),
    }


def read_regions(path):
    regs = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#") or ln.startswith("track"):
                continue
            f = ln.rstrip("\n").split("\t")
            # Format 2: regions_manuscript.tsv  (name, chrom:start-end, ...)
            if len(f) >= 2 and ":" in f[1] and "-" in f[1]:
                try:
                    c, se = f[1].split(":", 1)
                    s, e = [int(x.replace(",", "")) for x in se.split("-", 1)]
                except ValueError:
                    continue
                regs.append((c, s, e, f[0]))
                continue
            # Format 1: BED (chrom, start, end, [name])
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue          # header line
            regs.append((f[0], s, e, f[3] if len(f) > 3 else "%s:%d-%d" % (f[0], s, e)))
    if not regs:
        sys.exit("no regions parsed from %s (need BED or name<TAB>chr:start-end)"
                 % path)
    return regs


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--cond-a", required=True)
    ap.add_argument("--cond-b", required=True)
    ap.add_argument("--chroms", default=None, help="comma list; default all")
    ap.add_argument("--min-signal", type=float, default=0.5,
                    help="floor on |d| for observed AND null (production 0.5)")
    ap.add_argument("--strata", type=int, default=10,
                    help="quantile strata of swap-invariant covariate; 1 = none")
    ap.add_argument("--swap-mode", default="both",
                    choices=["both", "swap1", "swap2"])
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--check-diff-qcat", default=None)
    ap.add_argument("--prod-stats", default=None)
    ap.add_argument("--prod-null-qcats", nargs="+", default=None)
    ap.add_argument("--regions", default=None)
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cq, bp = import_production(args.repo)
    chroms = [c.strip() for c in args.chroms.split(",")] if args.chroms else None
    od = os.path.dirname(os.path.abspath(args.out_prefix))
    os.makedirs(od, exist_ok=True)

    rows = read_sheet(args.sheet)
    A = pick_replicates(rows, args.cond_a)
    B = pick_replicates(rows, args.cond_b)
    samples = A + B
    log("A: %s   B: %s" % ([r["sample"] for r in A], [r["sample"] for r in B]))

    bl, vl = [], []
    ntr = set()
    for r in samples:
        log("loading %s (%s)" % (r["sample"], r["qcat"]))
        b, s, n = load_summed(cq, r["qcat"], chroms)
        bl.append(b); vl.append(s); ntr.add(n)
        log("  %s bins, %d tracks" % (format(len(b), ","), n))
    if len(ntr) != 1:
        sys.exit("samples disagree on track count: %s" % sorted(ntr))
    common, (a1, a2, b1, b2) = align(bl, vl)
    del bl, vl
    n_common = len(common)
    log("common bins: %s" % format(n_common, ","))

    d = (a1 + a2) / 2.0 - (b1 + b2) / 2.0
    sw1 = (a1 + b1) / 2.0 - (a2 + b2) / 2.0
    sw2 = (a1 + b2) / 2.0 - (a2 + b1) / 2.0
    cov = (a1 + a2 + b1 + b2) / 4.0
    summary = {"cond_a": args.cond_a, "cond_b": args.cond_b,
               "samples": ",".join(r["sample"] for r in samples),
               "chroms": args.chroms or "all", "common_bins": n_common,
               "min_signal": args.min_signal, "strata": args.strata,
               "swap_mode": args.swap_mode}

    if args.check_diff_qcat:
        log("checking reconstruction against %s" % args.check_diff_qcat)
        chk = check_diff(cq, args.check_diff_qcat, chroms, common, d)
        summary.update(chk)
        log("  max |ours - prod| = %.3g, within tol: %.4f"
            % (chk["check_max_abs_err"], chk["check_frac_within_tol"]))
        if chk["check_frac_within_tol"] < 0.999:
            log("  WARNING: reconstruction does not match production diff; "
                "do not quote numbers until this is resolved")

    # strata on the swap-invariant covariate (all common bins)
    edges = strata_edges(cov, args.strata)
    stratum = np.clip(np.searchsorted(edges, cov, side="right") - 1,
                      0, len(edges) - 2)
    n_strata = len(edges) - 1

    floor = args.min_signal
    obs_abs = np.abs(d)
    tested = obs_abs >= floor

    # optional: restrict observed test set to production tested bins
    prod_p = prod_q = None
    if args.prod_stats:
        log("loading production stats %s" % args.prod_stats)
        ps = load_prod_stats(args.prod_stats, chroms)
        key = {(c, int(s), int(e)): i for i, (c, s, e) in
               enumerate(zip(ps["chrom"], ps["start"], ps["end"]))}
        pi = np.fromiter((key.get(b, -1) for b in common), dtype=np.int64,
                         count=n_common)
        in_prod = pi >= 0
        summary["prod_tested_bins"] = int(len(ps))
        summary["prod_tested_in_common"] = int(in_prod.sum())
        summary["ours_tested_not_in_prod"] = int((tested & ~in_prod).sum())
        tested = tested & in_prod
        prod_p = np.full(n_common, np.nan)
        prod_q = np.full(n_common, np.nan)
        prod_p[in_prod] = ps["pval"].to_numpy()[pi[in_prod]]
        prod_q[in_prod] = ps["pval_adj_bh"].to_numpy()[pi[in_prod]]

    # null pool
    parts_abs, parts_str = [], []
    if args.swap_mode in ("both", "swap1"):
        m = np.abs(sw1) >= floor
        parts_abs.append(np.abs(sw1)[m]); parts_str.append(stratum[m])
    if args.swap_mode in ("both", "swap2"):
        m = np.abs(sw2) >= floor
        parts_abs.append(np.abs(sw2)[m]); parts_str.append(stratum[m])
    null_abs = np.concatenate(parts_abs)
    null_str = np.concatenate(parts_str)
    summary["null_values_above_floor"] = int(len(null_abs))
    summary["min_achievable_p_unstratified"] = 1.0 / (len(null_abs) + 1)

    ti = np.where(tested)[0]
    summary["tested_bins"] = int(len(ti))
    log("tested bins: %s   null values: %s"
        % (format(len(ti), ","), format(len(null_abs), ",")))
    p, null_n = stratified_pvals(bp, obs_abs[ti], stratum[ti],
                                 null_abs, null_str, n_strata)
    if null_n:
        summary["min_null_per_stratum"] = int(min(null_n.values()))
        # Resolution limit: the smallest p a bin can get in the smallest
        # stratum, and how many bins must reach it before BH can reject it.
        pmin = 1.0 / (min(null_n.values()) + 1)
        summary["min_achievable_p_stratified"] = pmin
        summary["bh_rank_needed_at_min_p"] = int(np.ceil(pmin * max(1, len(ti))
                                                         / args.fdr))
        log("  resolution: min p %.2e in smallest stratum; BH can reject it "
            "only if >= %d bins reach it (fewer strata = finer resolution)"
            % (pmin, summary["bh_rank_needed_at_min_p"]))
    rej, q = bp.bh_fdr(p, args.fdr)
    summary["bh_hits_swap"] = int(np.sum(q < args.fdr))
    summary["bh_hits_swap_DNdir"] = int(np.sum((q < args.fdr) & (d[ti] > 0)))
    summary["pi0_swap"] = storey_pi0(p)
    summary["frac_p_lt_0.05_swap"] = float(np.mean(p < 0.05)) if len(p) else float("nan")
    if prod_p is not None:
        pp = prod_p[ti]; pq = prod_q[ti]
        summary["bh_hits_prod_same_bins"] = int(np.sum(pq < args.fdr))
        summary["pi0_prod_same_bins"] = storey_pi0(pp[np.isfinite(pp)])
        summary["frac_p_lt_0.05_prod"] = float(np.mean(pp < 0.05))
        both = np.isfinite(pp)
        if both.sum() > 2:
            from scipy.stats import spearmanr
            summary["spearman_p_swap_vs_prod"] = float(
                spearmanr(p[both], pp[both]).correlation)

    # batch diagnostic
    a_sw1 = np.abs(sw1); a_sw2 = np.abs(sw2)
    for qv in (0.9, 0.99, 0.999):
        v1 = float(np.quantile(a_sw1, qv)); v2 = float(np.quantile(a_sw2, qv))
        summary["swap1_over_swap2_q%s" % qv] = (v1 / v2) if v2 > 0 else float("nan")

    # quantiles
    QS = (50, 90, 99, 99.9, 99.99)
    qrows = [quantile_row("observed_abs_d_all", obs_abs, QS),
             quantile_row("swap1_abs_all", a_sw1, QS),
             quantile_row("swap2_abs_all", a_sw2, QS),
             quantile_row("observed_abs_d_above_floor", obs_abs[obs_abs >= floor], QS),
             quantile_row("swap_null_above_floor", null_abs, QS)]
    prod_null = None
    if args.prod_null_qcats:
        log("loading %d production null qcat(s) (genome-wide files)"
            % len(args.prod_null_qcats))
        res = bp.collect_null_scores(args.prod_null_qcats,
                                     min_signal=floor, diff_mode=True)
        prod_null = np.asarray(res[0] if isinstance(res, tuple) else res)
        qrows.append(quantile_row("prod_shift_null_above_floor", prod_null, QS))
        summary["prod_null_values_above_floor"] = int(len(prod_null))
        if chroms:
            summary["note_prod_null"] = ("production null is genome-wide; "
                                         "swap null is restricted to --chroms")

    # ---- write outputs ----
    P = args.out_prefix
    with open(P + ".summary.tsv", "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in summary.items():
            fh.write("%s\t%s\n" % (k, v))
    with open(P + ".quantiles.tsv", "w") as fh:
        cols = ["distribution", "n"] + ["q%s" % q for q in QS]
        fh.write("\t".join(cols) + "\n")
        for r in qrows:
            fh.write("\t".join(str(r[c]) if c in ("distribution", "n")
                               else "%.6g" % r[c] for c in cols) + "\n")
    with gzip.open(P + ".perbin.tsv.gz", "wt") as fh:
        cols = ["chrom", "start", "end", "d", "cov", "stratum", "p_swap", "q_swap",
                "a1", "a2", "b1", "b2"]
        if prod_p is not None:
            cols += ["p_prod", "q_prod"]
        fh.write("\t".join(cols) + "\n")
        for j, i in enumerate(ti):
            c, s, e = common[i]
            row = [c, str(s), str(e), "%.6g" % d[i], "%.6g" % cov[i],
                   str(int(stratum[i])), "%.4g" % p[j], "%.4g" % q[j],
                   "%.4g" % a1[i], "%.4g" % a2[i], "%.4g" % b1[i], "%.4g" % b2[i]]
            if prod_p is not None:
                row += ["%.4g" % prod_p[i], "%.4g" % prod_q[i]]
            fh.write("\t".join(row) + "\n")

    if args.regions:
        regs = read_regions(args.regions)
        tb = [common[i] for i in ti]
        tchrom = np.array([b[0] for b in tb])
        tstart = np.array([b[1] for b in tb], dtype=np.int64)
        with open(P + ".regions.tsv", "w") as fh:
            cols = ["region", "chrom", "start", "end", "tested_bins",
                    "hits_swap", "hits_swap_Adir", "min_q_swap"]
            if prod_p is not None:
                cols += ["hits_prod", "min_q_prod"]
            fh.write("\t".join(cols) + "\n")
            for (c, s, e, nm) in regs:
                m = (tchrom == c) & (tstart >= s) & (tstart < e)
                n = int(m.sum())
                row = [nm, c, str(s), str(e), str(n),
                       str(int(np.sum(q[m] < args.fdr))),
                       str(int(np.sum((q[m] < args.fdr) & (d[ti][m] > 0)))),
                       ("%.3g" % q[m].min()) if n else "NA"]
                if prod_p is not None:
                    pq = prod_q[ti][m]
                    row += [str(int(np.sum(pq < args.fdr))),
                            ("%.3g" % np.nanmin(pq)) if n else "NA"]
                fh.write("\t".join(row) + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))

        def surv(v, lab, col, ls="-"):
            v = np.sort(np.asarray(v))
            if len(v) == 0:
                return
            y = 1.0 - np.arange(len(v)) / float(len(v))
            step = max(1, len(v) // 20000)
            ax[0].plot(v[::step], y[::step], color=col, ls=ls, lw=1.4, label=lab)
        surv(obs_abs[obs_abs >= floor], "observed |d|", "#222222")
        surv(null_abs, "swap null (replicate)", "#1f77b4")
        if prod_null is not None:
            surv(prod_null, "circular-shift null (production)", "#d62728", "--")
        ax[0].set_yscale("log")
        ax[0].set_xlabel("|summed differential|  (>= floor %.2g)" % floor)
        ax[0].set_ylabel("fraction >= x")
        ax[0].set_title("%s vs %s: null width" % (args.cond_a, args.cond_b),
                        fontsize=10)
        ax[0].legend(fontsize=8)
        bins = np.linspace(0, 1, 41)
        ax[1].hist(p, bins=bins, color="#1f77b4", alpha=0.6,
                   label="swap null", density=True)
        if prod_p is not None:
            pp = prod_p[ti]
            ax[1].hist(pp[np.isfinite(pp)], bins=bins, color="#d62728",
                       alpha=0.45, label="production", density=True)
        ax[1].axhline(1.0, color="gray", lw=0.8, ls=":")
        ax[1].set_xlabel("p-value (tested bins)")
        ax[1].set_ylabel("density")
        ax[1].set_title("BH hits: swap %d%s" % (
            summary["bh_hits_swap"],
            ("  production %d" % summary["bh_hits_prod_same_bins"])
            if prod_p is not None else ""), fontsize=10)
        ax[1].legend(fontsize=8)
        for a in ax:
            for sp in ("top", "right"):
                a.spines[sp].set_visible(False)
        plt.tight_layout()
        plt.savefig(P + ".pdf", bbox_inches="tight")
        plt.close()
    except Exception as e:
        log("plot skipped: %s" % e)

    for ext in (".summary.tsv", ".quantiles.tsv", ".perbin.tsv.gz",
                ".regions.tsv" if args.regions else None, ".pdf"):
        if ext:
            log("wrote " + P + ext)


if __name__ == "__main__":
    main()
