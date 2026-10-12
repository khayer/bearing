#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : benchmark/complementarity_benchmark.py
# Version  : 1.0.1
# Date     : 2026-10-12 (L1 is now bearing_pvalue.py --diff-stat l1)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""
Two-condition complementarity benchmark: BEARING vs a per-assay test.

QUESTION
--------
The Discussion (and new Table S5) state that BEARING and per-assay
differential tests (edgeR / csaw / limma on one assay at a time) answer
different questions: BEARING scores a change in the COMPOSITION (balance) of
the assays at a bin, a per-assay test scores a change in the MAGNITUDE of one
assay. This script tests that statement on simulated data with a known truth,
instead of asserting it.

SIMULATED DATA
--------------
One synthetic chromosome, six tracks in the production order
[ATAC, RNA+, RNA-, CTCF, RAD21, H3K27ac], two conditions A and B with
--reps replicates each. Shared across all samples (genomic structure):
  * a per-bin coverage factor C_b ~ log-normal(0, --coverage-sigma)
  * a per-bin background composition Pbg_b ~ Dirichlet(--bg-kappa * Q), so
    null bins are not all exactly Q (set --bg-kappa 0 for P = Q)
Per sample: a library factor ~ log-normal(0, --lib-sigma), negative-binomial
noise (--nb-dispersion) on a count scale of --counts-per-unit counts per
signal unit, then Gaussian smoothing (--smooth-sigma-bins). Expected total
signal of a null bin is --depth (BigWig units, the units of min_signal and of
the prior alpha).

Planted blocks (--copies copies of every configuration, --block-bins bins
each, separated by --gap-bins null bins). Feature blocks carry --feature-mult
times the null depth. Composition P = (1 - s) Q + s a, with the allocation
vectors a of simulate_bearing_tracks.py.

  type          what differs between A and B                 truth
  feature_null  nothing (an enriched feature in both)        no change
  magnitude     total signal x f, composition identical      magnitude only
  switch        composition a_A -> a_B, same total           balance only
  rebalance     every track x (1 +/- delta), renormalised,   balance only
                same total (small, distributed shift)
  both          composition a_A -> a_B AND total x f         both

THE TWO TESTS (same simulated BigWig values)
-------------------------------------------
BEARING, as in production (results_prior1_n100), using the production code:
  * per-sample score: bigwig_to_qcat.signals_to_prob, Q = mean of P over all
    bins, bigwig_to_qcat.kl_scores_per_bin with --min-signal and the
    composition prior --prior-strength (alpha)
  * differential: mean(A reps) - mean(B reps) per track (compare_qcat
    diff_qcat); bin statistic = |sum over tracks| (bearing_pvalue --diff);
    tested when that is >= --pvalue-min-signal
  * null: --n-perms rounds of independent circular shifts of every track of
    every sample by U[min_shift, L - min_shift] (shift_bigwig.py), rescored
    and differenced the same way; pooled; bearing_pvalue.empirical_pvals;
    bearing_pvalue.bh_fdr over the tested bins
Per-assay comparator (a limma-trend style moderated t, one test per track):
  * per sample and track, signal is library-normalised (divided by the
    track's genome mean in that sample), log2(x + --assay-pseudo)
  * tested when the mean raw signal of that track across samples is
    >= --assay-min-signal (filterByExpr analogue)
  * moderated t with empirical-Bayes variance shrinkage (Smyth 2004),
    hyperparameters estimated within 20 abundance strata (the trend)
  * bin call: Simes combination of the tested tracks' p-values, BH over
    bins (bin-level FDR); the tracks that moved are those passing BH over
    all bin x track tests (sign pattern below)
  This is a replicate-based test, the way edgeR/limma/csaw are used; it does
  not need the circular-shift null.

OUTPUTS (--outdir)
  per_config.tsv   one row per block configuration: interior bins, fraction
                   called by BEARING and by the per-assay test, BEARING tested
                   fraction, the expected (noise-free) BEARING |differential|,
                   whether BEARING's top track matches the expected top
                   track, and the per-assay sign pattern (concordant = all
                   called tracks move the same way, the signature of a
                   magnitude change)
  per_type.tsv     the same pooled by block type, plus null-bin false-call
                   rates for both tests
  summary.txt      human-readable headline table
  per_bin.tsv.gz   (with --per-bin) every bin with truth and both results

ASCII only. Simulated data only; nothing here is experimental data.
"""

import argparse
import gzip
import math
import os
import sys

import numpy as np

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy import special, stats
except Exception:  # pragma: no cover
    sys.exit("ERROR: scipy is required")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import bigwig_to_qcat as b2q            # noqa: E402  production scorer
import bearing_pvalue as bpv           # noqa: E402  production p-values
from simulate_bearing_tracks import (  # noqa: E402
    DEFAULT_Q, TRACK_NAMES, N_TRACKS, named_allocations, grid_allocations)


# ---------------------------------------------------------------------------
# Block plan
# ---------------------------------------------------------------------------

def _comp(q, a, s):
    p = (1.0 - s) * q + s * np.asarray(a, dtype=np.float64)
    return p / p.sum()


def _rebalance(p, delta):
    w = np.array([1 + delta, 1 - delta] * (N_TRACKS // 2), dtype=np.float64)
    out = p * w
    return out / out.sum()


def block_configs(q):
    """List of configurations
        (config_id, type, P_A, P_B, rel_A, rel_B, s, f, delta)
    P_A / P_B None = that condition has only the per-bin background there
    (depth multiplier 1); otherwise the depth multiplier is
    --feature-mult * rel."""
    named = named_allocations()
    grid = grid_allocations()
    arch, txn = named["architecture"], named["transcription"]
    act = grid["arch_plus_active"]
    ctcf, rad21 = grid["dom_ctcf"], grid["dom_rad21"]
    cfg = []
    for nm, a in (("arch", arch), ("txn", txn)):
        p = _comp(q, a, 0.6)
        cfg.append(("feature_null_%s" % nm, "feature_null", p, p, 1.0, 1.0,
                    0.6, 1.0, 0.0))
    for nm, a in (("arch", arch), ("txn", txn)):
        p = _comp(q, a, 0.6)
        for f in (0.25, 0.5, 2.0, 4.0):
            cfg.append(("magnitude_%s_f%g" % (nm, f), "magnitude", p, p, 1.0, f,
                        0.6, f, 0.0))
    for s in (0.3, 0.6, 0.9):
        cfg.append(("switch_arch_to_txn_s%g" % s, "switch",
                    _comp(q, arch, s), _comp(q, txn, s), 1.0, 1.0, s, 1.0, 0.0))
        cfg.append(("switch_ctcf_to_rad21_s%g" % s, "switch",
                    _comp(q, ctcf, s), _comp(q, rad21, s), 1.0, 1.0, s, 1.0, 0.0))
    base = _comp(q, act, 0.6)
    for d in (0.15, 0.3, 0.5):
        cfg.append(("rebalance_d%g" % d, "rebalance", base, _rebalance(base, d),
                    1.0, 1.0, 0.6, 1.0, d))
    for f in (0.25, 4.0):
        cfg.append(("both_arch_to_txn_f%g" % f, "both",
                    _comp(q, arch, 0.6), _comp(q, txn, 0.6), 1.0, f, 0.6, f, 0.0))
    # on/off: a feature present in one condition only (background in the
    # other), the typical cross-cell-type event; magnitude AND composition.
    # "solo" = a single-assay feature (e.g. strand-specific transcription with
    # little else), the most extreme composition BEARING can see.
    solo = np.zeros(N_TRACKS)
    solo[1] = 1.0
    for nm, a in (("arch", arch), ("txn", txn), ("solo_rnaplus", solo)):
        for s in (0.6, 0.9):
            p = _comp(q, a, s)
            cfg.append(("onoff_gain_%s_s%g" % (nm, s), "onoff", None, p, 1.0, 1.0,
                        s, 1.0, 0.0))
            cfg.append(("onoff_loss_%s_s%g" % (nm, s), "onoff", p, None, 1.0, 1.0,
                        s, 1.0, 0.0))
    # composition switch between two near single-assay features
    cfg.append(("switch_solo_rnaplus_to_ctcf_s0.9", "switch",
                _comp(q, solo, 0.9), _comp(q, np.eye(N_TRACKS)[3], 0.9),
                1.0, 1.0, 0.9, 1.0, 0.0))
    return cfg


def lay_out(cfgs, copies, block_bins, gap_bins, offset, n_bins, rng):
    """Interleave copies (shuffled order) so no configuration sits in one
    region of the chromosome. Returns list of dicts with bin_start/bin_end."""
    order = [(c, k) for k in range(copies) for c in range(len(cfgs))]
    rng.shuffle(order)
    plan, cur = [], offset
    for ci, k in order:
        if cur + block_bins > n_bins:
            raise ValueError("block plan does not fit: increase --length or "
                             "reduce --copies/--block-bins/--gap-bins")
        plan.append(dict(cfg=ci, copy=k, bin_start=cur, bin_end=cur + block_bins))
        cur += block_bins + gap_bins
    return plan


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def simulate(args, rng):
    n = args.length // args.bin_size
    q = DEFAULT_Q / DEFAULT_Q.sum()
    cfgs = block_configs(q)
    plan = lay_out(cfgs, args.copies, args.block_bins, args.gap_bins,
                   args.start_offset_bins, n, rng)

    cov = rng.lognormal(0.0, args.coverage_sigma, size=n)
    if args.bg_kappa > 0:
        pbg = rng.dirichlet(args.bg_kappa * q, size=n)
    else:
        pbg = np.tile(q, (n, 1))
    PA = pbg.copy()
    PB = pbg.copy()
    mA = np.ones(n)
    mB = np.ones(n)
    cfg_idx = np.full(n, -1, dtype=int)
    bclass = np.array(["null"] * n, dtype=object)
    for blk in plan:
        c = cfgs[blk["cfg"]]
        bs, be = blk["bin_start"], blk["bin_end"]
        if c[2] is not None:
            PA[bs:be] = c[2]
            mA[bs:be] = args.feature_mult * c[4]
        if c[3] is not None:
            PB[bs:be] = c[3]
            mB[bs:be] = args.feature_mult * c[5]
        cfg_idx[bs:be] = blk["cfg"]
        bclass[bs:be] = "interior"
        e = max(1, int(math.ceil(2 * args.smooth_sigma_bins)))
        bclass[bs:min(be, bs + e)] = "edge"
        bclass[max(bs, be - e):be] = "edge"
    # bins next to a block (smoothing bleed) are not clean nulls
    e = max(1, int(math.ceil(3 * args.smooth_sigma_bins)))
    near = np.zeros(n, dtype=bool)
    for blk in plan:
        near[max(0, blk["bin_start"] - e):blk["bin_start"]] = True
        near[blk["bin_end"]:min(n, blk["bin_end"] + e)] = True
    bclass[near & (bclass == "null")] = "flank"

    expA = (args.depth * cov * mA)[:, None] * PA       # expected, BigWig units
    expB = (args.depth * cov * mB)[:, None] * PB
    samples = {}
    for cond, ex in (("A", expA), ("B", expB)):
        for r in range(args.reps):
            lib = rng.lognormal(0.0, args.lib_sigma)
            mu = np.clip(ex * lib * args.counts_per_unit, 1e-9, None)
            rr = float(args.nb_dispersion)
            x = rng.negative_binomial(rr, rr / (rr + mu)).astype(np.float64)
            x /= args.counts_per_unit
            if args.smooth_sigma_bins > 0:
                for t in range(N_TRACKS):
                    x[:, t] = gaussian_filter1d(x[:, t], args.smooth_sigma_bins,
                                                mode="wrap")
            samples["%s_rep%d" % (cond, r + 1)] = np.clip(x, 0.0, None)
    truth = dict(cfgs=cfgs, plan=plan, cfg_idx=cfg_idx, bclass=bclass,
                 expA=expA, expB=expB, q=q, n=n)
    return samples, truth


# ---------------------------------------------------------------------------
# BEARING (production functions)
# ---------------------------------------------------------------------------

def bearing_scores(x, min_signal, alpha, q=None):
    """Score one sample exactly as bigwig_to_qcat does for an un-normalised
    run: P = signals_to_prob(raw), Q = mean of P over all bins (unless q is
    given), clamped KL with the low-signal mask and the composition prior."""
    P = b2q.signals_to_prob(x)
    Q = P.mean(axis=0) if q is None else q
    s, _ = b2q.kl_scores_per_bin(P, Q, raw_signal_matrix=x, min_signal=min_signal,
                                 prior_strength=alpha, prior_source=x)
    return s, Q


def bearing_diff(sample_mats, reps, min_signal, alpha, qs=None):
    a = [bearing_scores(sample_mats["A_rep%d" % (r + 1)], min_signal, alpha,
                        None if qs is None else qs[0])[0] for r in range(reps)]
    b = [bearing_scores(sample_mats["B_rep%d" % (r + 1)], min_signal, alpha,
                        None if qs is None else qs[1])[0] for r in range(reps)]
    return np.mean(a, axis=0) - np.mean(b, axis=0)


def shifted(x, rng, min_shift_bins):
    n = x.shape[0]
    out = np.empty_like(x)
    hi = max(min_shift_bins + 1, n - min_shift_bins)
    for t in range(x.shape[1]):
        out[:, t] = np.roll(x[:, t], int(rng.integers(min_shift_bins, hi)))
    return out


BEARING_STATS = {
    # production: bearing_pvalue.py --diff tests |sum over tracks of the
    # signed per-track differential| = |total KL_A - total KL_B|
    "sum": lambda d: np.abs(d.sum(axis=1)),
    # bearing_pvalue.py --diff-stat l1 (config diff_stat: l1): L1 norm of the
    # per-track differential, sum over tracks of |KL_A,i - KL_B,i|; does not
    # cancel when one track's gain is another track's loss
    "l1": lambda d: np.abs(d).sum(axis=1),
}


def run_bearing(samples, args, rng, stat_names=("sum", "l1")):
    """Return {stat_name: result dict}. One set of permutations serves every
    statistic (same shifted samples), so the arms differ only in the
    statistic."""
    d = bearing_diff(samples, args.reps, args.min_signal, args.prior_strength)
    obs = {k: BEARING_STATS[k](d) for k in stat_names}
    msb = args.min_shift // args.bin_size
    null = {k: [] for k in stat_names}
    for _ in range(args.n_perms):
        sh = {nm: shifted(x, rng, msb) for nm, x in samples.items()}
        dn = bearing_diff(sh, args.reps, args.min_signal, args.prior_strength)
        for k in stat_names:
            v = BEARING_STATS[k](dn)
            null[k].append(v[v >= args.pvalue_min_signal])
    out = {}
    for k in stat_names:
        st = obs[k]
        nl = np.sort(np.concatenate(null[k])) if null[k] else np.array([])
        tested = st >= args.pvalue_min_signal
        p = np.ones(len(st))
        q = np.ones(len(st))
        if tested.any():
            if nl.size == 0:
                sys.exit("ERROR: empty BEARING null for statistic %s (no permuted "
                         "bin reached --pvalue-min-signal); lower it or raise "
                         "--n-perms" % k)
            p[tested] = bpv.empirical_pvals(st[tested], nl)
            q[tested] = bpv.bh_fdr(p[tested], args.fdr)[1]
        out[k] = dict(diff=d, total=d.sum(axis=1), stat=st, tested=tested, p=p,
                      q=q, called=tested & (q <= args.fdr), n_null=int(nl.size))
    return out


def expected_bearing(truth, args):
    """BEARING differential on the noise-free expectation (library 1), with the
    background Q of the noisy samples replaced by the expectation's own Q.
    Tells whether BEARING SHOULD respond to a block in principle."""
    exp = {}
    for r in range(args.reps):
        exp["A_rep%d" % (r + 1)] = truth["expA"]
        exp["B_rep%d" % (r + 1)] = truth["expB"]
    return bearing_diff(exp, args.reps, args.min_signal, args.prior_strength)


# ---------------------------------------------------------------------------
# Per-assay moderated t (limma-trend style)
# ---------------------------------------------------------------------------

def _trigamma_inverse(x):
    # Newton iteration of Smyth (2004) / limma::trigammaInverse
    x = np.asarray(x, dtype=np.float64)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(np.abs(-dif / y) < 1e-8):
            break
    return y


def fit_f_dist(s2, df):
    """Smyth (2004) moment estimate of the prior (d0, s0^2) for variances s2
    with residual df. Returns (d0, s02); d0 = inf when there is no excess
    variance."""
    s2 = np.asarray(s2, dtype=np.float64)
    s2 = s2[np.isfinite(s2) & (s2 > 0)]
    if s2.size < 10:
        return np.inf, float(np.median(s2)) if s2.size else 1.0
    z = np.log(s2)
    e = z - special.digamma(df / 2.0) + math.log(df / 2.0)
    emean = float(e.mean())
    evar = float(e.var(ddof=1)) - float(special.polygamma(1, df / 2.0))
    if evar > 0:
        d0 = float(2.0 * _trigamma_inverse(evar))
        s02 = math.exp(emean + special.digamma(d0 / 2.0) - math.log(d0 / 2.0))
    else:
        d0 = np.inf
        s02 = math.exp(emean)
    return d0, s02


def run_per_assay(samples, args, n_strata=20):
    A = [samples["A_rep%d" % (r + 1)] for r in range(args.reps)]
    B = [samples["B_rep%d" % (r + 1)] for r in range(args.reps)]
    allm = A + B
    n = allm[0].shape[0]
    nA, nB = len(A), len(B)
    df = nA + nB - 2
    raw_mean = np.mean(allm, axis=0)                      # (n, 6)
    tested = raw_mean >= args.assay_min_signal
    # library normalisation per sample per track: median-of-ratios size
    # factors (DESeq), robust to the planted blocks changing a track's total
    sf = np.ones((len(allm), N_TRACKS))
    stack = np.stack(allm)                                # (samples, n, 6)
    for t in range(N_TRACKS):
        v = stack[:, :, t]
        ok = (v > 0).all(axis=0)
        if ok.sum() >= 100:
            lg = np.log(v[:, ok])
            sf[:, t] = np.exp(np.median(lg - lg.mean(axis=0), axis=1))
    norm = [allm[i] / sf[i][None, :] for i in range(len(allm))]
    y = [np.log2(m + args.assay_pseudo) for m in norm]
    yA = np.stack(y[:nA])
    yB = np.stack(y[nA:])
    mA_, mB_ = yA.mean(axis=0), yB.mean(axis=0)
    ss = ((yA - mA_) ** 2).sum(axis=0) + ((yB - mB_) ** 2).sum(axis=0)
    s2 = ss / df
    ave = np.concatenate([yA, yB]).mean(axis=0)
    tstat = np.zeros((n, N_TRACKS))
    pval = np.ones((n, N_TRACKS))
    for t in range(N_TRACKS):
        idx = np.where(tested[:, t])[0]
        if idx.size == 0:
            continue
        qs = np.quantile(ave[idx, t], np.linspace(0, 1, n_strata + 1))
        strata = np.clip(np.searchsorted(qs, ave[idx, t], side="right") - 1,
                         0, n_strata - 1)
        for k in range(n_strata):
            sub = idx[strata == k]
            if sub.size == 0:
                continue
            d0, s02 = fit_f_dist(s2[sub, t], df)
            if np.isinf(d0):
                post = np.full(sub.size, s02)
                dft = 1e6
            else:
                post = (d0 * s02 + df * s2[sub, t]) / (d0 + df)
                dft = d0 + df
            se = np.sqrt(np.clip(post, 1e-300, None) * (1.0 / nA + 1.0 / nB))
            tt = (mA_[sub, t] - mB_[sub, t]) / se
            tstat[sub, t] = tt
            pval[sub, t] = 2.0 * stats.t.sf(np.abs(tt), dft)
    # bin-level call: Simes combination of the tested tracks' p-values, BH
    # over bins (controls the bin-level FDR; BH over bin x track tests and
    # "any track" would not, because a true bin usually has several
    # significant tracks and a false one only one)
    pm = np.where(tested, pval, np.nan)
    ntest = tested.sum(axis=1)
    ps = np.sort(np.where(tested, pval, np.inf), axis=1)
    ranks = np.arange(1, N_TRACKS + 1)[None, :]
    simes = np.min(np.where(np.isfinite(ps), ps * ntest[:, None] / ranks, np.inf),
                   axis=1)
    simes = np.clip(simes, 0.0, 1.0)
    anyt = ntest > 0
    qbin = np.ones(n)
    if anyt.any():
        qbin[anyt] = bpv.bh_fdr(simes[anyt], args.fdr)[1]
    called = anyt & (qbin <= args.fdr)
    # which tracks moved, within called bins: track-level BH over all tested
    # bin x track pairs; if none passes, the bin's smallest-p track
    q = np.ones((n, N_TRACKS))
    flat = tested.ravel()
    if flat.any():
        q.ravel()[flat] = bpv.bh_fdr(pval.ravel()[flat], args.fdr)[1]
    sig = tested & (q <= args.fdr) & called[:, None]
    best = np.nanargmin(np.where(tested, pm, np.inf), axis=1)
    none = called & ~sig.any(axis=1)
    sig[np.where(none)[0], best[none]] = True
    pos = (sig & (tstat > 0)).any(axis=1)
    neg = (sig & (tstat < 0)).any(axis=1)
    pattern = np.where(~called, "", np.where(pos & neg, "discordant", "concordant"))
    return dict(t=tstat, p=pval, q=q, tested=tested, sig=sig, called=called,
                q_bin=qbin, pattern=pattern, n_sig_tracks=sig.sum(axis=1))


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def _frac(x):
    x = np.asarray(x)
    return float(x.mean()) if x.size else float("nan")


TYPES = ["feature_null", "magnitude", "switch", "rebalance", "both", "onoff"]
TRUE_CHANGE = {"feature_null": False, "magnitude": True, "switch": True,
               "rebalance": True, "both": True, "onoff": True}


def evaluate(truth, bea, pa, expd, args):
    """bea: {stat: result}; 'sum' and 'l1' are the two --diff-stat options."""
    cfgs = truth["cfgs"]
    ci = truth["cfg_idx"]
    interior = truth["bclass"] == "interior"
    exp_stat = {k: BEARING_STATS[k](expd) for k in bea}
    exp_top = np.argmax(np.abs(expd), axis=1)
    obs_top = np.argmax(np.abs(bea["sum"]["diff"]), axis=1)
    rows = []
    for k, c in enumerate(cfgs):
        m = interior & (ci == k)
        r = dict(config=c[0], type=c[1], s=c[6], fold=c[7], delta=c[8],
                 n_bins=int(m.sum()))
        for st, res in bea.items():
            cb = res["called"][m]
            r["E_stat_%s" % st] = float(np.median(exp_stat[st][m]))
            r["bearing_%s_tested" % st] = _frac(res["tested"][m])
            r["bearing_%s_called" % st] = _frac(cb)
        cb = bea["sum"]["called"][m]
        r["expected_top_track"] = TRACK_NAMES[int(np.bincount(
            exp_top[m], minlength=N_TRACKS).argmax())]
        r["bearing_sum_top_track_correct"] = (
            _frac(obs_top[m][cb] == exp_top[m][cb]) if cb.any() else float("nan"))
        r["per_assay_called"] = _frac(pa["called"][m])
        r["per_assay_concordant"] = _frac(pa["pattern"][m] == "concordant")
        r["per_assay_discordant"] = _frac(pa["pattern"][m] == "discordant")
        r["per_assay_median_sig_tracks"] = float(np.median(pa["n_sig_tracks"][m]))
        r["bearing_sum_only"] = _frac(cb & ~pa["called"][m])
        r["per_assay_only"] = _frac(~cb & pa["called"][m])
        rows.append(r)

    nullm = truth["bclass"] == "null"
    groups = [(ty, interior & np.isin(ci, [k for k, c in enumerate(cfgs)
                                            if c[1] == ty])) for ty in TYPES]
    groups.append(("null_bins", nullm))
    trows = []
    for ty, m in groups:
        r = dict(type=ty, n_bins=int(m.sum()))
        for st, res in bea.items():
            r["bearing_%s_called" % st] = _frac(res["called"][m])
        r["per_assay_called"] = _frac(pa["called"][m])
        r["bearing_sum_only"] = _frac(bea["sum"]["called"][m] & ~pa["called"][m])
        r["per_assay_only"] = _frac(~bea["sum"]["called"][m] & pa["called"][m])
        r["per_assay_concordant"] = _frac(pa["pattern"][m] == "concordant")
        r["per_assay_discordant"] = _frac(pa["pattern"][m] == "discordant")
        trows.append(r)

    # empirical FDR: calls in no-change bins (clean null + feature_null
    # interiors) over all calls in scored bins (those + true-change interiors)
    fn_k = [k for k, c in enumerate(cfgs) if not TRUE_CHANGE[c[1]]]
    nochange = nullm | (interior & np.isin(ci, fn_k))
    scored = nochange | (interior & (ci >= 0))
    fdr = {}
    for st, res in bea.items():
        tot = int((res["called"] & scored).sum())
        fdr["bearing_%s" % st] = (int((res["called"] & nochange).sum()), tot)
    tot = int((pa["called"] & scored).sum())
    fdr["per_assay"] = (int((pa["called"] & nochange).sum()), tot)
    return rows, trows, fdr


def _fmt(v):
    if isinstance(v, float):
        return "NA" if not np.isfinite(v) else "%.4g" % v
    return str(v)


def write_tsv(path, rows):
    cols = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r[c]) for c in cols) + "\n")


def write_per_bin(path, truth, bea, pa, expd, args):
    cfgs = truth["cfgs"]
    s, l1 = bea["sum"], bea["l1"]
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(["chrom", "start", "end", "config", "type", "bin_class",
                            "exp_sum_diff", "sum_diff", "sum_tested", "sum_p",
                            "sum_q", "l1_stat", "l1_q", "top_track",
                            "per_assay_q_bin", "per_assay_n_sig",
                            "per_assay_pattern"]) + "\n")
        for i in range(truth["n"]):
            k = truth["cfg_idx"][i]
            cname, ctype = (cfgs[k][0], cfgs[k][1]) if k >= 0 else ("null", "null")
            fh.write("\t".join([
                args.chrom, str(i * args.bin_size), str((i + 1) * args.bin_size),
                cname, ctype, truth["bclass"][i],
                "%.4g" % expd[i].sum(), "%.4g" % s["total"][i],
                "1" if s["tested"][i] else "0", "%.3g" % s["p"][i],
                "%.3g" % s["q"][i], "%.4g" % l1["stat"][i], "%.3g" % l1["q"][i],
                TRACK_NAMES[int(np.argmax(np.abs(s["diff"][i])))],
                "%.3g" % pa["q_bin"][i], str(int(pa["n_sig_tracks"][i])),
                pa["pattern"][i] or "."]) + "\n")


def summary_text(args, bea, trows, rows, fdr):
    out = []
    out.append("complementarity_benchmark.py 1.0.0  (SIMULATED data)")
    out.append("depth %g  feature_mult %g  reps %d  prior_strength %g  min_signal %g  "
               "pvalue_min_signal %g  n_perms %d  seed %d"
               % (args.depth, args.feature_mult, args.reps, args.prior_strength,
                  args.min_signal, args.pvalue_min_signal, args.n_perms, args.seed))
    for st, lab in (("sum", "BEARING production |sum d|"),
                    ("l1", "BEARING --diff-stat l1 sum|d|")):
        res = bea[st]
        out.append("%-30s %7d bins tested, %6d called (q <= %g), pooled null %d"
                   % (lab, int(res["tested"].sum()), int(res["called"].sum()),
                      args.fdr, res["n_null"]))
    out.append("Empirical FDR (calls in no-change bins / all calls): " + "; ".join(
        "%s %d/%d = %s" % (k, v[0], v[1], ("%.3f" % (v[0] / v[1])) if v[1] else "NA")
        for k, v in fdr.items()))
    out.append("")
    out.append("Call rate by block type (fraction of block-interior bins; null_bins = "
               "false-call rate)")
    out.append("%-13s %7s %8s %8s %9s %8s %8s %9s %9s"
               % ("type", "bins", "B_sum", "B_L1", "perAssay", "Bsum_only",
                  "PA_only", "PA_concord", "PA_discord"))
    for r in trows:
        out.append("%-13s %7d %8.3f %8.3f %9.3f %8.3f %8.3f %9.3f %9.3f"
                   % (r["type"], r["n_bins"], r["bearing_sum_called"],
                      r["bearing_l1_called"], r["per_assay_called"],
                      r["bearing_sum_only"], r["per_assay_only"],
                      r["per_assay_concordant"], r["per_assay_discordant"]))
    out.append("")
    out.append("%-26s %7s %7s %8s %7s %7s %9s %7s"
               % ("config", "E|sum|", "E_L1", "Bsum_tst", "B_sum", "B_L1",
                  "perAssay", "topOK"))
    for r in rows:
        out.append("%-26s %7.3f %7.3f %8.3f %7.3f %7.3f %9.3f %7s"
                   % (r["config"], r["E_stat_sum"], r["E_stat_l1"],
                      r["bearing_sum_tested"], r["bearing_sum_called"],
                      r["bearing_l1_called"], r["per_assay_called"],
                      _fmt(r["bearing_sum_top_track_correct"])))
    out.append("")
    out.append("E|sum| / E_L1: the statistic on the noise-free expectation (testing "
               "floor %g). B_sum = BEARING differential with --diff-stat sum (|sum over "
               "tracks|, the original); B_L1 = --diff-stat l1 (sum of |per-track "
               "differential|). topOK = BEARING's largest |per-track differential| is the "
               "expected track (B_sum calls)." % args.pvalue_min_signal)
    return "\n".join(out) + "\n"


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chrom", default="chrSim")
    ap.add_argument("--length", type=int, default=20_000_000)
    ap.add_argument("--bin-size", type=int, default=200)
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--copies", type=int, default=6)
    ap.add_argument("--block-bins", type=int, default=25)
    ap.add_argument("--gap-bins", type=int, default=175)
    ap.add_argument("--start-offset-bins", type=int, default=500)
    ap.add_argument("--depth", type=float, default=2.0,
                    help="expected total signal of a null bin, BigWig units (default 2)")
    ap.add_argument("--feature-mult", type=float, default=4.0,
                    help="depth multiplier of planted feature blocks (default 4)")
    ap.add_argument("--bg-kappa", type=float, default=50.0,
                    help="Dirichlet concentration of per-bin background "
                         "composition around Q; 0 = every null bin is Q (default 50)")
    ap.add_argument("--coverage-sigma", type=float, default=0.5)
    ap.add_argument("--lib-sigma", type=float, default=0.1)
    ap.add_argument("--nb-dispersion", type=float, default=8.0)
    ap.add_argument("--counts-per-unit", type=float, default=10.0,
                    help="read counts per BigWig signal unit (sets Poisson "
                         "granularity; default 10)")
    ap.add_argument("--smooth-sigma-bins", type=float, default=1.0)
    # BEARING, production defaults (results_prior1_n100)
    ap.add_argument("--prior-strength", type=float, default=1.0)
    ap.add_argument("--min-signal", type=float, default=0.1)
    ap.add_argument("--pvalue-min-signal", type=float, default=0.5)
    ap.add_argument("--n-perms", type=int, default=20)
    ap.add_argument("--min-shift", type=int, default=1_500_000)
    ap.add_argument("--fdr", type=float, default=0.05)
    # per-assay comparator
    ap.add_argument("--assay-min-signal", type=float, default=0.1)
    ap.add_argument("--assay-pseudo", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--per-bin", action="store_true")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.reps < 2:
        sys.exit("ERROR: --reps must be >= 2 (the per-assay test needs replicates)")
    if 2 * args.min_shift >= args.length:
        sys.exit("ERROR: --min-shift must be < half of --length")
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    samples, truth = simulate(args, rng)
    bea = run_bearing(samples, args, np.random.default_rng(args.seed + 1))
    pa = run_per_assay(samples, args)
    expd = expected_bearing(truth, args)
    rows, trows, fdr = evaluate(truth, bea, pa, expd, args)
    write_tsv(os.path.join(args.outdir, "per_config.tsv"), rows)
    write_tsv(os.path.join(args.outdir, "per_type.tsv"), trows)
    txt = summary_text(args, bea, trows, rows, fdr)
    with open(os.path.join(args.outdir, "summary.txt"), "w") as fh:
        fh.write(txt)
    if args.per_bin:
        write_per_bin(os.path.join(args.outdir, "per_bin.tsv.gz"), truth, bea, pa,
                      expd, args)
    sys.stdout.write(txt)
    return rows, trows, fdr


if __name__ == "__main__":
    main()
