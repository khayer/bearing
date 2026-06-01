#!/usr/bin/env python3
"""
bearing_calibration.py
======================
Within-condition replicate-differential calibration check for the BEARING
differential FDR pipeline.

PURPOSE
-------
Two biological replicates of the SAME condition should differ only by noise;
there is no true differential signal between them. If you run the exact same
differential significance machinery used for real comparisons
(compare_qcat.py diff -> bearing_pvalue.py --diff) on a within-condition
rep1-vs-rep2 differential, the resulting per-bin p-values should be
approximately UNIFORM on [0, 1]. Departure from uniformity (inflation toward
small p, lambda >> 1, an excess of BH-significant bins) means the differential
FDR is mis-calibrated and would over-call false positives in real comparisons
too. This is the cleanest empirical calibration check for the FDR claims that
the regional differential statistics depend on.

This script does NOT reimplement scoring or p-value computation. It imports
bearing_pvalue.py and reuses, unchanged:
    parse_qcat, collect_null_scores, empirical_pvals, bh_fdr
so the p-values it evaluates are the same statistic bearing_pvalue.py --diff
reports. A hard preflight verifies the call signatures and fails loudly if
they have drifted (it will not guess a different interface).

"two-sided" == the --diff folded-magnitude test
-----------------------------------------------
bearing_pvalue.py --diff tests P(|X| >= |s|): both the observed signed
differential scores and the permuted null are folded to absolute value and
compared with the one-sided empirical survival function (the sign is recorded
separately). That folded magnitude test IS the differential statistic. This
calibration therefore evaluates empirical_pvals(|obs|, sorted|null|), not an
independent two-tailed estimator.

INPUTS
------
--diff-qcat   observed within-condition rep1-vs-rep2 diff qcat.bgz.
              Produce it with compare_qcat.py by listing rep1 and rep2 as two
              single-replicate "conditions", so the emitted diff is
              rep-vs-rep, not condition-mean-vs-condition-mean.
--null-qcat   one or more permuted rep-vs-rep diff qcat.bgz files. Produce
              them with generate_perm_nulls.py --diff-sheet on the same
              rep-as-condition sheet.

The --min-signal mask is applied IDENTICALLY to observed and null (matching
the BEARING principle that blacklist and low-signal masks are applied the same
way to real and permuted runs).

OUTPUTS (--out-prefix PREFIX)
-----------------------------
PREFIX_calibration.pdf   four-panel figure: p-value histogram, midrank uniform
                         QQ, observed-vs-null magnitude distribution, p-value
                         ECDF vs uniform.
PREFIX_perbin.tsv.gz     per-bin chrom/start/end/bearing_score/
                         bearing_score_tested/direction/pval/pval_adj_bh/
                         significant_fdr<alpha>.
PREFIX_summary.txt       one-line verdict plus the calibration numbers.

USAGE
-----
    python bearing_calibration.py \
        --diff-qcat   DN_rep1_vs_rep2.qcat.bgz \
        --null-qcat   perm1/diff_DN_rep1_vs_rep2.qcat.bgz \
                      perm2/diff_DN_rep1_vs_rep2.qcat.bgz \
                      perm3/diff_DN_rep1_vs_rep2.qcat.bgz \
        --out-prefix  calib/DN \
        --min-signal  0.0 \
        --fdr 0.05 \
        --label "DN rep1 vs rep2"

DEPENDENCIES
------------
    Python 3.8+
    numpy, scipy, matplotlib
    bearing_pvalue.py (and its deps, including pyBigWig) importable
"""

import argparse
import gzip
import importlib.util
import inspect
import sys
from pathlib import Path

import numpy as np
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Reused, unchanged, from bearing_pvalue.py. The preflight checks that each of
# these exists and exposes the parameter names this script depends on.
REQUIRED_SIGNATURES = {
    "parse_qcat":          ("path", "min_signal", "diff_mode"),
    "collect_null_scores": ("null_paths", "min_signal", "diff_mode"),
    "empirical_pvals":     ("observed_scores", "null_scores_sorted"),
    "bh_fdr":              ("pvals", "alpha"),
}

CHI2_MEDIAN_DF1 = float(stats.chi2.ppf(0.5, df=1))  # 0.4549363...


# ---------------------------------------------------------------------------
# Import + preflight (no scoring is ever reimplemented here)
# ---------------------------------------------------------------------------

def import_bearing_pvalue(path_hint):
    """
    Import bearing_pvalue.py as a module. Resolution order:
      1. explicit --bearing-pvalue path (a file, or a directory containing it)
      2. the directory of this script
      3. the current working directory
    """
    candidates = []
    if path_hint:
        p = Path(path_hint)
        candidates.append(p / "bearing_pvalue.py" if p.is_dir() else p)
    here = Path(__file__).resolve().parent
    candidates.append(here / "bearing_pvalue.py")
    candidates.append(Path.cwd() / "bearing_pvalue.py")

    for cand in candidates:
        if cand.is_file():
            spec = importlib.util.spec_from_file_location("bearing_pvalue",
                                                          str(cand))
            mod = importlib.util.module_from_spec(spec)
            sys.modules["bearing_pvalue"] = mod
            spec.loader.exec_module(mod)
            print("Imported bearing_pvalue from: %s" % cand, file=sys.stderr)
            return mod

    sys.exit(
        "ERROR: could not locate bearing_pvalue.py. Tried:\n  "
        + "\n  ".join(str(c) for c in candidates)
        + "\nPass its path with --bearing-pvalue."
    )


def preflight(bp):
    """
    Verify the reused functions exist with the parameter names this script
    depends on. Fail loudly (exit 2) on any mismatch.
    """
    problems = []
    for name, required_params in REQUIRED_SIGNATURES.items():
        fn = getattr(bp, name, None)
        if fn is None or not callable(fn):
            problems.append("  missing function: %s" % name)
            continue
        try:
            sig = inspect.signature(fn)
        except (TypeError, ValueError):
            problems.append("  cannot introspect signature: %s" % name)
            continue
        have = set(sig.parameters.keys())
        for rp in required_params:
            if rp not in have:
                problems.append(
                    "  %s: required parameter '%s' not found "
                    "(actual signature: %s)" % (name, rp, sig)
                )
    if problems:
        sys.exit(
            "PREFLIGHT FAILED -- bearing_pvalue.py interface has drifted:\n"
            + "\n".join(problems)
            + "\n\nThis script reuses these functions unchanged and will not "
              "guess a different interface. Update bearing_calibration.py to "
              "match, or pass the correct --bearing-pvalue."
        )
    print("Preflight OK: parse_qcat / collect_null_scores / empirical_pvals / "
          "bh_fdr present with expected parameters.", file=sys.stderr)


# ---------------------------------------------------------------------------
# Load observed (reuses bp.parse_qcat in diff mode)
# ---------------------------------------------------------------------------

def load_observed(bp, diff_qcat, min_signal):
    """
    Parse the observed within-condition rep-vs-rep diff qcat in diff mode
    (filters by abs(score) >= min_signal, keeps the signed total).

    Returns (rows, signed, tested):
      rows   : list of (chrom, start, end)
      signed : np.float64 array of signed differential scores
      tested : np.float64 array of |signed| (the value the null tests)
    """
    rows = []
    signed = []
    for chrom, start, end, score_total, _per_track in bp.parse_qcat(
            diff_qcat, min_signal=min_signal, diff_mode=True):
        rows.append((chrom, start, end))
        signed.append(score_total)
    if not rows:
        sys.exit("ERROR: no observed bins passed min_signal=%.4g in %s"
                 % (min_signal, diff_qcat))
    signed = np.asarray(signed, dtype=np.float64)
    return rows, signed, np.abs(signed)


# ---------------------------------------------------------------------------
# Calibration statistics
# ---------------------------------------------------------------------------

def genomic_inflation(pvals, p_floor):
    """
    Genomic-inflation lambda from one-df chi-square quantiles:
        lambda = median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1)
    Under Uniform(0,1) p, the median of chi2.isf(p,1) is chi2.ppf(0.5,1),
    so lambda -> 1.

    Guard: if more than 50% of bins sit at the permutation floor (the minimum
    achievable p), the median is pinned at the floor and lambda is just a
    function of the null size, not a calibration measure. Return NaN and
    flag saturation in that case.

    Returns (lambda_value, frac_at_floor, saturated_flag).
    """
    n = len(pvals)
    at_floor = int(np.count_nonzero(pvals <= p_floor * (1.0 + 1e-9)))
    frac_floor = at_floor / n if n else float("nan")
    if frac_floor > 0.5:
        return float("nan"), frac_floor, True
    p = np.clip(pvals, np.finfo(np.float64).tiny, 1.0)
    chi2_stats = stats.chi2.isf(p, df=1)
    lam = float(np.median(chi2_stats) / CHI2_MEDIAN_DF1)
    return lam, frac_floor, False


def midrank_qq_points(pvals, max_points):
    """
    Midrank uniform-QQ positions on the -log10 scale. Average ranks give tied
    p-values (from the discrete empirical null) a shared plotting position.

    Returns (exp_neglog10, obs_neglog10), subsampled to <= max_points with the
    most-significant tail kept in full.
    """
    n = len(pvals)
    ranks = stats.rankdata(pvals, method="average")
    expected = ranks / (n + 1.0)
    tiny = np.finfo(np.float64).tiny
    x = -np.log10(np.clip(expected, tiny, 1.0))
    y = -np.log10(np.clip(pvals, tiny, 1.0))
    if n > max_points:
        order = np.argsort(pvals)
        n_tail = min(2000, max(1, max_points // 4))
        tail = order[:n_tail]
        rest = order[n_tail:]
        step = max(1, len(rest) // max(1, (max_points - n_tail)))
        keep = np.concatenate([tail, rest[::step]])
        x, y = x[keep], y[keep]
    return x, y


def subsample(arr, max_points, seed=42):
    """Uniform random subsample of a 1-D array (for plotting only)."""
    if len(arr) <= max_points:
        return arr
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(arr), size=max_points, replace=False)
    return arr[idx]


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(out_pdf, label, pvals, tested, null_sorted, p_floor,
                lam, frac_floor, saturated, n_sig, fdr, max_qq_points):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    title = "BEARING differential FDR calibration"
    if label:
        title += "  --  %s" % label
    fig.suptitle(title, fontsize=13, fontweight="bold")

    # Panel A: p-value histogram (uniform under a calibrated null)
    ax = axes[0, 0]
    ax.hist(pvals, bins=50, range=(0.0, 1.0), density=True,
            color="#0075C9", edgecolor="white", linewidth=0.3)
    ax.axhline(1.0, color="#FFC72C", linestyle="--", linewidth=1.5,
               label="uniform (calibrated)")
    ax.set_xlabel("empirical p-value")
    ax.set_ylabel("density")
    ax.set_title("A. P-value histogram")
    ax.legend(fontsize=8, loc="upper right")

    # Panel B: midrank uniform QQ on -log10 scale
    ax = axes[0, 1]
    xq, yq = midrank_qq_points(pvals, max_qq_points)
    lim = max(xq.max() if len(xq) else 1.0, yq.max() if len(yq) else 1.0)
    ax.plot([0, lim], [0, lim], color="#FFC72C", linestyle="--",
            linewidth=1.5, label="y = x")
    ax.scatter(xq, yq, s=4, alpha=0.4, color="#003087", edgecolors="none")
    if saturated:
        lam_txt = "lambda = undefined (floor pileup)"
    else:
        lam_txt = "lambda = %.3f" % lam
    ax.set_xlabel("expected -log10(p) [uniform]")
    ax.set_ylabel("observed -log10(p)")
    ax.set_title("B. Midrank QQ  (%s)" % lam_txt)
    ax.legend(fontsize=8, loc="upper left")

    # Panel C: observed |diff| vs pooled null |diff|
    ax = axes[1, 0]
    null_plot = subsample(null_sorted, 2_000_000)
    hi = float(max(np.quantile(tested, 0.9995) if len(tested) else 1.0,
                   np.quantile(null_sorted, 0.9995) if len(null_sorted) else 1.0))
    hi = hi if hi > 0 else 1.0
    bins = np.linspace(0.0, hi, 80)
    ax.hist(null_plot, bins=bins, density=True, histtype="step",
            color="#8b0000", linewidth=1.5, label="permutation null |diff|")
    ax.hist(tested, bins=bins, density=True, histtype="step",
            color="#00A9CE", linewidth=1.5, label="observed |diff|")
    ax.set_yscale("log")
    ax.set_xlabel("|differential bearing score|")
    ax.set_ylabel("density (log)")
    ax.set_title("C. Observed vs null magnitude")
    ax.legend(fontsize=8, loc="upper right")

    # Panel D: p-value ECDF vs uniform
    ax = axes[1, 1]
    p_sorted = np.sort(subsample(pvals, 200_000))
    ecdf_y = np.arange(1, len(p_sorted) + 1) / len(p_sorted)
    ax.plot([0, 1], [0, 1], color="#FFC72C", linestyle="--", linewidth=1.5,
            label="uniform")
    ax.plot(p_sorted, ecdf_y, color="#003087", linewidth=1.5, label="observed")
    ax.set_xlabel("p-value")
    ax.set_ylabel("empirical CDF")
    ax.set_title("D. P-value ECDF vs uniform")
    ax.legend(fontsize=8, loc="lower right")

    foot = ("n bins = %s   |   at floor (p = %.2e) = %s (%.2f%%)   |   "
            "BH-significant at FDR %.3g = %s") % (
                "{:,}".format(len(pvals)), p_floor,
                "{:,}".format(int(round(frac_floor * len(pvals)))),
                100.0 * frac_floor, fdr, "{:,}".format(int(n_sig)))
    fig.text(0.5, 0.005, foot, ha="center", fontsize=9)

    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(out_pdf, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-bin TSV
# ---------------------------------------------------------------------------

def write_perbin_tsv(out_path, rows, signed, tested, pvals, padj, rejected, fdr):
    sig_col = "significant_fdr%g" % fdr
    with gzip.open(out_path, "wt") as fh:
        fh.write("\t".join([
            "chrom", "start", "end", "bearing_score", "bearing_score_tested",
            "direction", "pval", "pval_adj_bh", sig_col]) + "\n")
        for i, (chrom, start, end) in enumerate(rows):
            s = signed[i]
            direction = "+" if s > 0 else ("-" if s < 0 else "0")
            fh.write("\t".join([
                chrom, str(start), str(end),
                "%.6g" % s, "%.6g" % tested[i], direction,
                "%.6e" % pvals[i], "%.6e" % padj[i],
                "1" if rejected[i] else "0"]) + "\n")


# ---------------------------------------------------------------------------
# Verdict + summary
# ---------------------------------------------------------------------------

def build_verdict(lam, frac_floor, saturated, n_sig, n_bins, fdr):
    """
    Heuristic verdict. The brief's principle holds: judge calibration on the
    bulk of the distribution, not the saturated tail. Thresholds below are a
    starting read, not a substitute for inspecting panels A, B and D.
    """
    sig_frac = n_sig / n_bins if n_bins else float("nan")
    if saturated:
        return ("SATURATED -- more than 50%% of bins sit at the permutation "
                "floor, so lambda is undefined. The null is too sparse for "
                "this observed differential: add permutation rounds (raise the "
                "null bin count) and/or raise --min-signal. Calibration cannot "
                "be read until this is resolved. BH-significant: %s / %s." % (
                    "{:,}".format(int(n_sig)), "{:,}".format(int(n_bins))))
    parts = []
    if 0.90 <= lam <= 1.15:
        parts.append("CALIBRATED -- lambda = %.3f is close to 1" % lam)
    elif lam > 1.15:
        parts.append("INFLATED -- lambda = %.3f exceeds 1; the differential "
                     "FDR likely over-calls in real comparisons" % lam)
    else:
        parts.append("CONSERVATIVE -- lambda = %.3f is below 1; p-values skew "
                     "large (under-powered, not anti-conservative)" % lam)
    parts.append("BH-significant at FDR %.3g: %s / %s (%.4f%%); a calibrated "
                 "within-replicate null expects near zero" % (
                     fdr, "{:,}".format(int(n_sig)), "{:,}".format(int(n_bins)),
                     100.0 * sig_frac))
    parts.append("bins at permutation floor: %.2f%%" % (100.0 * frac_floor))
    return "  ".join(parts)


def write_summary(out_path, label, n_bins, n_null, p_floor, lam, frac_floor,
                  saturated, n_sig, fdr, ks_stat, ks_caveat, verdict):
    lam_txt = "undefined (floor pileup)" if saturated else "%.4f" % lam
    lines = [
        "BEARING differential FDR calibration summary",
        "label                 : %s" % (label or "(unlabeled)"),
        "observed bins (n)     : %s" % "{:,}".format(int(n_bins)),
        "pooled null bins (n)  : %s" % "{:,}".format(int(n_null)),
        "permutation floor p   : %.6e  (= 1 / (n_null + 1))" % p_floor,
        "genomic-inflation lam : %s" % lam_txt,
        "bins at floor         : %s (%.4f%%)" % (
            "{:,}".format(int(round(frac_floor * n_bins))), 100.0 * frac_floor),
        "BH-significant (FDR %g): %s (%.6f%%)" % (
            fdr, "{:,}".format(int(n_sig)), 100.0 * n_sig / n_bins),
        "KS vs uniform (stat)  : %.4f   [%s]" % (ks_stat, ks_caveat),
        "",
        "VERDICT: %s" % verdict,
        "",
        "Note: p-values here are the folded-magnitude --diff statistic from "
        "bearing_pvalue.py (P(|X| >= |s|)), evaluated through the same "
        "parse_qcat / collect_null_scores / empirical_pvals / bh_fdr functions. "
        "Read calibration from the bulk (panels A, B, D), not the saturated tail.",
    ]
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Within-condition replicate-differential calibration check "
                    "for the BEARING differential FDR pipeline. Reuses "
                    "bearing_pvalue.py; does not reimplement scoring.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff-qcat", required=True, metavar="BGZ",
                    help="Observed within-condition rep1-vs-rep2 diff qcat.bgz "
                         "(from compare_qcat.py with reps as single-replicate "
                         "conditions).")
    ap.add_argument("--null-qcat", required=True, nargs="+", metavar="BGZ",
                    help="One or more permuted rep-vs-rep diff qcat.bgz files "
                         "(from generate_perm_nulls.py --diff-sheet).")
    ap.add_argument("--out-prefix", required=True, metavar="PREFIX",
                    help="Output prefix; writes PREFIX_calibration.pdf, "
                         "PREFIX_perbin.tsv.gz, PREFIX_summary.txt.")
    ap.add_argument("--min-signal", type=float, default=0.0, metavar="FLOAT",
                    help="abs(score) threshold applied IDENTICALLY to observed "
                         "and null (default 0.0).")
    ap.add_argument("--fdr", type=float, default=0.05, metavar="FLOAT",
                    help="Benjamini-Hochberg FDR level for the significant-bin "
                         "count (default 0.05).")
    ap.add_argument("--bearing-pvalue", default=None, metavar="PATH",
                    help="Path to bearing_pvalue.py (file or its directory). "
                         "Defaults to this script's directory, then CWD.")
    ap.add_argument("--label", default="", metavar="STR",
                    help="Label for figure titles / summary (e.g. "
                         "'DN rep1 vs rep2').")
    ap.add_argument("--max-qq-points", type=int, default=50000, metavar="N",
                    help="Max points drawn in the QQ scatter (statistics use "
                         "all bins; default 50000).")
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    if out_prefix.parent and not out_prefix.parent.exists():
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 64)
    print("  bearing_calibration.py")
    print("  observed diff : %s" % args.diff_qcat)
    print("  null diff (%d) : %s" % (len(args.null_qcat), args.null_qcat[0]
                                     + (" ..." if len(args.null_qcat) > 1 else "")))
    print("  min-signal    : %g" % args.min_signal)
    print("  FDR           : %g" % args.fdr)
    print("=" * 64)

    bp = import_bearing_pvalue(args.bearing_pvalue)
    preflight(bp)

    print("Loading observed differential...", file=sys.stderr)
    rows, signed, tested = load_observed(bp, args.diff_qcat, args.min_signal)
    print("  observed bins: %s" % "{:,}".format(len(rows)), file=sys.stderr)

    print("Pooling permutation null...", file=sys.stderr)
    null_sorted, _ = bp.collect_null_scores(
        args.null_qcat, min_signal=args.min_signal, diff_mode=True)
    n_null = len(null_sorted)
    p_floor = 1.0 / (n_null + 1.0)

    # Identical statistic to bearing_pvalue.py --diff: P(|X| >= |s|).
    pvals = bp.empirical_pvals(tested, null_sorted)
    rejected, padj = bp.bh_fdr(pvals, args.fdr)
    n_sig = int(np.count_nonzero(rejected))

    lam, frac_floor, saturated = genomic_inflation(pvals, p_floor)

    # KS vs uniform is informational only: at genome scale it is hugely
    # over-powered and will reject for trivially small departures.
    ks_stat = float(stats.kstest(pvals, "uniform").statistic)
    ks_caveat = ("informational only; over-powered at n=%s, do not threshold on it"
                 % "{:,}".format(len(pvals)))

    pdf_path = str(out_prefix) + "_calibration.pdf"
    tsv_path = str(out_prefix) + "_perbin.tsv.gz"
    txt_path = str(out_prefix) + "_summary.txt"

    make_figure(pdf_path, args.label, pvals, tested, null_sorted, p_floor,
                lam, frac_floor, saturated, n_sig, args.fdr, args.max_qq_points)
    write_perbin_tsv(tsv_path, rows, signed, tested, pvals, padj, rejected,
                     args.fdr)
    verdict = build_verdict(lam, frac_floor, saturated, n_sig, len(rows),
                            args.fdr)
    summary_text = write_summary(txt_path, args.label, len(rows), n_null,
                                 p_floor, lam, frac_floor, saturated, n_sig,
                                 args.fdr, ks_stat, ks_caveat, verdict)

    print("\n" + summary_text)
    print("\nWrote:")
    for f in (pdf_path, tsv_path, txt_path):
        print("  %s" % f)


if __name__ == "__main__":
    main()
