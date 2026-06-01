#!/usr/bin/env python3
"""
bearing_pvalue.py -- Statistical significance scoring for Bearing qcat output.

Computes per-bin p-values for Bearing scores (sum of per-state KL contributions
per 200 bp bin) corrected for multiple testing (Benjamini-Hochberg FDR).

Outputs:
  1. A BigWig of -log10(p-value) per bin  (<prefix>.neglog10p.bw)
  2. A TSV with full per-bin statistics   (<prefix>.stats.tsv)
  3. Optionally a score distribution PDF  (<prefix>.score_dist.pdf)

NULL MODEL OPTIONS
------------------
Three null model strategies are available, in order of recommendation:

  1. --null-qcat FILE [FILE ...]   [RECOMMENDED]
     Empirical permutation null. Supply one or more qcat.bgz files produced
     by running bigwig_to_qcat.py on circularly-shifted input BigWigs (see
     below). P-values are empirical survival probabilities under the pooled
     null score distribution. This approach makes no parametric assumptions
     and is valid for any assay combination, any cell type, and any score
     distribution shape -- including the flat/uniform distributions observed
     in broadly active cell types such as 3T3 fibroblasts and DN thymocytes,
     where parametric Gamma-based methods fail.

  2. Default (no flag)
     2-component Gamma mixture EM. Appropriate when the score distribution
     is bimodal with a clear low-score background peak. Fails silently on
     flat distributions -- always check --score-plot.

  3. --fit-quantile FLOAT
     Fit Gamma to the lower FRACTION of bins. Simpler fallback that avoids
     the EM boundary problem but still assumes a Gamma shape. Use when the
     EM gives unexpected results but the distribution has a visible left peak.

HOW TO GENERATE PERMUTATION NULL BIGWIGS
-----------------------------------------
For each input BigWig, apply a random circular shift before running scoring:

  python shift_bigwig.py --bw atac.bw --shift 5000000 --out atac_perm.bw
  # (repeat for all tracks with independent random shifts)
  python bigwig_to_qcat.py --bw atac_perm.bw ctcf_perm.bw ... --out perm1.qcat.bgz

Then pass the permuted qcat(s) as the null:

  python bearing_pvalue.py \\
    --qcat sample.qcat.bgz \\
    --null-qcat perm1.qcat.bgz perm2.qcat.bgz \\
    --out-prefix results/sample \\
    --score-plot

Multiple permuted qcats can be pooled for a denser null distribution.
The more permutations, the lower the minimum achievable p-value:
  min_p = 1 / (n_null_bins + 1)

SIGNED DIFF QCAT
-----------------
For signed diff qcats (from compare_qcat.py), use --diff to test both
directions using |bearing_score|. The output BigWig is signed:
positive -log10(p) for A > B bins, negative for B > A bins.

CATEGORY NAMES
--------------
The companion *_cats.json produced by bigwig_to_qcat.py is auto-detected
from the qcat filename. Per-state TSV columns are then labelled by assay
name (e.g. kl_ATAC, kl_CTCF) instead of kl_1, kl_2, etc.

USAGE
-----
  # Recommended: empirical permutation null
  python bearing_pvalue.py \\
    --qcat sample.qcat.bgz \\
    --null-qcat perm1.qcat.bgz perm2.qcat.bgz \\
    --out-prefix results/sample \\
    --fdr 0.05 \\
    --score-plot

  # Parametric Gamma (EM -- check score plot)
  python bearing_pvalue.py \\
    --qcat sample.qcat.bgz \\
    --out-prefix results/sample \\
    --score-plot

  # Parametric Gamma (manual quantile fallback)
  python bearing_pvalue.py \\
    --qcat sample.qcat.bgz \\
    --out-prefix results/sample \\
    --fit-quantile 0.40

  # Signed diff qcat with permutation null
  python bearing_pvalue.py \\
    --qcat diff_A_vs_B.qcat.bgz \\
    --null-qcat diff_perm1.qcat.bgz \\
    --out-prefix results/diff_A_vs_B \\
    --diff

OUTPUT COLUMNS (TSV)
--------------------
  chrom, start, end       -- 200 bp bin coordinates
  bearing_score           -- sum of KL contributions across all states
  null_method             -- empirical, gamma_em, or gamma_quantile
  pval                    -- one-sided p-value: P(X >= score) under null
  pval_adj_bh             -- Benjamini-Hochberg adjusted p-value
  significant_fdrN        -- 1 if pval_adj_bh <= FDR threshold, else 0
  kl_<name>               -- per-state KL score (one column per assay)

In --diff mode, additional columns:
  bearing_score_tested    -- |bearing_score| used for null test
  direction               -- +, -, or 0

DEPENDENCIES
------------
  Python 3.8+
  pip install numpy scipy pyBigWig
  pip install matplotlib   # only for --score-plot
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyBigWig


# ---------------------------------------------------------------------------
# Built-in mm10 chromosome sizes (GRCm38)
# ---------------------------------------------------------------------------
MM10_CHROM_SIZES = {
    "chr1":  195471971, "chr2":  182113224, "chr3":  160039680,
    "chr4":  156508116, "chr5":  151834684, "chr6":  149736546,
    "chr7":  145441459, "chr8":  129401213, "chr9":  124595110,
    "chr10": 130694993, "chr11": 122082543, "chr12": 120129022,
    "chr13": 120421639, "chr14": 124902244, "chr15": 104043685,
    "chr16":  98207768, "chr17":  94987271, "chr18":  90702639,
    "chr19":  61431566, "chrX":  171031299, "chrY":   91744698,
}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_qcat(path, min_signal=0.0, diff_mode=False):
    """
    Yield (chrom, start, end, score_total, per_track_dict) from a qcat.bgz.

    In diff_mode, bins are filtered by abs(score) >= min_signal.
    Returns score_total with its original sign; filtering uses abs value.
    """
    import gzip
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            qcat_col = parts[3]
            if qcat_col.startswith("{"):
                payload = json.loads(qcat_col)
                pairs = payload.get("qcat", [])
            else:
                qcat_start = qcat_col.find("qcat:")
                if qcat_start == -1:
                    continue
                raw_start = qcat_col.find(",raw:", qcat_start)
                if raw_start >= 0:
                    qcat_payload = qcat_col[qcat_start + 5:raw_start]
                else:
                    qcat_payload = qcat_col[qcat_start + 5:]
                pairs = json.loads(qcat_payload)
            per_track = {int(state_idx): float(score)
                         for score, state_idx in pairs}
            score_total = sum(per_track.values())
            tested = abs(score_total) if diff_mode else score_total
            if tested < min_signal:
                continue
            yield chrom, start, end, score_total, per_track


def collect_null_scores(null_paths, min_signal=0.0, diff_mode=False,
                        max_per_file=None, per_track_pvals=False):
    """
    Pool scores from one or more permuted null qcat files.

    Returns a sorted numpy array of null scores (absolute values if diff_mode).
    Subsamples to max_per_file per file if specified, to manage memory.
    """
    null_scores = []
    null_per_track = defaultdict(list) if per_track_pvals else None
    for path in null_paths:
        file_scores = []
        print(f"  Loading null: {path}", file=sys.stderr)
        for _, _, _, score, per_track in parse_qcat(path, min_signal=min_signal,
                                                    diff_mode=diff_mode):
            s = abs(score) if diff_mode else score
            if s >= min_signal:
                file_scores.append(s)
                if per_track_pvals:
                    for t_idx, t_score in per_track.items():
                        null_per_track[int(t_idx)].append(float(abs(t_score) if diff_mode else t_score))
        print(f"    {len(file_scores):,} null bins loaded", file=sys.stderr)
        if max_per_file and len(file_scores) > max_per_file:
            rng = np.random.default_rng(42)
            file_scores = rng.choice(file_scores, size=max_per_file,
                                     replace=False).tolist()
            print(f"    Subsampled to {max_per_file:,}", file=sys.stderr)
        null_scores.extend(file_scores)

    if not null_scores:
        sys.exit("ERROR: no null scores collected from --null-qcat files.")

    arr = np.sort(np.array(null_scores, dtype=np.float64))
    print(f"  Total null bins: {len(arr):,}  "
          f"(min p achievable: {1/(len(arr)+1):.2e})", file=sys.stderr)
    if not per_track_pvals:
        return arr, None

    per_track_sorted = {}
    for t_idx, vals in null_per_track.items():
        if vals:
            per_track_sorted[t_idx] = np.sort(np.array(vals, dtype=np.float64))
    return arr, per_track_sorted


def empirical_pvals(observed_scores, null_scores_sorted):
    """
    Compute empirical one-sided p-values: P(X >= s) under null.

    Uses the sorted null array for fast searchsorted lookup.

    p(s) = (number of null scores >= s) / (total null scores + 1)

    The +1 in the denominator (Davison-Hinkley convention) ensures
    p > 0 always and is conservative for scores above the null maximum.
    """
    n_null = len(null_scores_sorted)
    # searchsorted gives the index of the first null score >= s
    idx = n_null - np.searchsorted(null_scores_sorted, observed_scores,
                                   side="left")
    pvals = (idx + 1) / (n_null + 1)  # +1 numerator: count s itself
    # Conservative clip: p cannot exceed 1, cannot be exactly 0
    return np.clip(pvals, 1.0 / (n_null + 1), 1.0)


def empirical_pvals_per_track(observed_per_track_matrix, null_per_track_dict):
    """
    Compute empirical p-values per track independently.

    observed_per_track_matrix : (n_bins, n_tracks)
    null_per_track_dict       : {track_idx (1-based): sorted null scores}
    """
    n_bins, n_tracks = observed_per_track_matrix.shape
    out = np.ones((n_bins, n_tracks), dtype=np.float64)
    for col in range(n_tracks):
        t_idx = col + 1
        null_sorted = null_per_track_dict.get(t_idx)
        if null_sorted is None or len(null_sorted) == 0:
            continue
        out[:, col] = empirical_pvals(observed_per_track_matrix[:, col], null_sorted)
    return out


def self_null_pvals(score_matrix):
    """
    Per-track self-null p-values from within-track genome-wide ranks.

    p(b,i) = (n_bins - rank_i + 1) / (n_bins + 1), where rank_i is ascending
    rank (1-based) of score s(b,i) within track i.
    """
    n_bins, n_tracks = score_matrix.shape
    out = np.ones((n_bins, n_tracks), dtype=np.float64)
    for i in range(n_tracks):
        col = score_matrix[:, i]
        order = np.argsort(col, kind="mergesort")
        ranks = np.empty(n_bins, dtype=np.int64)
        ranks[order] = np.arange(1, n_bins + 1)
        out[:, i] = (n_bins - ranks + 1) / (n_bins + 1)
    return out


def sanitize_track_name(name):
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")
    return s or "track"


def write_pvalue_bigwig(out_path, chrom_sizes, rows, pvals, diff_mode=False):
    """Write -log10(p) BigWig for provided p-values aligned to rows."""
    neglog10_p = -np.log10(np.clip(np.asarray(pvals, dtype=np.float64), 1e-300, 1.0))

    chrom_list = list(chrom_sizes.items())
    bw = pyBigWig.open(out_path, "w")
    bw.addHeader(chrom_list)

    by_chrom = defaultdict(list)
    for i, (chrom, start, end, score, _) in enumerate(rows):
        bw_val = float(neglog10_p[i])
        if diff_mode and score < 0:
            bw_val = -bw_val
        by_chrom[chrom].append((start, end, bw_val))

    skipped_unknown = skipped_invalid = skipped_nonmono = 0
    for chrom, chrom_size in chrom_list:
        if chrom not in by_chrom:
            continue
        if chrom_size is None:
            skipped_unknown += len(by_chrom[chrom])
            continue
        entries = sorted(by_chrom[chrom])
        starts, ends, vals = [], [], []
        last_start = -1
        for s, e, v in entries:
            s = max(0, int(s))
            e = min(int(e), int(chrom_size))
            if e <= s:
                skipped_invalid += 1
                continue
            if s <= last_start:
                skipped_nonmono += 1
                continue
            starts.append(s)
            ends.append(e)
            vals.append(float(v))
            last_start = s
        if starts:
            bw.addEntries([chrom] * len(starts), starts, ends=ends, values=vals)
    bw.close()
    return skipped_unknown, skipped_invalid, skipped_nonmono


# ---------------------------------------------------------------------------
# Parametric Gamma methods (kept as fallback)
# ---------------------------------------------------------------------------

def fit_gamma(scores, fit_quantile=1.0):
    """Fit Gamma to the lower fit_quantile fraction (manual override)."""
    from scipy.stats import gamma as gamma_dist
    arr = np.array(scores, dtype=np.float64)
    arr = arr[arr > 0]
    if fit_quantile < 1.0:
        cap = np.quantile(arr, fit_quantile)
        arr = arr[arr <= cap]
    mean, var = arr.mean(), arr.var(ddof=1)
    k, theta = mean**2 / var, var / mean
    return gamma_dist(a=k, scale=theta), k, theta, float(arr.max())


def fit_gamma_mixture(scores, max_iter=300, tol=1e-6):
    """
    Fit a 2-component Gamma mixture by EM to isolate the background null.

    WARNING: fails silently on flat/uniform distributions (broadly active
    cell types). Always check --score-plot. Use --null-qcat instead when
    the score distribution lacks a visible low-score background peak.
    """
    from scipy.stats import gamma as gamma_dist
    arr = np.array(scores, dtype=np.float64)
    arr = arr[arr > 0]

    def _wmom(x, w):
        W = w.sum()
        if W < 10:
            return None, None
        m1 = np.dot(w, x) / W
        m2 = np.dot(w, x**2) / W
        var = m2 - m1**2
        if var <= 0 or m1 <= 0:
            return None, None
        return m1**2 / var, var / m1

    med = np.median(arr)
    k_bg, th_bg = _wmom(arr, (arr <= med).astype(float))
    k_sg, th_sg = _wmom(arr, (arr >  med).astype(float))
    if k_bg is None: k_bg, th_bg = 5.0, 0.2
    if k_sg is None: k_sg, th_sg = 2.0, 0.5
    pi_bg = 0.6

    prev_ll = -np.inf
    for it in range(1, max_iter + 1):
        pdf_bg = pi_bg * gamma_dist.pdf(arr, a=k_bg, scale=th_bg)
        pdf_sg = (1.0 - pi_bg) * gamma_dist.pdf(arr, a=k_sg, scale=th_sg)
        denom = pdf_bg + pdf_sg
        valid = denom > 0
        w_bg = np.where(valid, pdf_bg / denom, 0.5)
        ll = np.log(np.where(valid, denom, 1e-300)).sum()
        pi_bg = w_bg.mean()
        k_n, th_n = _wmom(arr, w_bg)
        k_s, th_s = _wmom(arr, 1.0 - w_bg)
        if k_n is not None: k_bg, th_bg = k_n, th_n
        if k_s is not None: k_sg, th_sg = k_s, th_s
        if k_bg * th_bg > k_sg * th_sg:
            k_bg, th_bg, k_sg, th_sg = k_sg, th_sg, k_bg, th_bg
            pi_bg = 1.0 - pi_bg
        if abs(ll - prev_ll) < tol * (abs(prev_ll) + 1.0):
            print(f"  Mixture EM converged after {it} iterations (ll={ll:.1f})",
                  file=sys.stderr)
            break
        prev_ll = ll
    else:
        print(f"  WARNING: EM reached max_iter={max_iter}.",
              file=sys.stderr)

    if pi_bg < 0.3 or pi_bg > 0.97:
        print(f"  WARNING: pi_bg={pi_bg:.3f} -- mixture may not have "
              f"converged cleanly. Consider --null-qcat (empirical null) "
              f"or --fit-quantile as a fallback.", file=sys.stderr)

    print(f"  Background: k={k_bg:.4f}, theta={th_bg:.4f}, "
          f"mean={k_bg*th_bg:.4f}, pi={pi_bg:.3f}", file=sys.stderr)
    print(f"  Signal:     k={k_sg:.4f}, theta={th_sg:.4f}, "
          f"mean={k_sg*th_sg:.4f}, pi={1-pi_bg:.3f}", file=sys.stderr)

    fitted = gamma_dist(a=k_bg, scale=th_bg)
    bg_cap = float(fitted.ppf(0.999))
    return fitted, k_bg, th_bg, bg_cap


# ---------------------------------------------------------------------------
# BH FDR correction
# ---------------------------------------------------------------------------

def bh_fdr(pvals, alpha):
    """
    Benjamini-Hochberg FDR correction.
    Returns (rejected, pvals_adj) arrays in input order.
    """
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    pvals_adj = np.minimum(1.0, pvals * n / ranks)
    # enforce monotonicity (cumulative minimum from right)
    pvals_adj_sorted = pvals_adj[order]
    for i in range(n - 2, -1, -1):
        pvals_adj_sorted[i] = min(pvals_adj_sorted[i], pvals_adj_sorted[i + 1])
    pvals_adj[order] = pvals_adj_sorted
    return pvals_adj <= alpha, pvals_adj


# ---------------------------------------------------------------------------
# Score distribution plot
# ---------------------------------------------------------------------------

def plot_score_distribution(scores_arr, out_path, sample_name="",
                             null_method="gamma_em",
                             score_normalised=False,
                             # parametric args (None for empirical)
                             fitted_dist=None, background_cap=None,
                             fdr_score=None,
                             # empirical null args
                             null_scores_sorted=None,
                             fdr_threshold_score=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib not available, skipping score plot",
              file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    cap = float(np.percentile(scores_arr, 99.9))
    ax.hist(scores_arr[scores_arr <= cap], bins=200, density=True,
            color="#aac4e0", alpha=0.70, label="Observed scores", zorder=2)

    xs = np.linspace(0, cap, 400)

    if null_method == "empirical" and null_scores_sorted is not None:
        # Plot the smoothed empirical null density
        null_cap = float(np.percentile(null_scores_sorted,
                                       min(99.9, 100 * len(null_scores_sorted)
                                           / (len(null_scores_sorted) + 1))))
        ax.hist(null_scores_sorted[null_scores_sorted <= null_cap],
                bins=200, density=True,
                color="#f4a261", alpha=0.45,
                label="Permutation null", zorder=3)
        # Mark FDR threshold
        if fdr_threshold_score is not None:
            ax.axvline(fdr_threshold_score, color="#2ca02c", lw=1.5,
                       label=f"FDR threshold ({fdr_threshold_score:.2f})")
        null_label = "Empirical permutation null"
    else:
        # Parametric Gamma overlay
        if fitted_dist is not None:
            ax.plot(xs, fitted_dist.pdf(xs), color="#e05c3a", lw=1.8,
                    label="Fitted Gamma null", zorder=4)
        if background_cap is not None:
            ax.axvline(background_cap, color="#aaa", lw=1.0, ls="--",
                       label=f"Background cap ({background_cap:.2f})")
        if fdr_score is not None:
            ax.axvline(fdr_score, color="#2ca02c", lw=1.5,
                       label=f"FDR threshold ({fdr_score:.2f})")
        null_label = ("Gamma mixture EM" if null_method == "gamma_em"
                      else f"Gamma (quantile fit)")

    ax.set_yscale("log")
    ax.set_xlabel("Bearing score (sum of KL contributions)")
    ax.set_ylabel("Density (log scale)")
    norm_note = " [0-1 normalised]" if score_normalised else ""
    title = f"Bearing score distribution vs null model  [{sample_name}]{norm_note}"
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
    ax.text(0.98, 0.04, f"Null: {null_label}",
            transform=ax.transAxes, fontsize=7,
            ha="right", va="bottom", color="#666")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Score plot: {out_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    ap.add_argument("--qcat",         required=True,
                    help="Input qcat.bgz produced by bigwig_to_qcat.py.")
    ap.add_argument("--null-qcat",    nargs="+", default=None,
                    metavar="FILE",
                    help="[RECOMMENDED] One or more qcat.bgz files from "
                         "circularly-shifted (permuted) BigWigs. Enables "
                         "the empirical permutation null. Multiple files "
                         "are pooled for a denser null distribution.")
    ap.add_argument("--cats-json",    default=None,
                    help="Categories JSON from bigwig_to_qcat.py. "
                         "Auto-detected from the qcat path if not given.")
    ap.add_argument("--chrom-sizes",  default=None,
                    help="Two-column chrom.sizes file. "
                         "Defaults to built-in mm10 sizes.")
    ap.add_argument("--out-prefix",   required=True,
                    help="Output path prefix.")
    ap.add_argument("--fdr",          type=float, default=0.05,
                    help="FDR threshold for BH correction (default 0.05).")
    ap.add_argument("--min-signal",   type=float, default=0.5,
                    help="Minimum score (or |score| in --diff mode) to "
                         "include in null fitting and significance testing "
                         "(default 0.5).")
    ap.add_argument("--fit-quantile", type=float, default=None,
                    metavar="FLOAT",
                    help="Parametric fallback: fit the Gamma null to the "
                         "lower FRACTION of bins. E.g. --fit-quantile 0.40. "
                         "Ignored if --null-qcat is given.")
    ap.add_argument("--score-plot",   action="store_true",
                    help="Save a score distribution plot "
                         "(<out-prefix>.score_dist.pdf).")
    ap.add_argument("--diff",         action="store_true",
                    help="Diff mode: filter and test using |bearing_score|. "
                         "Output BigWig is signed: positive -log10(p) for "
                         "A > B bins, negative for B > A bins.")
    ap.add_argument("--null-subsample", type=int, default=None,
                    metavar="N",
                    help="Subsample at most N bins per null qcat file "
                         "to manage memory (default: use all bins).")
    ap.add_argument("--per-track-pvals", action="store_true",
                    help=(
                        "In addition to the overall p-value, compute a separate p-value for "
                        "each track i testing whether s(b,i) is unusual under the permutation "
                        "null. Requires --null-qcat. Adds N extra columns to the TSV "
                        "(pval_<track> and pval_adj_<track>) and writes N extra BigWig files "
                        "(<prefix>.<track>.neglog10p.bw)."
                    ))
    ap.add_argument("--self-null", action="store_true",
                    help=(
                        "Compute per-track p-values using the genome-wide distribution of "
                        "each track's own KL scores as the null (no permutation files needed). "
                        "p(b,i) = rank(s(b,i)) / n_bins. Produces per-track p-values only; "
                        "the overall S(b) p-value still requires --null-qcat or Gamma null. "
                        "Caveat: tests 'is this bin unusual for track i genome-wide' not "
                        "'is this bin unusual relative to random co-enrichment'."
                    ))
    ap.add_argument("--sort-output", action="store_true",
                    help=(
                        "Write the stats TSV sorted by ascending pval_adj_bh (most "
                        "significant first) rather than genomic order."
                    ))
    ap.add_argument("--sig-bed", metavar="FILE", default=None,
                    help=(
                        "Write a BED4 file of significant bins (pval_adj_bh <= --fdr) sorted "
                        "by ascending pval_adj_bh. Columns: chrom, start, end, -log10(pval_adj). "
                        "Suitable for use as a ranked feature list."
                    ))
    args = ap.parse_args()

    # ── 0. Determine null method ───────────────────────────────────────────
    if args.null_qcat:
        null_method = "empirical"
    elif args.fit_quantile is not None:
        null_method = "gamma_quantile"
    else:
        null_method = "gamma_em"

    print(f"Null method: {null_method}", file=sys.stderr)
    if null_method == "gamma_em":
        print("  NOTE: Gamma EM can fail on flat distributions (broadly "
              "active cell types). Always check --score-plot. Consider "
              "--null-qcat for robustness.", file=sys.stderr)

    # ── 0b. Load category names ────────────────────────────────────────────
    cats_json_path = args.cats_json
    if cats_json_path is None:
        base = str(args.qcat)
        for suffix in (".qcat.bgz", ".bgz"):
            if base.endswith(suffix):
                base = base[:-len(suffix)]
                break
        candidate = base + "_cats.json"
        if Path(candidate).exists():
            cats_json_path = candidate
            print(f"Auto-detected categories: {cats_json_path}", file=sys.stderr)

    state_names: dict = {}
    if cats_json_path:
        with open(cats_json_path) as fh:
            cats_doc = json.load(fh)
        cats = cats_doc.get("categories", {})
        for key, val in cats.items():
            name = val[0] if isinstance(val, list) else val.get("name", key)
            state_names[int(key)] = name
        score_normalised = bool(cats_doc.get("normalize_score", False))
        print(f"Loaded {len(state_names)} category names.", file=sys.stderr)
    else:
        score_normalised = False
        print("No categories JSON found -- TSV columns will be kl_1, kl_2, ...",
              file=sys.stderr)

    # ── 1. Load chrom sizes ────────────────────────────────────────────────
    if args.chrom_sizes:
        chrom_sizes = {}
        with open(args.chrom_sizes) as fh:
            for line in fh:
                c, s = line.split()
                chrom_sizes[c] = int(s)
    else:
        chrom_sizes = MM10_CHROM_SIZES
        print("Using built-in mm10 chrom sizes.", file=sys.stderr)

    # ── 2. Collect observed scores ─────────────────────────────────────────
    print("Pass 1: collecting observed scores...", file=sys.stderr)
    all_bins = []
    all_scores = []
    n_pos = n_neg = n_zero = 0

    for chrom, start, end, score, per_track in parse_qcat(
            args.qcat, min_signal=args.min_signal, diff_mode=args.diff):
        tested = abs(score) if args.diff else score
        all_scores.append(tested)
        all_bins.append((chrom, start, end, score, per_track))
        if score > 0: n_pos += 1
        elif score < 0: n_neg += 1
        else: n_zero += 1

    if args.diff:
        print(f"  Diff mode: {len(all_bins):,} bins kept "
              f"(+:{n_pos:,}  -:{n_neg:,}  0:{n_zero:,})",
              file=sys.stderr)

    if len(all_scores) < 100:
        sys.exit("ERROR: too few bins with signal to compute p-values.")

    scores_arr = np.array(all_scores, dtype=np.float64)
    print(f"  Observed bins: {len(scores_arr):,}  "
          f"mean={scores_arr.mean():.4f}  "
          f"p50={np.percentile(scores_arr,50):.4f}  "
          f"p99={np.percentile(scores_arr,99):.4f}", file=sys.stderr)

    # Build aligned per-track score matrix for optional per-track p-values.
    track_keys = sorted({k for _, _, _, _, per_track in all_bins for k in per_track.keys()})
    per_track_obs = None
    if track_keys:
        per_track_obs = np.array(
            [[float(per_track.get(t, 0.0)) for t in track_keys] for _, _, _, _, per_track in all_bins],
            dtype=np.float64,
        )

    # ── 3. Fit null model and compute p-values ─────────────────────────────
    fitted_dist = None
    background_cap = None
    k = theta = None
    null_scores_sorted = None
    fdr_threshold_score = None

    if null_method == "empirical":
        # ── Empirical permutation null ─────────────────────────────────────
        print("Loading permutation null scores...", file=sys.stderr)
        null_scores_sorted, null_per_track_sorted = collect_null_scores(
            args.null_qcat,
            min_signal=args.min_signal,
            diff_mode=args.diff,
            max_per_file=args.null_subsample,
            per_track_pvals=args.per_track_pvals,
        )
        pvals = empirical_pvals(scores_arr, null_scores_sorted)

        # Diagnostics
        print("  Score distribution (empirical p-values):", file=sys.stderr)
        for pct in [50, 75, 90, 95, 99, 99.9]:
            s = float(np.percentile(scores_arr, pct))
            p = float(empirical_pvals(np.array([s]), null_scores_sorted)[0])
            print(f"    p{pct:5.1f} = {s:.4f}  ->  "
                  f"p = {p:.2e}  (-log10p = {-np.log10(p+1e-300):.2f})",
                  file=sys.stderr)

        # BH threshold score (approximate: find score where empirical p
        # equals fdr/n_tests)
        n = len(scores_arr)
        bh_p_threshold = args.fdr / n
        # Find the score at which empirical p would equal bh_p_threshold
        # p(s) = (n_null >= s + 1) / (n_null_total + 1)
        # => n_null >= s = bh_p_threshold * (n_null_total + 1) - 1
        n_null_total = len(null_scores_sorted)
        target_count = max(0, int(np.ceil(
            bh_p_threshold * (n_null_total + 1) - 1)))
        if target_count < n_null_total:
            fdr_threshold_score = float(
                null_scores_sorted[n_null_total - target_count - 1])
        else:
            fdr_threshold_score = float(null_scores_sorted[-1])
        print(f"  For FDR={args.fdr} with {n:,} tests, "
              f"need p < {bh_p_threshold:.2e} "
              f"-> score > {fdr_threshold_score:.4f}", file=sys.stderr)

    else:
        # ── Parametric Gamma null ──────────────────────────────────────────
        if null_method == "gamma_quantile":
            print(f"Fitting Gamma to lower {args.fit_quantile:.0%} of bins...",
                  file=sys.stderr)
            fitted_dist, k, theta, background_cap = fit_gamma(
                all_scores, fit_quantile=args.fit_quantile)
            n_bg = int((scores_arr <= background_cap).sum())
            print(f"  k={k:.4f}, theta={theta:.4f}  "
                  f"background cap={background_cap:.4f}  "
                  f"({n_bg:,}/{len(scores_arr):,} background bins)",
                  file=sys.stderr)
        else:
            print("Fitting 2-component Gamma mixture (EM)...", file=sys.stderr)
            fitted_dist, k, theta, background_cap = fit_gamma_mixture(
                all_scores)

        from scipy.stats import gamma as gamma_dist_scipy
        pvals = 1.0 - fitted_dist.cdf(scores_arr)
        pvals = np.clip(pvals, 1e-300, 1.0)
        null_per_track_sorted = None

        print("  Score distribution:", file=sys.stderr)
        for pct in [50, 75, 90, 95, 99, 99.9]:
            s = float(np.percentile(scores_arr, pct))
            p = float(1.0 - fitted_dist.cdf(s))
            print(f"    p{pct:5.1f} = {s:.4f}  ->  "
                  f"p = {p:.2e}  (-log10p = {-np.log10(p+1e-300):.2f})",
                  file=sys.stderr)
        n = len(scores_arr)
        bh_p_threshold = args.fdr / n
        fdr_threshold_score = float(fitted_dist.ppf(1.0 - bh_p_threshold))
        print(f"  For FDR={args.fdr} with {n:,} tests, "
              f"need p < {bh_p_threshold:.2e} "
              f"-> score > {fdr_threshold_score:.4f}", file=sys.stderr)
        if args.diff:
            print("  Diff mode threshold applies to |bearing_score|.",
                  file=sys.stderr)

    # ── 4. Score distribution plot ─────────────────────────────────────────
    if args.score_plot:
        plot_score_distribution(
            scores_arr,
            out_path=args.out_prefix + ".score_dist.pdf",
            sample_name=Path(args.out_prefix).name,
            null_method=null_method,
            score_normalised=score_normalised,
            fitted_dist=fitted_dist,
            background_cap=background_cap,
            fdr_score=fdr_threshold_score if null_method != "empirical" else None,
            null_scores_sorted=null_scores_sorted,
            fdr_threshold_score=fdr_threshold_score if null_method == "empirical" else None,
        )

    # ── 5. BH FDR correction ───────────────────────────────────────────────
    rejected, pvals_adj = bh_fdr(pvals, alpha=args.fdr)
    neglog10_p = -np.log10(np.clip(pvals, 1e-300, 1.0))

    # Optional per-track p-values
    perm_track_enabled = bool(args.per_track_pvals and args.null_qcat and per_track_obs is not None)
    if args.per_track_pvals and not args.null_qcat:
        print("WARNING: --per-track-pvals requested without --null-qcat; skipping permutation per-track p-values.",
              file=sys.stderr)
    self_track_enabled = bool(args.self_null and per_track_obs is not None)

    perm_track_p = perm_track_adj = None
    self_track_p = self_track_adj = None

    if perm_track_enabled:
        perm_track_p = empirical_pvals_per_track(per_track_obs, null_per_track_sorted or {})
        perm_track_adj = np.zeros_like(perm_track_p)
        for j in range(perm_track_p.shape[1]):
            _, adj = bh_fdr(perm_track_p[:, j], alpha=args.fdr)
            perm_track_adj[:, j] = adj

    if self_track_enabled:
        self_track_p = self_null_pvals(np.abs(per_track_obs) if args.diff else per_track_obs)
        self_track_adj = np.zeros_like(self_track_p)
        for j in range(self_track_p.shape[1]):
            _, adj = bh_fdr(self_track_p[:, j], alpha=args.fdr)
            self_track_adj[:, j] = adj

    # Reflect self-null mode in TSV metadata while preserving the overall null
    # method context when both permutation and self-null are present.
    if self_track_enabled and perm_track_enabled:
        tsv_null_method = f"{null_method}+self_null"
    elif self_track_enabled:
        tsv_null_method = "self_null"
    else:
        tsv_null_method = null_method

    # ── 6. Write BigWig ────────────────────────────────────────────────────
    bw_path = args.out_prefix + ".neglog10p.bw"
    print(f"Writing BigWig: {bw_path}", file=sys.stderr)
    skipped_unknown, skipped_invalid, skipped_nonmono = write_pvalue_bigwig(
        bw_path,
        chrom_sizes,
        all_bins,
        pvals,
        diff_mode=args.diff,
    )

    if any([skipped_unknown, skipped_invalid, skipped_nonmono]):
        print(f"  WARNING: skipped {skipped_unknown} unknown-chrom, "
              f"{skipped_invalid} invalid, {skipped_nonmono} non-monotonic bins.",
              file=sys.stderr)

    # Write optional per-track BigWigs
    if track_keys and (perm_track_enabled or self_track_enabled):
        for col, t_idx in enumerate(track_keys):
            t_name = state_names.get(t_idx, str(t_idx))
            safe_t = sanitize_track_name(t_name)

            if perm_track_enabled:
                perm_bw = (
                    f"{args.out_prefix}.perm.{safe_t}.neglog10p.bw"
                    if self_track_enabled else
                    f"{args.out_prefix}.{safe_t}.neglog10p.bw"
                )
                print(f"Writing per-track BigWig: {perm_bw}", file=sys.stderr)
                write_pvalue_bigwig(
                    perm_bw,
                    chrom_sizes,
                    all_bins,
                    perm_track_p[:, col],
                    diff_mode=args.diff,
                )

            if self_track_enabled:
                self_bw = (
                    f"{args.out_prefix}.self.{safe_t}.neglog10p.bw"
                    if perm_track_enabled else
                    f"{args.out_prefix}.{safe_t}.neglog10p.bw"
                )
                print(f"Writing self-null per-track BigWig: {self_bw}", file=sys.stderr)
                write_pvalue_bigwig(
                    self_bw,
                    chrom_sizes,
                    all_bins,
                    self_track_p[:, col],
                    diff_mode=args.diff,
                )

    # ── 7. Write TSV ───────────────────────────────────────────────────────
    tsv_path = args.out_prefix + ".stats.tsv"
    print(f"Writing TSV: {tsv_path}", file=sys.stderr)

    def _col(idx):
        return "kl_" + state_names.get(idx, str(idx))

    def _track_label(idx):
        return state_names.get(idx, str(idx))

    perm_track_prefix = "pval_perm_" if (perm_track_enabled and self_track_enabled) else "pval_"
    self_track_prefix = "pval_self_" if perm_track_enabled else "pval_"

    with open(tsv_path, "w") as fh:
        header = (
            ["chrom", "start", "end", "bearing_score"]
            + (["bearing_score_tested"] if args.diff else [])
            + ["null_method"]
            + ["score_normalised"]
            + (["gamma_k", "gamma_theta"] if k is not None else [])
            + ["pval", "pval_adj_bh",
               "significant_fdr" + str(args.fdr)]
            + (["direction"] if args.diff else [])
            + [_col(t) for t in track_keys]
        )
        if perm_track_enabled:
            header += [f"{perm_track_prefix}{_track_label(t)}" for t in track_keys]
            header += [f"{perm_track_prefix}adj_{_track_label(t)}" for t in track_keys]
        if self_track_enabled:
            header += [f"{self_track_prefix}{_track_label(t)}" for t in track_keys]
            header += [f"{self_track_prefix}adj_{_track_label(t)}" for t in track_keys]
        fh.write("\t".join(header) + "\n")

        rows_for_write = []
        for i, (chrom, start, end, score, per_track) in enumerate(all_bins):
            row = (
                [chrom, start, end, f"{score:.6f}"]
                + ([f"{abs(score):.6f}"] if args.diff else [])
                + [tsv_null_method]
                + ["true" if score_normalised else "false"]
                + ([f"{k:.6f}", f"{theta:.6f}"] if k is not None else [])
                + [f"{pvals[i]:.6e}",
                   f"{pvals_adj[i]:.6e}",
                   "1" if rejected[i] else "0"]
                + (["+" if score > 0 else "-" if score < 0 else "0"]
                   if args.diff else [])
                + [f"{per_track.get(t, 0.0):.6f}" for t in track_keys]
            )
            if perm_track_enabled:
                row += [f"{perm_track_p[i, col]:.6e}" for col in range(len(track_keys))]
                row += [f"{perm_track_adj[i, col]:.6e}" for col in range(len(track_keys))]
            if self_track_enabled:
                row += [f"{self_track_p[i, col]:.6e}" for col in range(len(track_keys))]
                row += [f"{self_track_adj[i, col]:.6e}" for col in range(len(track_keys))]
            rows_for_write.append(row)

        if args.sort_output:
            rows_for_write.sort(key=lambda r: float(r[header.index("pval_adj_bh")]))

        for row in rows_for_write:
            fh.write("\t".join(str(x) for x in row) + "\n")

    if args.sig_bed:
        sig_rows = [
            (chrom, start, end, float(padj))
            for (chrom, start, end, _, _), padj, rej in zip(all_bins, pvals_adj, rejected)
            if bool(rej)
        ]
        sig_rows.sort(key=lambda x: x[3])
        with open(args.sig_bed, "w") as bf:
            for chrom, start, end, padj in sig_rows:
                bf.write(
                    f"{chrom}\t{start}\t{end}\t"
                    f"{-np.log10(max(padj, 1e-300)):.6f}\n"
                )

    # ── 8. Summary ─────────────────────────────────────────────────────────
    n_sig = int(rejected.sum())
    pct_sig = 100 * n_sig / len(all_bins) if all_bins else 0
    print(f"\nDone. {n_sig:,}/{len(all_bins):,} bins significant "
          f"({pct_sig:.1f}%)  FDR < {args.fdr}  [{null_method}]",
          file=sys.stderr)
    print(f"  BigWig : {bw_path}", file=sys.stderr)
    print(f"  TSV    : {tsv_path}", file=sys.stderr)
    if args.sig_bed:
        print(f"  SigBED : {args.sig_bed}", file=sys.stderr)


if __name__ == "__main__":
    main()
