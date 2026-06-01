#!/usr/bin/env python3
"""
evaluate_bearing_recovery.py
============================
Ground-truth evaluator for the BEARING pipeline.

Runs the REAL pipeline (bigwig_to_qcat.py) on synthetic BigWigs produced by
simulate_bearing_tracks.py, then checks that the recovered per-track allocation
matches the planted allocation, and that the recovered composition is stable
across the three free parameters (low-signal mask, zero-clamp epsilon, and
normalization quantile).

This evaluator does NOT reimplement scoring. The whole-genome runs go through
bigwig_to_qcat.py as a subprocess. The parameter sweep over the zero-clamp
epsilon reuses the REAL scoring functions (signals_to_prob, kl_scores_per_bin)
imported directly from bigwig_to_qcat.py, because that epsilon is the module
constant PSEUDOCOUNT and is not exposed as a CLI flag.

WHAT "RECOVERY" MEANS HERE (important)
--------------------------------------
BEARING's per-track score is a CLAMPED KL contribution:

    score_{b,i} = max(0, P_{b,i} * log2(P_{b,i} / Q_i))

The clamp to zero makes each per-track score a one-sided ENRICHMENT metric, not
a true KL divergence: tracks whose planted proportion is BELOW background
(a_i < Q_i) are suppressed to zero and carry no recovered mass. Recovery is
therefore evaluated on the ENRICHED SUBSET of each block (tracks with a_i > Q_i,
i.e. the tracks the block actually pushes up), and the suppressed tracks are
reported separately as expected-zero rather than counted as error. This is a
property of the metric, not a failure of recovery, and the asymmetry is shown
explicitly in its own panel.

Because the recovered allocation lives on the probability simplex (it sums to
1, so the per-track values are not independent), raw Pearson/Spearman/OLS on
closed compositional data can be misleading. We report those (reviewers expect
them) AND an Aitchison-style compositional error (centered-log-ratio distance
on the enriched subset).

INPUTS
    --truth        <prefix>_truth.tsv      from simulate_bearing_tracks.py
    --bw           six BigWigs in track order ATAC RNA+ RNA- CTCF RAD21 H3K27ac
    --chrom-sizes  <prefix>.chrom.sizes
    --scorer       path to bigwig_to_qcat.py (default: ./bigwig_to_qcat.py)

OUTPUTS (under --outdir)
    fig_a_recovered_vs_planted.png        panel a: pooled per-track scatter
    fig_b_named_blocks.png                panel b: planted vs recovered bars
    fig_d_null_vs_planted.png             panel d: score + dominant-call null behavior
    fig_e_parameter_sensitivity.png       panel e: stability across mask/eps/quantile
    recovery_per_block.tsv                planted vs recovered allocations + errors
    recovery_summary.txt                  metrics (Pearson, Spearman, slope,
                                          Aitchison, rank-recovery, null stats)

DEPENDENCIES
    numpy, scipy, matplotlib, pysam (for reading qcat.bgz), and the scorer's
    own deps (pyBigWig, pysam). matplotlib uses the Agg backend.
"""

import argparse
import importlib.util
import json
import os
import subprocess
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Track order MUST match simulate_bearing_tracks.py and the first six categories
# of bigwig_to_qcat.py (state index 1..6). RAD21 is "Cohesin" in the scorer.
TRACK_NAMES = ["ATAC", "RNA+", "RNA-", "CTCF", "RAD21", "H3K27ac"]
N_TRACKS = 6
TRUTH_A_COLS = ["a_ATAC", "a_RNAplus", "a_RNAminus", "a_CTCF", "a_RAD21",
                "a_H3K27ac"]


# -----------------------------------------------------------------------------
# Preflight: load the scorer and verify its interface matches our expectations.
# Fail loudly with a clear message if anything has drifted.
# -----------------------------------------------------------------------------
def load_scorer(scorer_path):
    if not os.path.isfile(scorer_path):
        sys.exit("ERROR: scorer not found at %s (use --scorer)" % scorer_path)
    spec = importlib.util.spec_from_file_location("b2q_scorer", scorer_path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as exc:
        sys.exit("ERROR: failed to import scorer %s: %r" % (scorer_path, exc))
    return mod


def preflight(mod, scorer_path):
    """
    Verify the scorer exposes the interface this evaluator depends on. Returns a
    dict of the discovered constants. Exits with a diagnostic on mismatch.
    """
    problems = []

    # required functions for the in-process epsilon sweep
    for fn in ("signals_to_prob", "kl_scores_per_bin"):
        if not hasattr(mod, fn):
            problems.append("missing function: %s" % fn)

    # required constants
    consts = {}
    for c in ("PSEUDOCOUNT", "MIN_SIGNAL", "BIN_SIZE", "ALL_CATEGORIES"):
        if not hasattr(mod, c):
            problems.append("missing constant: %s" % c)
        else:
            consts[c] = getattr(mod, c)

    # category names: first six must map to our track order
    if "ALL_CATEGORIES" in consts:
        cat_names = [c[0] for c in consts["ALL_CATEGORIES"][:N_TRACKS]]
        expected = ["ATAC", "RNAseq +", "RNAseq -", "CTCF", "Cohesin", "H3K27ac"]
        if cat_names != expected:
            problems.append(
                "first six scorer categories %r != expected %r; track-to-state "
                "mapping would be wrong" % (cat_names, expected))

    # CLI signature check: scorer must accept the flags we pass
    try:
        help_txt = subprocess.run(
            [sys.executable, scorer_path, "--help"],
            capture_output=True, text=True, timeout=60).stdout
    except Exception as exc:
        help_txt = ""
        problems.append("could not run '%s --help': %r" % (scorer_path, exc))
    for flag in ("--bw", "--out", "--chrom-sizes", "--min-signal",
                 "--normalize-method", "--no-extras"):
        if help_txt and flag not in help_txt:
            problems.append("scorer --help does not mention required flag %s"
                            % flag)

    if problems:
        sys.stderr.write("PREFLIGHT FAILED -- scorer interface mismatch:\n")
        for p in problems:
            sys.stderr.write("  - %s\n" % p)
        sys.stderr.write(
            "\nThis evaluator was written against a bigwig_to_qcat.py that:\n"
            "  * exposes signals_to_prob() and kl_scores_per_bin()\n"
            "  * defines PSEUDOCOUNT, MIN_SIGNAL, BIN_SIZE, ALL_CATEGORIES\n"
            "  * has first six categories ATAC, RNAseq +, RNAseq -, CTCF, "
            "Cohesin, H3K27ac\n"
            "  * accepts --bw --out --chrom-sizes --min-signal "
            "--normalize-method --no-extras\n"
            "Update the evaluator (or the scorer) so they agree, then re-run.\n")
        sys.exit(2)

    print("preflight OK: scorer interface matches "
          "(PSEUDOCOUNT=%g, MIN_SIGNAL=%g, BIN_SIZE=%d)"
          % (consts["PSEUDOCOUNT"], consts["MIN_SIGNAL"], consts["BIN_SIZE"]))
    return consts


# -----------------------------------------------------------------------------
# Truth-table loading.
# -----------------------------------------------------------------------------
def load_truth(path):
    """
    Return a structured dict with arrays indexed by bin_index:
        chrom (str, single value), starts, ends, block_id, bin_class,
        strength, A (n_bins x 6), is_null (bool).
    """
    import csv
    starts, ends, block_id, bin_class, strength, is_null = [], [], [], [], [], []
    A = []
    chrom = None
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            chrom = r["chrom"]
            starts.append(int(r["start"]))
            ends.append(int(r["end"]))
            block_id.append(r["block_id"])
            bin_class.append(r["bin_class"])
            strength.append(float(r["strength_s"]))
            is_null.append(r["is_null"] == "1")
            A.append([float(r[c]) for c in TRUTH_A_COLS])
    return {
        "chrom": chrom,
        "starts": np.array(starts, dtype=np.int64),
        "ends": np.array(ends, dtype=np.int64),
        "block_id": np.array(block_id, dtype=object),
        "bin_class": np.array(bin_class, dtype=object),
        "strength": np.array(strength, dtype=np.float64),
        "A": np.array(A, dtype=np.float64),
        "is_null": np.array(is_null, dtype=bool),
    }


# -----------------------------------------------------------------------------
# Run the REAL scorer via subprocess, parse the qcat output into per-bin
# score vectors keyed by genomic start coordinate.
# -----------------------------------------------------------------------------
def run_scorer(scorer_path, bw_paths, chrom_sizes, out_path,
               min_signal=None, normalize_method=None, normalize_tracks=False,
               pseudocount=None):
    cmd = [sys.executable, scorer_path,
           "--bw"] + list(bw_paths) + [
           "--out", out_path,
           "--chrom-sizes", chrom_sizes,
           "--no-extras"]
    if min_signal is not None:
        cmd += ["--min-signal", str(min_signal)]
    if normalize_tracks:
        cmd += ["--normalize-tracks"]
    if normalize_method is not None:
        cmd += ["--normalize-method", normalize_method]
    if pseudocount is not None:
        cmd += ["--pseudocount", str(pseudocount)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout + "\n" + res.stderr + "\n")
        sys.exit("ERROR: scorer failed (cmd: %s)" % " ".join(cmd))
    return out_path


def parse_qcat(qcat_bgz, chrom, n_states=N_TRACKS):
    """
    Parse a qcat.bgz into a dict: start_coord -> score vector (n_states,) in
    track order (state index i -> position i-1). Uses pysam to read the bgzf.
    Lines look like:
      chrom  start  end  id:N,qcat:[[score,state],...],raw:[...]
    """
    import pysam
    scores_by_start = {}
    tb = pysam.TabixFile(qcat_bgz)
    for line in tb.fetch(chrom):
        parts = line.rstrip("\n").split("\t")
        start = int(parts[1])
        col = parts[3]
        # extract the qcat:[...] JSON array
        key = "qcat:"
        i = col.find(key)
        if i < 0:
            continue
        j = i + len(key)
        # the array is well-formed JSON ending at the matching bracket; find it
        depth = 0
        end = None
        for k in range(j, len(col)):
            if col[k] == "[":
                depth += 1
            elif col[k] == "]":
                depth -= 1
                if depth == 0:
                    end = k + 1
                    break
        if end is None:
            continue
        pairs = json.loads(col[j:end])
        vec = np.zeros(n_states, dtype=np.float64)
        for score, state in pairs:
            vec[int(state) - 1] = float(score)
        scores_by_start[start] = vec
    return scores_by_start


def scores_matrix_for_truth(scores_by_start, truth):
    """
    Align parsed scores to truth bins by start coordinate.
    Returns S (n_bins x 6); bins missing from the qcat (e.g. masked/blacklisted)
    are filled with zeros and flagged in a boolean 'present' array.
    """
    n = len(truth["starts"])
    S = np.zeros((n, N_TRACKS), dtype=np.float64)
    present = np.zeros(n, dtype=bool)
    for i, st in enumerate(truth["starts"]):
        v = scores_by_start.get(int(st))
        if v is not None:
            S[i] = v
            present[i] = True
    return S, present


# -----------------------------------------------------------------------------
# Recovery metrics.
# -----------------------------------------------------------------------------
def recovered_allocation(score_vec):
    """
    Normalize a per-track score vector to a recovered allocation r (sums to 1).
    If all scores are zero (fully masked / null), return None.
    """
    tot = score_vec.sum()
    if tot <= 0:
        return None
    return score_vec / tot


def enriched_mask(a_vec, q_vec):
    """Boolean mask of tracks the block pushes UP relative to background."""
    return a_vec > q_vec


def clr(vec):
    """Centered log-ratio of a strictly positive composition."""
    v = np.clip(vec, 1e-12, None)
    logv = np.log(v)
    return logv - logv.mean()


def aitchison_distance(p, q):
    """Aitchison distance between two compositions (Euclidean in CLR space)."""
    return float(np.linalg.norm(clr(p) - clr(q)))


def estimate_Q_from_truth(truth):
    """
    The background Q used by the scorer is the mean of P across all bins. We do
    not have direct access to the scorer's internal Q, but the NULL allocation
    in the truth table equals the planted Q exactly, so use that.
    """
    null_rows = truth["A"][truth["is_null"]]
    if len(null_rows) == 0:
        # fall back to global mean of planted P
        return truth["A"].mean(axis=0)
    return null_rows[0]


# -----------------------------------------------------------------------------
# Expected recovery under the clamped-KL metric.
#
# BEARING does NOT recover the planted allocation a directly. The per-track
# score is the clamped KL contribution
#       contrib_i = max(0, P_i * log2(P_i / Q_i)),   P = (1 - s) Q + s a
# so the recovered allocation r (= contrib normalized to sum 1) is a
# KL-reweighted, one-sided version of a, not a itself. Comparing r against
# raw a therefore puts points off the diagonal BY CONSTRUCTION (the reweighting
# over-weights tracks that are both high in P and far above Q, and zeros every
# track at or below Q). The correct ground-truth reference is the analytic
# expected allocation computed from the KNOWN (a, s, Q). If BEARING is faithful,
# recovered r lands on the diagonal against this expected allocation; residual
# scatter is then genuine NB-noise/smoothing recovery error, not expected
# reweighting.
# -----------------------------------------------------------------------------
def expected_clamped_kl_allocation(a, s, q):
    """
    Analytic allocation BEARING is expected to recover for a bin planted with
    allocation a at strength s against background q. Returns a 6-vector summing
    to 1 (or all zeros if the clamp kills everything, which should not happen
    for s > 0 with a != q).
    """
    P = (1.0 - s) * q + s * a
    P = P / P.sum()
    contrib = P * np.log2(P / (q + 1e-300) + 1e-300)
    contrib = np.clip(contrib, 0.0, None)
    tot = contrib.sum()
    if tot <= 0:
        return np.zeros_like(contrib)
    return contrib / tot


# -----------------------------------------------------------------------------
# Panel a: pooled recovered allocation per track, against BOTH references.
#   (left)  recovered r vs EXPECTED clamped-KL allocation  -- the correct
#           ground-truth test; should sit on the diagonal.
#   (right) recovered r vs RAW planted a -- shown for completeness and
#           explicitly labeled as including the expected KL reweighting, so the
#           systematic off-diagonal pull is understood rather than mistaken for
#           a recovery failure.
# Both use the survivor set (tracks with positive expected contribution), the
# same set BEARING can populate, so the comparison denominator is fixed.
# -----------------------------------------------------------------------------
def panel_a(truth, S, present, q, outpath):
    from scipy.stats import pearsonr, spearmanr

    # ---- collect per-bin points (survivor subset) AND per-block aggregates ----
    # per-bin arrays
    exp_vals, raw_vals, rec_vals = [], [], []
    track_idx = []
    # per-block accumulators: block_id -> list over bins of (expected, raw, rec)
    # restricted to that block's survivor set, plus the track indices.
    from collections import defaultdict
    block_exp = defaultdict(list)   # block -> list of e_s arrays
    block_raw = defaultdict(list)   # block -> list of a_s arrays
    block_rec = defaultdict(list)   # block -> list of r_s arrays
    block_surv = {}                 # block -> survivor track indices

    interior = (truth["bin_class"] == "interior") & present & (~truth["is_null"])
    for i in np.where(interior)[0]:
        r = recovered_allocation(S[i])
        if r is None:
            continue
        a = truth["A"][i]
        s = truth["strength"][i]
        bid = str(truth["block_id"][i])
        expected = expected_clamped_kl_allocation(a, s, q)
        surv = expected > 0
        if surv.sum() == 0:
            continue
        e_s = expected[surv] / expected[surv].sum()
        r_s = r[surv] / r[surv].sum() if r[surv].sum() > 0 else r[surv]
        a_s = a[surv] / a[surv].sum()
        surv_idx = np.where(surv)[0]

        for k, ti in enumerate(surv_idx):
            exp_vals.append(e_s[k]); raw_vals.append(a_s[k])
            rec_vals.append(r_s[k]); track_idx.append(ti)

        block_exp[bid].append(e_s)
        block_raw[bid].append(a_s)
        block_rec[bid].append(r_s)
        block_surv[bid] = surv_idx

    exp_vals = np.array(exp_vals); raw_vals = np.array(raw_vals)
    rec_vals = np.array(rec_vals); track_idx = np.array(track_idx)

    # per-block means and SDs (one row of points per block, per surviving track)
    blk_exp_pts, blk_rec_pts, blk_rec_sd, blk_track = [], [], [], []
    blk_raw_pts = []
    for bid in block_exp:
        E = np.vstack(block_exp[bid])   # (n_bins, n_surv)
        R = np.vstack(block_rec[bid])
        A = np.vstack(block_raw[bid])
        e_mean = E.mean(axis=0)         # expected is identical across bins; mean == value
        r_mean = R.mean(axis=0)
        r_sd = R.std(axis=0)
        a_mean = A.mean(axis=0)
        for k, ti in enumerate(block_surv[bid]):
            blk_exp_pts.append(e_mean[k]); blk_rec_pts.append(r_mean[k])
            blk_rec_sd.append(r_sd[k]); blk_track.append(ti)
            blk_raw_pts.append(a_mean[k])
    blk_exp_pts = np.array(blk_exp_pts); blk_rec_pts = np.array(blk_rec_pts)
    blk_rec_sd = np.array(blk_rec_sd); blk_track = np.array(blk_track)
    blk_raw_pts = np.array(blk_raw_pts)

    def stats(x, y):
        if len(x) <= 2:
            return float("nan"), float("nan"), float("nan"), float("nan")
        pe = pearsonr(x, y)[0]
        sp = spearmanr(x, y)[0]
        sl, ic = np.polyfit(x, y, 1)
        return pe, sp, sl, ic

    # headline stats are on the per-block means vs expected
    pe_be, sp_be, sl_be, ic_be = stats(blk_exp_pts, blk_rec_pts)
    pe_br, sp_br, sl_br, ic_br = stats(blk_raw_pts, blk_rec_pts)
    # per-bin stats (the noisy view)
    pe_e, sp_e, sl_e, ic_e = stats(exp_vals, rec_vals)
    pe_r, sp_r, sl_r, ic_r = stats(raw_vals, rec_vals)

    colors = plt.cm.tab10(np.linspace(0, 1, N_TRACKS))
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))

    def scatter_panel(ax, x, y, ti_arr, sd=None, xlabel="", title="",
                      pe=None, sp=None, sl=None, ic=None):
        for ti in range(N_TRACKS):
            m = ti_arr == ti
            if m.any():
                if sd is not None:
                    ax.errorbar(x[m], y[m], yerr=sd[m], fmt="o", ms=5,
                                color=colors[ti], ecolor=colors[ti],
                                elinewidth=0.8, capsize=2, alpha=0.85,
                                label=TRACK_NAMES[ti])
                else:
                    ax.scatter(x[m], y[m], s=10, alpha=0.4, color=colors[ti],
                               label=TRACK_NAMES[ti])
        ax.plot([0, 1], [0, 1], "k--", lw=1, label="identity")
        if sl is not None and np.isfinite(sl):
            xx = np.linspace(0, 1, 50)
            ax.plot(xx, sl * xx + ic, "r-", lw=1.2, label="fit slope=%.3f" % sl)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("recovered proportion")
        ax.set_title("%s\nPearson=%.3f  Spearman=%.3f  slope=%.3f"
                     % (title, pe, sp, sl))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.legend(fontsize=7, loc="upper left")

    # --- TOP ROW: per-block means (headline ground-truth statement) ---
    scatter_panel(axes[0, 0], blk_exp_pts, blk_rec_pts, blk_track, sd=blk_rec_sd,
                  xlabel="expected clamped-KL proportion (from known a, s, Q)",
                  title="Per-block mean recovered vs EXPECTED\n"
                        "(headline ground-truth test; error bars = bin-to-bin SD)",
                  pe=pe_be, sp=sp_be, sl=sl_be, ic=ic_be)
    scatter_panel(axes[0, 1], blk_raw_pts, blk_rec_pts, blk_track, sd=blk_rec_sd,
                  xlabel="planted proportion a (raw)",
                  title="Per-block mean recovered vs RAW planted a\n"
                        "(off-diagonal = expected KL reweighting, not error)",
                  pe=pe_br, sp=sp_br, sl=sl_br, ic=ic_br)

    # --- BOTTOM ROW: per-bin scatter (bin-level noise, honest) ---
    scatter_panel(axes[1, 0], exp_vals, rec_vals, track_idx, sd=None,
                  xlabel="expected clamped-KL proportion (from known a, s, Q)",
                  title="Per-bin recovered vs EXPECTED\n"
                        "(bin-level precision; scatter = NB noise + smoothing)",
                  pe=pe_e, sp=sp_e, sl=sl_e, ic=ic_e)
    scatter_panel(axes[1, 1], raw_vals, rec_vals, track_idx, sd=None,
                  xlabel="planted proportion a (raw)",
                  title="Per-bin recovered vs RAW planted a",
                  pe=pe_r, sp=sp_r, sl=sl_r, ic=ic_r)

    fig.suptitle("Recovery of planted composition "
                 "(top: per-block means; bottom: per-bin)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    return {
        "per_block_vs_expected": {"pearson": pe_be, "spearman": sp_be,
                                  "slope": sl_be, "intercept": ic_be,
                                  "n_points": int(len(blk_exp_pts))},
        "per_block_vs_raw_a": {"pearson": pe_br, "spearman": sp_br,
                               "slope": sl_br, "intercept": ic_br,
                               "n_points": int(len(blk_raw_pts))},
        "per_bin_vs_expected": {"pearson": pe_e, "spearman": sp_e,
                                "slope": sl_e, "intercept": ic_e,
                                "n_points": int(len(exp_vals))},
        "per_bin_vs_raw_a": {"pearson": pe_r, "spearman": sp_r,
                             "slope": sl_r, "intercept": ic_r,
                             "n_points": int(len(raw_vals))},
    }


# -----------------------------------------------------------------------------
# Panel b: named blocks planted vs recovered bars.
# -----------------------------------------------------------------------------
def block_mean_recovered(truth, S, present, block_substr, interior_only=True):
    """Mean recovered allocation across the interior bins of a named block."""
    sel = np.array([block_substr in str(b) for b in truth["block_id"]])
    if interior_only:
        sel = sel & (truth["bin_class"] == "interior")
    sel = sel & present & (~truth["is_null"])
    rs = []
    for i in np.where(sel)[0]:
        r = recovered_allocation(S[i])
        if r is not None:
            rs.append(r)
    if not rs:
        return None, 0
    return np.mean(rs, axis=0), len(rs)


def panel_b(truth, S, present, outpath):
    named = {"architecture": "named_architecture_s0.80",
             "transcription": "named_transcription_s0.80"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(N_TRACKS)
    width = 0.38
    results = {}
    for ax, (label, block_id) in zip(axes, named.items()):
        # planted allocation a for this block (from any of its bins)
        idx = np.where(truth["block_id"] == block_id)[0]
        a = truth["A"][idx[0]] if len(idx) else np.full(N_TRACKS, np.nan)
        r_mean, n = block_mean_recovered(truth, S, present, block_id)
        if r_mean is None:
            r_mean = np.full(N_TRACKS, np.nan)
        ax.bar(x - width / 2, a, width, label="planted a", color="#888888")
        ax.bar(x + width / 2, r_mean, width, label="recovered r", color="#1f77b4")
        ax.set_xticks(x)
        ax.set_xticklabels(TRACK_NAMES, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("proportion")
        ax.set_title("%s block (n=%d interior bins)" % (label, n))
        ax.legend(fontsize=8)
        results[label] = {"planted": a.tolist(), "recovered": r_mean.tolist(),
                          "n_bins": int(n)}
    fig.suptitle("Planted vs recovered allocation, named dual-scale blocks")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return results


# -----------------------------------------------------------------------------
# Panel c (rank recovery) + d (null behavior).
# -----------------------------------------------------------------------------
def rank_recovery(truth, S, present, q):
    """
    Fraction of planted interior bins whose top recovered track equals the
    planted top ENRICHED track (argmax of a over the enriched subset).
    """
    interior = (truth["bin_class"] == "interior") & present & (~truth["is_null"])
    hits = 0
    total = 0
    for i in np.where(interior)[0]:
        r = recovered_allocation(S[i])
        if r is None:
            continue
        a = truth["A"][i]
        planted_top = int(np.argmax(a))
        # only meaningful if planted_top is actually enriched
        if a[planted_top] <= q[planted_top]:
            continue
        recovered_top = int(np.argmax(r))
        total += 1
        if recovered_top == planted_top:
            hits += 1
    return (hits / total if total else float("nan")), hits, total


def panel_d(truth, S, present, outpath):
    """
    Null vs planted behavior.

    Two things are checked:
      (1) The total-score distribution: planted bins should score higher than
          null bins. We report the medians and the separation.
      (2) Spurious-call check. Because of realistic NB noise and smoothing, NULL
          bins do not score exactly zero -- their observed composition deviates
          slightly from Q, producing small positive clamped-KL scores. A raw
          "score > 0" test would therefore flag almost every bin and is not the
          right question. The biologically meaningful test is whether null bins
          produce a CONFIDENT dominant call, i.e. a total score above the
          threshold that the permutation null would deem significant. We
          approximate that threshold empirically as a high percentile of the
          null-bin score distribution (the same logic BEARING's empirical
          circular-shift null uses) and report the fraction of null vs planted
          bins exceeding it. A well-behaved method keeps the null exceedance at
          the nominal false-positive rate while the planted exceedance is high.
    """
    total_score = S.sum(axis=1)
    null = truth["is_null"] & present
    planted = (~truth["is_null"]) & present & (truth["bin_class"] == "interior")

    null_scores = total_score[null]
    planted_scores = total_score[planted]

    # Null-derived significance threshold: 95th percentile of null scores.
    # By construction ~5% of null bins exceed it (the nominal FPR); the test is
    # whether planted bins exceed it far more often, and whether the null
    # exceedance stays at ~5% rather than blowing up.
    alpha = 0.05
    if len(null_scores):
        thr = float(np.quantile(null_scores, 1.0 - alpha))
    else:
        thr = float("nan")
    null_exceed = float(np.mean(null_scores > thr)) if len(null_scores) else float("nan")
    planted_exceed = float(np.mean(planted_scores > thr)) if len(planted_scores) else float("nan")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    hi = float(np.percentile(
        np.concatenate([null_scores, planted_scores])
        if (len(null_scores) + len(planted_scores)) else [0, 1], 99))
    bins = np.linspace(0, max(1e-6, hi), 60)
    ax.hist(null_scores, bins=bins, alpha=0.6, label="null bins", color="#aaaaaa",
            density=True)
    ax.hist(planted_scores, bins=bins, alpha=0.6, label="planted (interior)",
            color="#d62728", density=True)
    if np.isfinite(thr):
        ax.axvline(thr, color="k", ls="--", lw=1,
                   label="null p%d threshold" % int((1 - alpha) * 100))
    ax.set_xlabel("total BEARING score per bin")
    ax.set_ylabel("density")
    ax.set_title("Score distribution: null vs planted")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.bar([0, 1], [null_exceed, planted_exceed],
           color=["#aaaaaa", "#d62728"])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["null bins", "planted bins"])
    ax.set_ylabel("fraction above null p95 threshold")
    ax.set_ylim(0, 1.05)
    ax.axhline(alpha, color="k", ls=":", lw=1, label="nominal FPR (%.2f)" % alpha)
    ax.set_title("Confident-call check\n(null near nominal FPR, planted high)")
    for xi, v in enumerate([null_exceed, planted_exceed]):
        ax.text(xi, v + 0.02, "%.3f" % v, ha="center", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    return {
        "null_median_score": float(np.median(null_scores)) if len(null_scores) else float("nan"),
        "planted_median_score": float(np.median(planted_scores)) if len(planted_scores) else float("nan"),
        "null_p95_threshold": thr,
        "null_fraction_above_threshold": null_exceed,
        "planted_fraction_above_threshold": planted_exceed,
    }


# -----------------------------------------------------------------------------
# Panel e: parameter sensitivity (the key panel).
# Sweeps (1) min-signal mask, (2) zero-clamp epsilon = PSEUDOCOUNT,
# (3) normalization quantile, and plots the architecture block's recovered
# RAD21+CTCF fraction vs parameter value. Flat lines => not a parameter artifact.
# -----------------------------------------------------------------------------
def rad21_ctcf_fraction(truth, S, present, block_substr):
    """Recovered RAD21+CTCF fraction averaged over a block's interior bins."""
    r_mean, n = block_mean_recovered(truth, S, present, block_substr)
    if r_mean is None:
        return float("nan")
    # RAD21 = index 4, CTCF = index 3
    return float(r_mean[3] + r_mean[4])


def insample_rescore(mod, raw_signal, q, pseudocount, min_signal):
    """
    Re-run the REAL scoring functions in-process with a chosen pseudocount
    (zero-clamp epsilon) and min_signal. raw_signal is (n_bins x 6) observed
    signal (clipped >=0). Returns the per-track score matrix (n_bins x 6).

    We temporarily set mod.PSEUDOCOUNT so signals_to_prob uses the swept value,
    then restore it. This uses the scorer's own math, not a reimplementation.
    """
    saved = mod.PSEUDOCOUNT
    try:
        mod.PSEUDOCOUNT = pseudocount
        P = mod.signals_to_prob(raw_signal)
        scores, _ = mod.kl_scores_per_bin(
            P, q, raw_signal_matrix=raw_signal, min_signal=min_signal,
            normalize_score=False)
    finally:
        mod.PSEUDOCOUNT = saved
    return scores


def read_raw_from_bigwigs(bw_paths, chrom, chrom_len, bin_size):
    """Read mean signal per bin from each BigWig -> (n_bins x 6) raw matrix."""
    import pyBigWig
    n_bins = chrom_len // bin_size
    raw = np.zeros((n_bins, len(bw_paths)), dtype=np.float64)
    for ti, p in enumerate(bw_paths):
        bw = pyBigWig.open(p)
        vals = bw.stats(chrom, 0, n_bins * bin_size, type="mean", nBins=n_bins)
        bw.close()
        raw[:, ti] = [v if v is not None else 0.0 for v in vals]
    raw[:, 2] = np.abs(raw[:, 2])  # RNA- negative-strand: abs() like the scorer
    return np.clip(raw, 0.0, None)


def panel_e(mod, truth, bw_paths, chrom_sizes, scorer_path, q, consts, outpath,
            min_signal_grid, eps_grid, quantile_methods, tmpdir,
            pseudocount=None):
    chrom = truth["chrom"]
    chrom_len = int(truth["ends"].max())
    bin_size = consts["BIN_SIZE"]

    # Two target blocks:
    #  - SATURATED: the named architecture block. RAD21/CTCF are the only tracks
    #    above background, so the clamp pins RAD21+CTCF fraction at 1.0 for any
    #    parameter setting. Flat-at-1.0 is real but uninformative (a ceiling).
    #  - INFORMATIVE: arch_plus_active, where ATAC and H3K27ac also survive, so
    #    RAD21+CTCF sits near ~0.73. A flat line HERE means the parameters do not
    #    perturb a NON-trivial composition -- the claim relevant to the real
    #    55-71% cohesin/CTCF split. We sweep the s=0.60 instance.
    blocks = [("saturated (named architecture)", "named_architecture_s0.80"),
              ("informative (arch+active, s0.60)", "arch_plus_active_s0.60")]

    def fracs_for(S, present):
        return [rad21_ctcf_fraction(truth, S, present, bid)
                for _, bid in blocks]

    # --- (1) min-signal sweep via real CLI runs ---
    ms_x, ms_y = [], [[] for _ in blocks]
    for ms in min_signal_grid:
        out = os.path.join(tmpdir, "sweep_ms_%g.qcat.bgz" % ms)
        run_scorer(scorer_path, bw_paths, chrom_sizes, out, min_signal=ms,
                   pseudocount=pseudocount)
        sb = parse_qcat(out, chrom)
        S, present = scores_matrix_for_truth(sb, truth)
        ms_x.append(ms)
        for bi, f in enumerate(fracs_for(S, present)):
            ms_y[bi].append(f)

    # --- (3) normalization quantile sweep via real CLI runs ---
    qm_labels, qm_y = [], [[] for _ in blocks]
    for method in quantile_methods:
        out = os.path.join(tmpdir, "sweep_norm_%s.qcat.bgz" % method)
        run_scorer(scorer_path, bw_paths, chrom_sizes, out,
                   normalize_method=method, normalize_tracks=True,
                   pseudocount=pseudocount)
        sb = parse_qcat(out, chrom)
        S, present = scores_matrix_for_truth(sb, truth)
        qm_labels.append(method)
        for bi, f in enumerate(fracs_for(S, present)):
            qm_y[bi].append(f)

    # --- (2) zero-clamp epsilon (PSEUDOCOUNT) sweep via in-process real funcs ---
    raw = read_raw_from_bigwigs(bw_paths, chrom, chrom_len, bin_size)
    P_default = mod.signals_to_prob(raw)
    q_runtime = P_default.mean(axis=0)
    eps_x, eps_y = [], [[] for _ in blocks]
    for eps in eps_grid:
        S = insample_rescore(mod, raw, q_runtime, pseudocount=eps,
                             min_signal=consts["MIN_SIGNAL"])
        present = np.ones(S.shape[0], dtype=bool)
        eps_x.append(eps)
        for bi, f in enumerate(fracs_for(S, present)):
            eps_y[bi].append(f)

    # --- plot: three panels, two lines each (saturated vs informative block) ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    block_colors = ["#9467bd", "#1f77b4"]   # saturated, informative
    block_styles = ["o--", "o-"]

    ax = axes[0]
    for bi, (label, _) in enumerate(blocks):
        ax.plot(ms_x, ms_y[bi], block_styles[bi], color=block_colors[bi],
                label=label)
    ax.set_xscale("log")
    ax.set_xlabel("min-signal mask threshold")
    ax.set_ylabel("recovered RAD21+CTCF fraction")
    ax.set_title("(1) low-signal mask")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7, loc="lower left")

    ax = axes[1]
    for bi in range(len(blocks)):
        ax.plot(eps_x, eps_y[bi], block_styles[bi], color=block_colors[bi])
    ax.set_xscale("log")
    ax.set_xlabel("zero-clamp epsilon (PSEUDOCOUNT)")
    ax.set_title("(2) clamp epsilon")
    ax.set_ylim(0, 1.05)

    ax = axes[2]
    for bi in range(len(blocks)):
        ax.plot(range(len(qm_labels)), qm_y[bi], block_styles[bi],
                color=block_colors[bi])
    ax.set_xticks(range(len(qm_labels)))
    ax.set_xticklabels(qm_labels, rotation=20, ha="right", fontsize=8)
    ax.set_xlabel("normalization method")
    ax.set_title("(3) normalization quantile")
    ax.set_ylim(0, 1.05)

    fig.suptitle("Parameter sensitivity of recovered composition "
                 "(flat = not a parameter artifact; informative block is the "
                 "non-saturated test)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    def pack(xs, ys):
        return {"x": list(map(float, xs)),
                "saturated": list(map(float, ys[0])),
                "informative": list(map(float, ys[1]))}
    return {
        "blocks": [b[0] for b in blocks],
        "min_signal": pack(ms_x, ms_y),
        "epsilon": pack(eps_x, eps_y),
        "normalization": {"methods": qm_labels,
                          "saturated": list(map(float, qm_y[0])),
                          "informative": list(map(float, qm_y[1]))},
    }


# -----------------------------------------------------------------------------
# Per-block planted-vs-recovered TSV with errors.
# -----------------------------------------------------------------------------
def write_per_block_tsv(truth, S, present, q, outpath):
    blocks = [b for b in sorted(set(truth["block_id"])) if b != "null"]
    header = (["block_id", "strength_s", "n_interior_bins"]
              + ["planted_%s" % t for t in TRACK_NAMES]
              + ["expected_%s" % t for t in TRACK_NAMES]
              + ["recovered_%s" % t for t in TRACK_NAMES]
              + ["abs_err_vs_expected_%s" % t for t in TRACK_NAMES]
              + ["L1_vs_expected", "aitchison_vs_expected",
                 "L1_vs_raw_a", "top_track_match"])
    with open(outpath, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for b in blocks:
            sel = (truth["block_id"] == b) & (truth["bin_class"] == "interior") \
                & present & (~truth["is_null"])
            idx = np.where(sel)[0]
            if len(idx) == 0:
                continue
            a = truth["A"][idx[0]]
            s = truth["strength"][idx[0]]
            expected = expected_clamped_kl_allocation(a, s, q)
            rs = [recovered_allocation(S[i]) for i in idx]
            rs = [r for r in rs if r is not None]
            if not rs:
                continue
            r_mean = np.mean(rs, axis=0)
            # error vs EXPECTED (the meaningful recovery error)
            abs_err_exp = np.abs(expected - r_mean)
            surv = expected > 0
            if surv.sum() > 0:
                e_s = expected[surv] / expected[surv].sum()
                r_s = r_mean[surv] / r_mean[surv].sum() if r_mean[surv].sum() > 0 else r_mean[surv]
                l1_exp = float(np.abs(e_s - r_s).sum())
                ait_exp = (aitchison_distance(e_s, r_s)
                           if r_mean[surv].sum() > 0 else float("nan"))
                a_s = a[surv] / a[surv].sum()
                l1_raw = float(np.abs(a_s - r_s).sum())
            else:
                l1_exp, ait_exp, l1_raw = float("nan"), float("nan"), float("nan")
            top_match = int(np.argmax(expected) == np.argmax(r_mean))
            row = [b, "%.3f" % s, str(len(rs))]
            row += ["%.4f" % v for v in a]
            row += ["%.4f" % v for v in expected]
            row += ["%.4f" % v for v in r_mean]
            row += ["%.4f" % v for v in abs_err_exp]
            row += ["%.4f" % l1_exp, "%.4f" % ait_exp, "%.4f" % l1_raw,
                    str(top_match)]
            fh.write("\t".join(row) + "\n")
    return outpath


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Evaluate BEARING ground-truth recovery on synthetic data "
                    "using the real bigwig_to_qcat.py pipeline.")
    ap.add_argument("--truth", required=True, help="<prefix>_truth.tsv")
    ap.add_argument("--bw", nargs=6, required=True, metavar="BW",
                    help="six BigWigs in order ATAC RNA+ RNA- CTCF RAD21 H3K27ac")
    ap.add_argument("--chrom-sizes", required=True, help="<prefix>.chrom.sizes")
    ap.add_argument("--scorer", default="./bigwig_to_qcat.py",
                    help="path to bigwig_to_qcat.py (default: ./bigwig_to_qcat.py)")
    ap.add_argument("--outdir", default="bearing_recovery_eval",
                    help="output directory (default: bearing_recovery_eval)")
    ap.add_argument("--min-signal-grid", default="0.001,0.01,0.05,0.1,0.5",
                    help="comma-separated min-signal thresholds for panel e")
    ap.add_argument("--epsilon-grid", default="1e-9,1e-7,1e-6,1e-5,1e-3",
                    help="comma-separated zero-clamp epsilons (PSEUDOCOUNT)")
    ap.add_argument("--quantile-methods", default="nonzero-quantile,quantile",
                    help="comma-separated normalize-method values for panel e")
    ap.add_argument("--pseudocount", type=float, default=None,
                    help="zero-clamp epsilon passed to EVERY scorer run "
                         "(baseline panels a-d and the panel-e mask/normalization "
                         "runs). When set, all recovery metrics are computed at "
                         "this pseudocount, so a sweep over it (e.g. from "
                         "run_bearing_recovery_sweep.sh) genuinely recomputes "
                         "slope/Pearson/rank, not just panel e. The in-process "
                         "epsilon sweep (--epsilon-grid) is independent of this.")
    ap.add_argument("--normalize-tracks", action="store_true",
                    help="pass --normalize-tracks to EVERY scorer run, so the "
                         "headline recovery metrics (baseline panels a-d) are "
                         "computed with per-track (nonzero-quantile) "
                         "normalization applied before P is formed. When set, "
                         "the expected clamped-KL allocation is unchanged (it is "
                         "analytic from the planted a, s, Q), so recovered-vs-"
                         "expected remains a fair test of whether normalization "
                         "preserves faithful attribution. Lets a sweep treat "
                         "track normalization as an on/off knob.")
    ap.add_argument("--keep-tmp", action="store_true",
                    help="keep intermediate qcat files in the outdir")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = os.path.join(args.outdir, "tmp")
    os.makedirs(tmpdir, exist_ok=True)

    mod = load_scorer(args.scorer)
    consts = preflight(mod, args.scorer)

    truth = load_truth(args.truth)
    q = estimate_Q_from_truth(truth)
    print("background Q (from null bins): %s"
          % ", ".join("%s=%.3f" % (TRACK_NAMES[i], q[i]) for i in range(N_TRACKS)))

    if args.pseudocount is not None:
        print("pseudocount (zero-clamp epsilon) for all scorer runs: %g"
              % args.pseudocount)
    if args.normalize_tracks:
        print("per-track normalization (--normalize-tracks) applied to all "
              "scorer runs")

    # Baseline run (panels a-d). Uses --pseudocount if given, else scorer default.
    base_out = os.path.join(tmpdir, "baseline.qcat.bgz")
    run_scorer(args.scorer, args.bw, args.chrom_sizes, base_out,
               normalize_tracks=args.normalize_tracks,
               pseudocount=args.pseudocount)
    scores_by_start = parse_qcat(base_out, truth["chrom"])
    S, present = scores_matrix_for_truth(scores_by_start, truth)
    print("parsed %d scored bins; %d/%d truth bins present in qcat"
          % (len(scores_by_start), int(present.sum()), len(truth["starts"])))

    summary = {}
    summary["panel_a"] = panel_a(truth, S, present, q,
                                 os.path.join(args.outdir, "fig_a_recovered_vs_planted.png"))
    summary["panel_b"] = panel_b(truth, S, present,
                                 os.path.join(args.outdir, "fig_b_named_blocks.png"))
    rr, hits, tot = rank_recovery(truth, S, present, q)
    summary["panel_c_rank_recovery"] = {"fraction": rr, "hits": hits, "total": tot}
    summary["panel_d"] = panel_d(truth, S, present,
                                 os.path.join(args.outdir, "fig_d_null_vs_planted.png"))
    summary["panel_e"] = panel_e(
        mod, truth, args.bw, args.chrom_sizes, args.scorer, q, consts,
        os.path.join(args.outdir, "fig_e_parameter_sensitivity.png"),
        min_signal_grid=[float(x) for x in args.min_signal_grid.split(",")],
        eps_grid=[float(x) for x in args.epsilon_grid.split(",")],
        quantile_methods=[m.strip() for m in args.quantile_methods.split(",")],
        tmpdir=tmpdir, pseudocount=args.pseudocount)

    tsv = write_per_block_tsv(truth, S, present, q,
                              os.path.join(args.outdir, "recovery_per_block.tsv"))

    # summary text
    sumpath = os.path.join(args.outdir, "recovery_summary.txt")
    with open(sumpath, "w") as fh:
        fh.write("BEARING ground-truth recovery summary\n")
        fh.write("=====================================\n\n")
        a = summary["panel_a"]
        pbe = a["per_block_vs_expected"]; pbr = a["per_block_vs_raw_a"]
        pne = a["per_bin_vs_expected"]; pnr = a["per_bin_vs_raw_a"]
        fh.write("Panel a (recovery of planted composition):\n")
        fh.write("  PER-BLOCK mean vs EXPECTED (headline): "
                 "Pearson=%.3f Spearman=%.3f slope=%.3f n=%d\n"
                 % (pbe["pearson"], pbe["spearman"], pbe["slope"], pbe["n_points"]))
        fh.write("  PER-BLOCK mean vs RAW a (incl KL reweighting): "
                 "Pearson=%.3f Spearman=%.3f slope=%.3f n=%d\n"
                 % (pbr["pearson"], pbr["spearman"], pbr["slope"], pbr["n_points"]))
        fh.write("  per-bin vs EXPECTED (bin-level precision): "
                 "Pearson=%.3f Spearman=%.3f slope=%.3f n=%d\n"
                 % (pne["pearson"], pne["spearman"], pne["slope"], pne["n_points"]))
        fh.write("  per-bin vs RAW a: "
                 "Pearson=%.3f Spearman=%.3f slope=%.3f n=%d\n"
                 % (pnr["pearson"], pnr["spearman"], pnr["slope"], pnr["n_points"]))
        fh.write("\nPanel b (named blocks):\n")
        for label, d in summary["panel_b"].items():
            fh.write("  %s (n=%d): planted=%s recovered=%s\n"
                     % (label, d["n_bins"],
                        ["%.2f" % v for v in d["planted"]],
                        ["%.2f" % v for v in d["recovered"]]))
        c = summary["panel_c_rank_recovery"]
        fh.write("\nPanel c (rank recovery): %.3f (%d/%d planted bins)\n"
                 % (c["fraction"], c["hits"], c["total"]))
        d = summary["panel_d"]
        fh.write("\nPanel d (null behavior): null median score=%.4g, "
                 "planted median score=%.4g, null p95 threshold=%.4g, "
                 "null fraction above threshold=%.3f, "
                 "planted fraction above threshold=%.3f\n"
                 % (d["null_median_score"], d["planted_median_score"],
                    d["null_p95_threshold"],
                    d["null_fraction_above_threshold"],
                    d["planted_fraction_above_threshold"]))
        e = summary["panel_e"]
        fh.write("\nPanel e (parameter sensitivity, RAD21+CTCF recovered "
                 "fraction; reported for the SATURATED named block and the "
                 "INFORMATIVE non-saturated arch+active block):\n")
        for key, label in [("min_signal", "min-signal"),
                           ("epsilon", "epsilon"),
                           ("normalization", "normalize")]:
            xs = e[key].get("x", e[key].get("methods"))
            fh.write("  %s x=%s\n" % (label, xs))
            fh.write("    saturated   -> %s\n"
                     % ["%.3f" % v for v in e[key]["saturated"]])
            fh.write("    informative -> %s\n"
                     % ["%.3f" % v for v in e[key]["informative"]])
        # robustness: spread of the INFORMATIVE block across all sweeps
        infovals = (e["min_signal"]["informative"] + e["epsilon"]["informative"]
                    + e["normalization"]["informative"])
        infovals = [v for v in infovals if v == v]
        if infovals:
            fh.write("  INFORMATIVE-block RAD21+CTCF fraction range across all "
                     "sweeps: %.3f - %.3f (spread %.3f)\n"
                     % (min(infovals), max(infovals),
                        max(infovals) - min(infovals)))

    print("\nWrote outputs to %s/:" % args.outdir)
    for f in ("fig_a_recovered_vs_planted.png", "fig_b_named_blocks.png",
              "fig_d_null_vs_planted.png", "fig_e_parameter_sensitivity.png",
              "recovery_per_block.tsv", "recovery_summary.txt"):
        print("  %s" % f)

    if not args.keep_tmp:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
