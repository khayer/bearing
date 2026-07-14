"""Tier 2 statistics-core tests.

Unlike the import guard (structural), these exercise the actual scientific
machinery that the manuscript's significance claims rest on and that had no
unit coverage: empirical p-values, BH-FDR, JSD, the differential sign
convention, the circular-shift permutation null, and the regional-enrichment
batch path (manuscript Table S3).

All assertions are derived from the functions' documented behaviour, verified
against the source in bearing_pvalue.py / compare_qcat.py / regional_enrichment.py
/ shift_bigwig.py. ASCII only.
"""

import csv
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


# ===========================================================================
# bearing_pvalue.empirical_pvals -- empirical one-sided p, P(X >= s) under null
# ===========================================================================

def test_empirical_pvals_range_and_bounds():
    from bearing_pvalue import empirical_pvals

    null = np.sort(np.linspace(0.0, 1.0, 100))
    obs = np.array([-5.0, 0.0, 0.5, 1.0, 5.0])
    p = empirical_pvals(obs, null)
    n = len(null)
    # p is bounded in [1/(n+1), 1] by construction (Davison-Hinkley + clip).
    assert np.all(p >= 1.0 / (n + 1) - 1e-12)
    assert np.all(p <= 1.0 + 1e-12)
    # A score below the entire null gets p == 1.0 (all null >= s).
    assert p[0] == pytest.approx(1.0)
    # A score above the entire null hits the floor 1/(n+1).
    assert p[-1] == pytest.approx(1.0 / (n + 1))


def test_empirical_pvals_monotonic_nonincreasing():
    """Higher observed score => smaller-or-equal p-value. This is the core
    property a reviewer will check: significance must be monotone in signal."""
    from bearing_pvalue import empirical_pvals

    null = np.sort(np.random.default_rng(0).normal(size=500))
    obs = np.linspace(-3, 3, 50)
    p = empirical_pvals(obs, null)
    # non-increasing along increasing observed score
    assert np.all(np.diff(p) <= 1e-12)


# ===========================================================================
# bearing_pvalue.bh_fdr -- Benjamini-Hochberg
# ===========================================================================

def test_bh_fdr_known_example():
    from bearing_pvalue import bh_fdr

    # Hand-computable: p*n/rank then right-to-left cumulative min.
    pvals = np.array([0.001, 0.5, 0.9])
    rejected, adj = bh_fdr(pvals, alpha=0.05)
    # rank1 0.001*3/1=0.003 ; rank2 0.5*3/2=0.75 ; rank3 0.9*3/3=0.9
    assert adj[0] == pytest.approx(0.003, rel=1e-9)
    assert adj[1] == pytest.approx(0.75, rel=1e-9)
    assert adj[2] == pytest.approx(0.9, rel=1e-9)
    assert list(rejected) == [True, False, False]


def test_bh_fdr_all_reject_boundary():
    from bearing_pvalue import bh_fdr

    # p_i = i*alpha/n makes every adjusted p == alpha -> all rejected.
    n = 5
    alpha = 0.05
    pvals = np.array([(i + 1) * alpha / n for i in range(n)])
    rejected, adj = bh_fdr(pvals, alpha)
    assert np.all(adj <= alpha + 1e-12)
    assert np.all(rejected)


def test_bh_fdr_adjusted_in_unit_interval_and_capped():
    from bearing_pvalue import bh_fdr

    rng = np.random.default_rng(1)
    pvals = rng.uniform(size=200)
    _, adj = bh_fdr(pvals, alpha=0.05)
    assert np.all(adj >= 0.0)
    assert np.all(adj <= 1.0)


# ===========================================================================
# compare_qcat -- JSD symmetry/bounds + probability normalization
# ===========================================================================

def test_scores_to_prob_rows_sum_to_one():
    from compare_qcat import scores_to_prob

    m = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [5.0, 0.0, 0.0]])
    P = scores_to_prob(m)
    assert np.allclose(P.sum(axis=1), 1.0)
    assert np.all(P >= 0.0)


def test_normalize_dist_zero_vector_is_uniform():
    from compare_qcat import _normalize_dist

    out = _normalize_dist(np.zeros(4))
    assert np.allclose(out, 0.25)


def test_js_divergence_symmetry_and_zero():
    from compare_qcat import js_divergence

    rng = np.random.default_rng(2)
    P = rng.dirichlet(np.ones(6), size=40)
    Q = rng.dirichlet(np.ones(6), size=40)
    assert js_divergence(P, Q) == pytest.approx(js_divergence(Q, P))
    # identical distributions -> zero divergence
    assert js_divergence(P, P) == pytest.approx(0.0, abs=1e-12)


def test_js_divergence_bounds_and_maximum():
    from compare_qcat import js_divergence

    # Disjoint one-hot distributions -> maximal JSD == 1.0 (base-2).
    P = np.array([[1.0, 0.0]])
    Q = np.array([[0.0, 1.0]])
    jsd = js_divergence(P, Q)
    assert jsd == pytest.approx(1.0, abs=1e-9)
    assert 0.0 <= jsd <= 1.0


# ===========================================================================
# compare_qcat.diff_qcat -- differential sign convention (diff = A - B)
# ===========================================================================

def test_diff_qcat_sign_convention(tmp_path):
    """diff_i = KL_A_i - KL_B_i ; positive => state more active in A.
    Verified end to end by reading the written qcat back."""
    pytest.importorskip("pysam")
    import gzip

    from compare_qcat import diff_qcat

    bins = [("chr1", 0, 200), ("chr1", 200, 400)]
    # 3 states. bin0: A dominant in state0 ; bin1: B dominant in state0.
    A = np.array([[3.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    B = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

    out = tmp_path / "diff.qcat.bgz"
    diff_qcat(bins, A, B, str(out))
    assert out.exists()

    rows = {}
    with gzip.open(str(out), "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            chrom, start = parts[0], int(parts[1])
            qjson = parts[3].split("qcat:", 1)[1]
            import json
            rows[(chrom, start)] = json.loads(qjson)

    # bin0: top pair is state 0 (1-indexed -> 1) with POSITIVE score (A > B).
    top0 = rows[("chr1", 0)][0]
    assert top0[1] == 1 and top0[0] == pytest.approx(3.0)
    # bin1: top pair is state 0 with NEGATIVE score (B > A).
    top1 = rows[("chr1", 200)][0]
    assert top1[1] == 1 and top1[0] == pytest.approx(-3.0)


# ===========================================================================
# shift_bigwig.circular_shift_bigwig -- permutation-null shift properties
# ===========================================================================

def _write_bigwig(path, chrom, bin_size, values):
    """Write one value per bin over a chromosome of len(values)*bin_size."""
    import pyBigWig

    chrom_len = len(values) * bin_size
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, chrom_len)])
    starts = [i * bin_size for i in range(len(values))]
    ends = [s + bin_size for s in starts]
    bw.addEntries([chrom] * len(values), starts, ends=ends,
                  values=[float(v) for v in values])
    bw.close()
    return chrom_len


def _read_bin_values(path, chrom, bin_size, n_bins):
    import pyBigWig

    bw = pyBigWig.open(str(path))
    out = []
    for i in range(n_bins):
        v = bw.stats(chrom, i * bin_size, (i + 1) * bin_size, type="mean")[0]
        out.append(0.0 if v is None else float(v))
    bw.close()
    return np.array(out)


def test_circular_shift_full_period_is_identity(tmp_path):
    pytest.importorskip("pyBigWig")
    from shift_bigwig import circular_shift_bigwig

    bin_size = 200
    vals = [0.0, 1.0, 2.0, 0.0, 5.0, 0.0, 3.0, 0.0, 0.0, 7.0]
    n = len(vals)
    src = tmp_path / "in.bw"
    dst = tmp_path / "out.bw"
    _write_bigwig(src, "chrSim", bin_size, vals)

    # A shift of exactly n_bins * bin_size wraps to a no-op.
    circular_shift_bigwig(str(src), str(dst), shift=n * bin_size, bin_size=bin_size)
    got = _read_bin_values(dst, "chrSim", bin_size, n)
    assert np.allclose(got, vals)


def test_circular_shift_preserves_total_signal_and_wraps(tmp_path):
    pytest.importorskip("pyBigWig")
    from shift_bigwig import circular_shift_bigwig

    bin_size = 200
    vals = [5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0]
    n = len(vals)
    src = tmp_path / "in.bw"
    dst = tmp_path / "out.bw"
    _write_bigwig(src, "chrSim", bin_size, vals)

    # shift = 2 bins (np.roll right by 2). bin0 -> bin2, bin9 -> bin1 (wrap).
    circular_shift_bigwig(str(src), str(dst), shift=2 * bin_size, bin_size=bin_size)
    got = _read_bin_values(dst, "chrSim", bin_size, n)

    # total signal mass is conserved (roll permutes, never creates/destroys)
    assert got.sum() == pytest.approx(sum(vals))
    # explicit wrap check
    assert got[2] == pytest.approx(5.0)   # bin0 -> bin2
    assert got[1] == pytest.approx(9.0)   # bin9 wraps to bin1
    assert got[0] == pytest.approx(0.0)


# ===========================================================================
# regional_enrichment batch -- Table S3 code path (folded-in coverage)
# ===========================================================================

def _run_regional_batch(tmp_path, region_assign):
    """Run `regional_enrichment.py batch` on a synthetic diff-stats table with
    one sub-bin (CBE-sized) region, under the given --region-assign mode.
    Returns the parsed output row for the region."""
    locus_chrom, locus_start = "chr6", 40790000
    bin_size = 200
    n_bins = 40

    # Build a diff-stats TSV: header chrom start end pval bearing_score.
    stats = tmp_path / ("diff_%s.stats.tsv" % region_assign)
    with open(stats, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "start", "end", "pval", "bearing_score"])
        for i in range(n_bins):
            start = locus_start + 1000 + i * bin_size
            end = start + bin_size
            if i == 5 or i < 4:
                w.writerow([locus_chrom, start, end, 0.001, 1.5])   # sig, A>B
            elif i < 8:
                w.writerow([locus_chrom, start, end, 0.001, -1.2])  # sig, B>A
            else:
                w.writerow([locus_chrom, start, end, 0.8, 0.1])     # not sig

    # An 18 bp CBE-sized region sitting INSIDE bin index 5 -- no bin start lands
    # in it, so 'start' assignment finds zero bins; 'overlap' finds bin 5.
    bin5_start = locus_start + 1000 + 5 * bin_size
    region_start = bin5_start + 50
    region_end = region_start + 18
    regions = tmp_path / "cbe.bed"
    with open(regions, "w") as fh:
        fh.write("%s\t%d\t%d\t%s\n" % (locus_chrom, region_start, region_end, "CBE1"))

    out = tmp_path / ("enrich_%s.tsv" % region_assign)
    locus = "%s:%d-%d" % (locus_chrom, locus_start, locus_start + n_bins * bin_size + 2000)

    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "regional_enrichment.py"), "batch",
         "--diff-table", str(stats), "--regions", str(regions),
         "--region-assign", region_assign,
         "--locus", locus, "--out", str(out)],
        cwd=str(REPO_ROOT), capture_output=True, text=True,
    )
    assert proc.returncode == 0, "batch failed:\n%s\n%s" % (proc.stdout, proc.stderr)
    assert out.exists()
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) >= 1, "expected at least one region row"
    return rows[0]


def test_regional_enrichment_batch_table_s3_overlap(tmp_path):
    """Table S3 code path: `regional_enrichment.py batch --region-assign overlap`
    over a sub-bin CBE region produces a non-empty, well-formed table with the
    overlapping bin assigned."""
    r = _run_regional_batch(tmp_path, "overlap")
    assert r["region_name"] == "CBE1"
    assert int(r["L_region_bins"]) > 0          # overlap captures bin 5
    assert int(r["k"]) > 0                       # and it is significant
    assert "q_combined" in r and r["q_combined"] != ""


def test_regional_enrichment_start_assign_misses_subbin_region(tmp_path):
    """Motivation for --region-assign overlap: under the default 'start' mode a
    sub-bin CBE captures zero bins (no 200 bp grid start lands inside 18 bp)."""
    r = _run_regional_batch(tmp_path, "start")
    assert int(r["L_region_bins"]) == 0
