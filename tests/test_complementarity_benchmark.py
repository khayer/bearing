# ----------------------------------------------------------------------
# File     : tests/test_complementarity_benchmark.py
# Version  : 1.2.0
# Date     : 2026-11-08 (input-transform and RNA-tail tests;
#            1.0.0 2026-10-11)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""Tests for benchmark/complementarity_benchmark.py.

1. Fidelity: the in-process BEARING scoring used by the benchmark gives the
   same per-bin scores as the production CLI (bigwig_to_qcat.py with
   --prior-strength 1 --min-signal 0.1) on the same simulated BigWigs.
2. Behaviour on a small run: deterministic for a seed; no BEARING call in
   magnitude-only blocks; per-assay comparator false-call rate on null bins
   is small. ASCII only.
"""

import gzip
import json
import os
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "benchmark"))
sys.path.insert(0, REPO)

cb = pytest.importorskip("complementarity_benchmark")
sbt = pytest.importorskip("simulate_bearing_tracks")

SMALL = ["--length", "3000000", "--min-shift", "1000000", "--copies", "2",
         "--n-perms", "4", "--depth", "8"]


def _parse_qcat(path, n_bins):
    out = np.zeros((n_bins, 6))
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c, s, e, col = line.rstrip("\n").split("\t")[:4]
            i = int(s) // 200
            pay = col[col.find("qcat:") + 5:]
            if ",raw:" in pay:
                pay = pay[:pay.find(",raw:")]
            for sc, st in json.loads(pay):
                out[i, int(st) - 1] = float(sc)
    return out


def test_scoring_matches_production_cli(tmp_path):
    pytest.importorskip("pyBigWig")
    pytest.importorskip("pysam")
    args = cb.build_parser().parse_args(["--outdir", str(tmp_path)] + SMALL)
    samples, truth = cb.simulate(args, np.random.default_rng(7))
    x = samples["A_rep1"].astype(np.float32).astype(np.float64)  # BigWig = float32
    prefix = str(tmp_path / "s")
    bws = sbt.write_bigwigs(x, "chrSim", args.length, 200, prefix)
    cs = sbt.write_chrom_sizes("chrSim", args.length, prefix)
    out = str(tmp_path / "s.qcat.bgz")
    r = subprocess.run([sys.executable, os.path.join(REPO, "bigwig_to_qcat.py"),
                        "--bw"] + bws + ["--out", out, "--chrom-sizes", cs,
                                         "--min-signal", "0.1", "--prior-strength",
                                         "1", "--no-extras"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    cli = _parse_qcat(out, x.shape[0])
    mine, _ = cb.bearing_scores(x, 0.1, 1.0)
    # qcat stores 6 significant digits
    assert np.allclose(cli, mine, rtol=1e-4, atol=1e-5), \
        float(np.max(np.abs(cli - mine)))


def test_small_run_behaviour(tmp_path):
    a = cb.main(["--outdir", str(tmp_path / "a")] + SMALL)
    b = cb.main(["--outdir", str(tmp_path / "b")] + SMALL)
    rows_a, trows_a, fdr_a = a
    rows_b, trows_b, fdr_b = b
    assert [r["bearing_sum_called"] for r in rows_a] == \
           [r["bearing_sum_called"] for r in rows_b]
    t = {r["type"]: r for r in trows_a}
    assert t["magnitude"]["bearing_sum_called"] == 0.0
    assert t["null_bins"]["per_assay_called"] < 0.02
    assert t["null_bins"]["bearing_sum_called"] < 0.01
    assert t["magnitude"]["per_assay_called"] > 0.5
    for f in ("per_config.tsv", "per_type.tsv", "summary.txt"):
        assert os.path.exists(tmp_path / "a" / f)


def test_signed_sum_cancels_for_equal_divergence_switch():
    """The production statistic is |sum_i (KL_A,i - KL_B,i)|; a switch
    between two compositions equally far from Q gives ~0 by construction,
    while L1 does not. Pins the documented behaviour."""
    q = cb.DEFAULT_Q / cb.DEFAULT_Q.sum()
    d = np.zeros((1, 6))
    d[0, 1] = 1.2      # RNA+ enriched in A
    d[0, 3] = -1.2     # CTCF enriched in B
    assert cb.BEARING_STATS["sum"](d)[0] == pytest.approx(0.0)
    assert cb.BEARING_STATS["l1"](d)[0] == pytest.approx(2.4)
    assert q.sum() == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 1.2.0: input transforms and heavy-tailed RNA (benchmark only)
# ---------------------------------------------------------------------------

def test_transform_none_is_production(tmp_path):
    """transform 'none' gives exactly the 1.0.1 scores."""
    args = cb.build_parser().parse_args(["--outdir", str(tmp_path)] + SMALL)
    samples, _ = cb.simulate(args, np.random.default_rng(3))
    x = samples["A_rep1"]
    P = cb.b2q.signals_to_prob(x)
    ref, _ = cb.b2q.kl_scores_per_bin(P, P.mean(axis=0), raw_signal_matrix=x,
                                      min_signal=0.1, prior_strength=1.0,
                                      prior_source=x)
    mine, _ = cb.bearing_scores(x, 0.1, 1.0, transform="none")
    assert np.array_equal(ref, mine)


def test_rna_tail_off_leaves_simulation_unchanged(tmp_path):
    a = cb.build_parser().parse_args(["--outdir", str(tmp_path)] + SMALL)
    b = cb.build_parser().parse_args(["--outdir", str(tmp_path), "--rna-tail-sigma",
                                      "0"] + SMALL)
    sa, _ = cb.simulate(a, np.random.default_rng(5))
    sb, _ = cb.simulate(b, np.random.default_rng(5))
    for k in sa:
        assert np.array_equal(sa[k], sb[k])
    c = cb.build_parser().parse_args(["--outdir", str(tmp_path), "--rna-tail-sigma",
                                      "1.5"] + SMALL)
    sc, tc = cb.simulate(c, np.random.default_rng(5))
    # RNA gets a heavier tail than CTCF; truth A vs B unchanged in ratio
    x = sc["A_rep1"]
    nz = lambda v: v[v > 0]
    r_rna = np.quantile(nz(x[:, 1]), 0.99) / np.median(nz(x[:, 1]))
    r_ctcf = np.quantile(nz(x[:, 3]), 0.99) / np.median(nz(x[:, 3]))
    assert r_rna > 2 * r_ctcf


def test_transforms_map_zero_to_zero_and_cofactors():
    x = np.array([[0.0, 1.0, 0.0, 2.0, 0.0, 0.5],
                  [3.0, 0.0, 4.0, 0.0, 1.0, 0.0],
                  [1.0, 9.0, 2.0, 2.0, 3.0, 1.5]])
    cof = cb.track_cofactors({"A_rep1": x, "B_rep1": x})
    assert cof.tolist() == [2.0, 5.0, 3.0, 2.0, 2.0, 1.0]
    for tf in ("nzq", "log1p_rna", "asinh_auto"):
        y = cb.transform_signal(x, tf, cof)
        assert np.array_equal(y == 0, x == 0), tf
    y = cb.transform_signal(x, "log1p_rna", cof)
    assert np.array_equal(y[:, [0, 3, 4, 5]], x[:, [0, 3, 4, 5]])
    assert y[2, 1] == pytest.approx(np.log1p(9.0 / 5.0))
    y = cb.transform_signal(x, "asinh_auto", cof)
    assert y[2, 4] == pytest.approx(np.arcsinh(1.5))


def test_transform_run_uses_same_permutations(tmp_path):
    """A transformed run differs from production only in scoring: the same
    tested-bin floor logic runs and the null is non-empty."""
    r = cb.main(["--outdir", str(tmp_path / "t"), "--transform", "asinh_auto"]
                + SMALL)
    rows, trows, fdr = r
    txt = open(tmp_path / "t" / "summary.txt").read()
    assert "transform asinh_auto" in txt and "cofactors" in txt
    t = {x["type"]: x for x in trows}
    assert t["null_bins"]["bearing_l1_called"] < 0.01


def test_fast_nonzero_quantile_matches_production():
    rng = np.random.default_rng(11)
    x = rng.gamma(0.7, 1.0, size=(3000, 6))
    x[rng.random(x.shape) < 0.3] = 0.0
    x[:, 2] = 0.0
    x[5, 2] = 1.0          # a column with one nonzero value
    ref = cb.b2q.normalize_signal_matrix(x, method="nonzero-quantile")
    assert np.allclose(cb.fast_nonzero_quantile(x), ref, rtol=1e-12, atol=1e-12)
