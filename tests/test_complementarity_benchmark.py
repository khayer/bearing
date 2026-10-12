# ----------------------------------------------------------------------
# File     : tests/test_complementarity_benchmark.py
# Version  : 1.0.0
# Date     : 2026-10-11
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
