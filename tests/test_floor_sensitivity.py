"""Tests for the floor-sensitivity supplement (SuppFig S11).

Covers the sweep's per-bin filename tag, the plotting script's pure math
(Panel A BES / detection curve, checked against the spec's uniform-Q sanity
table), the input parsers, and PDF reproducibility (byte-identical re-run).
Does not need production data. ASCII only.
"""

import gzip
import os

import numpy as np
import pytest

import pvminsig_sweep as sweep
import plot_floor_sensitivity as pfs


def test_floor_tag_formats():
    assert sweep.floor_tag(0.25) == "0p25"
    assert sweep.floor_tag(1.0) == "1p0"
    assert pfs.floor_tag(0.5) == "0p5"
    assert pfs.floor_tag(0.1) == "0p1"


def test_bes_clamped_pure_track_ceiling():
    q = 0.0836
    # pure track A: BES([1,0]) = -log2(q) = ceiling
    assert pfs.bes_scalar_q([1.0, 0.0], q) == pytest.approx(-np.log2(q), rel=1e-9)
    # a state at exactly q contributes 0 (clamped)
    assert pfs.bes_scalar_q([q, 1 - q], q) == pytest.approx((1 - q) * np.log2((1 - q) / q), rel=1e-9)


def test_panel_a_sanity_table_uniform_q():
    """Spec sanity table (uniform Q=0.0836, ceiling 3.58): floor 0.5 -> ~89/11,
    floor 1.0 -> not reached before 50/50."""
    q = 0.0836
    xs = np.arange(0.0, 0.5 + 1e-9, 0.001)
    ceiling, _ = pfs.detection_curve(q, xs)
    assert ceiling == pytest.approx(3.58, abs=0.02)
    x05, ratio05 = pfs.smallest_detectable_split(q, 0.5, xs)
    # within one percent of the tabulated 89/11
    assert x05 == pytest.approx(0.11, abs=0.01)
    # floor 1.0 is reached exactly at 50/50 ("not reached before 50/50"): y(0.5)=1.0
    x10, ratio10 = pfs.smallest_detectable_split(q, 1.0, xs)
    assert x10 == pytest.approx(0.5, abs=0.01) and ratio10 == "50/50"
    # floor 1.5 is never reached on 0..0.5
    assert pfs.smallest_detectable_split(q, 1.5, xs) is None


def test_read_sweep_and_stats(tmp_path):
    sw = tmp_path / "sweep.tsv"
    sw.write_text(
        "comparison\tfloor\tbins_tested\tnull_scores\tbins_bh_significant\tmin_pval\n"
        "DN_vs_DP\t0.5\t4200000\t1100000000\t3232\t1.09e-08\n"
        "DN_vs_DP\t1.5\t1400000\t180000000\t3924\t6.6e-08\n")
    d = pfs.read_sweep(str(sw))
    assert d["DN_vs_DP"][0.5]["sig"] == 3232
    assert d["DN_vs_DP"][1.5]["null"] == 180000000

    st = tmp_path / "diff_DN_vs_DP.stats.tsv"
    st.write_text(
        "chrom\tstart\tend\tpval\tpval_adj\tbearing_score\tc\tc\tc\tsignificant\n"
        "chr6\t41198200\t41198400\t0.5\t0.9\t3.0\t0\t0\t0\t0\n")
    s = pfs.read_stats(str(st))
    assert s[("chr6", 41198200)] == (3.0, 0.5)


def _tiny_inputs(tmp_path):
    sw = tmp_path / "floor_sweep.tsv"
    with open(sw, "w") as fh:
        fh.write("comparison\tfloor\tbins_tested\tnull_scores\t"
                 "bins_bh_significant\tmin_pval\n")
        for cmp in ("DN_vs_DP", "DN_vs_S3T3"):
            for f, t, n, s in [(0.1, 5e6, 3e9, 3158), (0.25, 4.6e6, 2e9, 3169),
                               (0.5, 4.2e6, 1.1e9, 3232), (1.0, 2.5e6, 4e8, 3512),
                               (1.5, 1.4e6, 1.8e8, 3924)]:
                fh.write("%s\t%g\t%d\t%d\t%d\t%.6g\n"
                         % (cmp, f, int(t), int(n), s, 12.0 / (n + 1)))
    pbd = tmp_path / "per_bin"
    pbd.mkdir()
    for tag, sig in (("0p5", 0), ("1p5", 1)):
        with gzip.open(pbd / ("DN_vs_DP_floor%s.tsv.gz" % tag), "wt") as fh:
            fh.write("chrom\tstart\tend\tbearing_score\tpval\tpval_adj_bh\tsignificant\n")
            fh.write("chr6\t41500000\t41500200\t1.7\t7e-7\t0.06\t%d\n" % sig)
    std = tmp_path / "pvalue"
    std.mkdir()
    with open(std / "diff_DN_vs_DP.stats.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tpval\tpval_adj\tbearing_score\tc\tc\tc\tsignificant\n")
        fh.write("chr6\t41198200\t41198400\t2.5e-8\t0.5\t3.0\t0\t0\t0\t0\n")
    import json
    qj = tmp_path / "Q_DN.json"
    json.dump({"sample": "DN", "categories": {str(i): "t%d" % i for i in range(1, 7)},
               "Q": {str(i): v for i, v in enumerate([0.084, 0.15, 0.1, 0.18, 0.2, 0.16], 1)},
               "min_signal": 0.1, "score_method": "kl"}, open(qj, "w"))
    ex = tmp_path / "ex.tsv"
    ex.write_text("chrom\tstart\tlabel\nchr6\t41198200\tVb\n")
    return sw, pbd, std, qj, ex


def test_plot_runs_and_pdf_is_reproducible(tmp_path):
    import sys
    sw, pbd, std, qj, ex = _tiny_inputs(tmp_path)

    def run(suffix):
        argv = ["plot_floor_sensitivity.py",
                "--sweep-tsv", str(sw), "--per-bin-dir", str(pbd),
                "--stats-dir", str(std), "--q-json", str(qj),
                "--example-bins", str(ex), "--mechanism-comparison", "DN_vs_DP",
                "--out-pdf", str(tmp_path / ("s%s.pdf" % suffix)),
                "--out-png", str(tmp_path / ("s%s.png" % suffix)),
                "--source-data", str(tmp_path / ("s%s.tsv" % suffix))]
        old = sys.argv
        sys.argv = argv
        try:
            assert pfs.main() == 0
        finally:
            sys.argv = old

    run("1")
    run("2")
    assert (tmp_path / "s1.pdf").exists()
    # byte-identical PDFs on re-run (verification #6)
    assert (tmp_path / "s1.pdf").read_bytes() == (tmp_path / "s2.pdf").read_bytes()
    # source data lists a real flipping bin (the one that goes 0->1)
    src = (tmp_path / "s1.tsv").read_text()
    assert "flipping_bin=chr6:41500000" in src
