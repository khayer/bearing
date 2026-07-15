"""Tests for the |delta insulation| table export in bearing_hic_combined_plot.py.

Covers write_insul_sig_table (the --insul-sig-out feature): TSV + bedgraph
formats, the header provenance line, start/end derivation from center, and the
invariant that `sig == 1` is exactly `delta_abs >= threshold` (the same flag
that drives the figure's significance stars). Export-only: the source table is
unchanged by writing it. ASCII only.
"""

import numpy as np
import pandas as pd
import pytest

import bearing_hic_combined_plot as bhc


def _synth_insul(n=200, res=25000, chrom="chr6", seed=0, shift=0.0):
    rng = np.random.default_rng(seed)
    start = np.arange(n) * res
    end = start + res
    score = rng.normal(loc=shift, scale=1.0, size=n)
    return pd.DataFrame({
        "chrom": chrom, "start": start, "end": end,
        "center": (start + end) / 2, "score": score,
    })


def _read_tsv(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    assert lines[0].startswith("#")
    header = lines[1].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[2:] if ln]
    return lines[0], header, rows


def test_tsv_export_columns_and_sig_matches_threshold(tmp_path):
    res = 25000
    a = _synth_insul(n=200, res=res, seed=1)
    b = _synth_insul(n=200, res=res, seed=2)
    df = bhc.compute_insul_pctile_threshold(a, b, pctile=95.0)

    out = tmp_path / "insul_sig.tsv"
    bhc.write_insul_sig_table(df, str(out), res, "genome", "chr6",
                              "A.bm", "B.bm", 95.0)

    comment, header, rows = _read_tsv(out)
    assert header == ["chrom", "start", "end", "center", "score_A", "score_B",
                      "delta_abs", "threshold", "pctile", "sig"]
    assert len(rows) == len(df)
    # provenance line records the essentials
    for token in ["insul_A=A.bm", "insul_B=B.bm", "resolution=25000",
                  "percentile=95", "scope=genome-wide", "threshold="]:
        assert token in comment

    # start/end derived from center +/- res/2 ; delta_abs == |score_A - score_B|
    for r in rows:
        c = float(r["center"])
        assert int(r["start"]) == int(round(c - res / 2.0))
        assert int(r["end"]) == int(round(c + res / 2.0))
        d = abs(float(r["score_A"]) - float(r["score_B"]))
        assert float(r["delta_abs"]) == pytest.approx(d, rel=1e-4, abs=1e-6)

    # sig == 1 is exactly delta_abs >= threshold (the star-marker flag)
    thr = float(rows[0]["threshold"])
    for r in rows:
        expect = 1 if float(r["delta_abs"]) >= thr else 0
        assert int(r["sig"]) == expect


def test_sig_count_is_top_fraction(tmp_path):
    # ~ (100 - pctile)% of finite bins are flagged.
    res = 25000
    a = _synth_insul(n=1000, res=res, seed=10)
    b = _synth_insul(n=1000, res=res, seed=11)
    pctile = 95.0
    df = bhc.compute_insul_pctile_threshold(a, b, pctile=pctile)
    n_sig = int(df["sig"].sum())
    expected = (100 - pctile) / 100.0 * len(df)
    # percentile ties allow a little slack, but it must be close to 5%.
    assert abs(n_sig - expected) <= 0.02 * len(df)


def test_bedgraph_export(tmp_path):
    res = 25000
    a = _synth_insul(n=50, res=res, seed=3)
    b = _synth_insul(n=50, res=res, seed=4)
    df = bhc.compute_insul_pctile_threshold(a, b, pctile=90.0)

    out = tmp_path / "insul_sig.bg"
    bhc.write_insul_sig_table(df, str(out), res, "chrom", "chr6",
                              "A.bm", "B.bm", 90.0)
    with open(out) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    assert lines[0].startswith("#")
    assert "scope=chrom-restricted (chr6)" in lines[0]
    for ln in lines[1:]:
        f = ln.split("\t")
        assert len(f) == 4                      # chrom start end delta
        int(f[1]); int(f[2]); float(f[3])       # parse cleanly
    assert len(lines) - 1 == len(df)


def test_export_does_not_mutate_source(tmp_path):
    res = 25000
    a = _synth_insul(n=40, res=res, seed=5)
    b = _synth_insul(n=40, res=res, seed=6)
    df = bhc.compute_insul_pctile_threshold(a, b, pctile=95.0)
    before = df.copy(deep=True)
    bhc.write_insul_sig_table(df, str(tmp_path / "x.tsv"), res, "genome",
                              "chr6", "A.bm", "B.bm", 95.0)
    pd.testing.assert_frame_equal(df, before)
