# ----------------------------------------------------------------------
# File     : tests/test_rna_dominance_transforms.py
# Version  : 1.1.0
# Date     : 2026-11-08 (tail statistics test; 1.0.0 2026-11-07)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Tests for the 1.1.0 transform variants of dev/rna_dominance_diagnostic.py
# (log1p_rna, asinh_auto), the data-dependent cofactors, the lead share and
# the known-region check. Pure numpy; bigwig_to_qcat is only needed for nzq.
# ASCII only.
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "dev"))
import rna_dominance_diagnostic as rdd  # noqa: E402


def _data():
    a = np.array([[0.0, 50.0, 0.0, 1.0, 1.0, 1.0],
                  [0.2, 0.0, 0.0, 0.1, 0.3, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]], dtype=np.float32)
    b = a.copy()
    b[0, 1] = 2.0
    return {"A_rep1": ("A", {"chr6": ([(0, 200), (200, 400), (400, 600)], a)}),
            "B_rep1": ("B", {"chr6": ([(0, 200), (200, 400), (400, 600)], b)})}


def test_cofactors_are_pooled_nonzero_medians():
    cof = rdd.track_cofactors(_data())
    assert np.isclose(cof[1], np.median([50.0, 2.0]))
    assert np.isclose(cof[3], np.median([1.0, 0.1, 1.0, 0.1]))
    assert cof[2] == 1.0  # all-zero track


def test_transforms_keep_zero_and_the_mask_signal():
    cof = rdd.track_cofactors(_data())
    raw = _data()["A_rep1"][1]["chr6"][1]
    for v in ("log1p_rna", "asinh_auto"):
        x, pre = rdd.variant_matrix(raw, v, 1, 2, cof)
        assert np.all(x[2] == 0)                       # zero stays zero (on/off unchanged)
        assert np.allclose(pre, np.clip(raw, 0, None))  # mask uses untransformed signal


def test_log1p_rna_shrinks_rna_dominance_only_on_rna():
    cof = rdd.track_cofactors(_data())
    raw = _data()["A_rep1"][1]["chr6"][1]
    x, _ = rdd.variant_matrix(raw, "log1p_rna", 1, 2, cof)
    assert np.allclose(x[:, [0, 3, 4, 5]], raw[:, [0, 3, 4, 5]])
    share_before = raw[0, 1] / raw[0].sum()
    share_after = x[0, 1] / x[0].sum()
    assert share_after < share_before


def test_production_unchanged():
    raw = _data()["A_rep1"][1]["chr6"][1]
    x, pre = rdd.variant_matrix(raw, "production", 1, 2, None)
    assert np.allclose(x, np.clip(raw, 0, None)) and np.allclose(pre, x)


def test_lead_share_and_region_check():
    d = np.array([[0, 3.0, 0, 0, 0, 0], [0.5, 0, 0, 0.5, 0, 0], [0, 0, 0, 0, 0, 0]])
    assert np.isclose(rdd.median_lead_share(d, np.array([0])), 1.0)
    assert np.isclose(rdd.median_lead_share(d, np.array([1])), 0.5)
    l1 = np.abs(d).sum(axis=1)
    tested = np.where(l1 >= 0.5)[0]
    order = tested[np.argsort(-l1[tested], kind="stable")]
    regs = rdd.parse_regions("R=chr6:0-250")
    out = rdd.region_check(d, l1, tested, order, [("top0.5", order[:1])], regs,
                           ["chr6"], np.zeros(3, dtype=np.int32),
                           np.array([0, 200, 400]), np.array([200, 400, 600]),
                           ["t%d" % i for i in range(6)])
    assert out[0]["tested"] == 2 and out[0]["in_top0.5"] == 1
    assert out[0]["best_rank_pct"] == 0.0 and out[0]["lead"] == "t1"


def test_tail_stats_per_condition_and_pooled():
    rows = rdd.tail_stats(_data(), ["ATAC", "RNA+", "RNA-", "CTCF", "Cohesin",
                                    "H3K27ac"], "A", "B")
    r = {(x["track"], x["group"]): x for x in rows}
    # RNA+: A has one nonzero value (50), B one (2); 3 bins per sample
    assert r[("RNA+", "A")]["frac_nonzero"] == np.float64(1 / 3)
    assert r[("RNA+", "A")]["q50"] == 50.0
    assert r[("RNA+", "pooled")]["max"] == 50.0
    assert r[("RNA+", "pooled")]["bins"] == 6
    # an all-zero track gives NaN quantiles, not an error
    assert np.isnan(r[("RNA-", "pooled")]["q50"])
    # q99/q50 is the dynamic-range ratio
    c = r[("CTCF", "pooled")]
    assert c["q99_over_q50"] == c["q99"] / c["q50"]
    assert len(rdd.tail_text(rows)) == 2 + len(rows)
