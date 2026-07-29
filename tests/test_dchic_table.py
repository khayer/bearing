"""Unit tests for the dcHiC compartment-table Python consumers.

These lock the pure logic that decides whether the supplementary table is
CORRECT -- the compartment sign call, the padj column priority, and the
experiment-column detection -- none of which need dcHiC or real data to test.

The detect_experiment_columns test also guards a regression: that helper uses
re.search, and dchic_region_result.py once used it without importing re, which
crashed the region-extract step at runtime (main() reaches it unconditionally).
Importing and exercising it here fails loudly if that import ever goes missing.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import build_compartment_supp_table as B
import dchic_region_result as R


# ---- dchic_region_result.py ------------------------------------------------

def test_parse_region_strips_commas():
    assert R.parse_region("chr6:40,872,100-41,558,435") == ("chr6", 40872100, 41558435)
    assert R.parse_region("chr6:39000000-44000000") == ("chr6", 39000000, 44000000)


def test_detect_experiment_columns_drops_meta_and_replicates():
    # exercises re.search -- guards the missing-`import re` regression.
    # Header mirrors a raw dcHiC bedGraph: metadata + per-condition + per-rep.
    header = ["chr", "start", "end", "DN", "DP", "S3T3",
              "padj", "DN_rep1", "DP_rep2"]
    assert R.detect_experiment_columns(header) == ["DN", "DP", "S3T3"]


def test_detect_experiment_columns_two_condition_run():
    header = ["chr", "start", "end", "DN", "ProB", "padj"]
    cols = R.detect_experiment_columns(header)
    assert cols == ["DN", "ProB"]
    assert len(cols) == 2          # -> pairwise (not multi-condition/global)


# ---- build_compartment_supp_table.py ---------------------------------------

def test_call_sign_convention():
    # positive PC1 = A compartment, non-positive = B (zero -> B)
    assert B.call(0.42) == "A"
    assert B.call(-0.42) == "B"
    assert B.call(0.0) == "B"


def test_pc_col_found_and_missing():
    header = ["chr", "start", "end", "DN_PC", "S3T3_PC", "padj"]
    assert B.pc_col(header, "DN") == "DN_PC"
    assert B.pc_col(header, "S3T3") == "S3T3_PC"
    with pytest.raises(SystemExit):
        B.pc_col(header, "EBKO")          # absent -> fail loud


def test_padj_col_priority_and_missing():
    # priority: padj_global > padj_lymphoid > padj
    assert B.padj_col(["chr", "padj", "padj_global"]) == "padj_global"
    assert B.padj_col(["chr", "padj", "padj_lymphoid"]) == "padj_lymphoid"
    assert B.padj_col(["chr", "padj"]) == "padj"
    with pytest.raises(SystemExit):
        B.padj_col(["chr", "start", "end"])   # no padj column -> fail loud


def test_load_contrast_parses_tsv(tmp_path):
    p = tmp_path / "Tcrb_compartment_dchic_S3T3_100kb.tsv"
    p.write_text(
        "# comment line ignored\n"
        "chr\tstart\tend\tDN_PC\tS3T3_PC\tpadj_global\tsign_change\tglobal_significant\n"
        "chr6\t41000000\t41100000\t0.30\t-0.20\t0.01\t1\t1\n"
        "chr6\t41100000\t41200000\t-0.10\t-0.15\t0.80\t0\t0\n"
    )
    d = B.load_contrast(str(p), "S3T3")
    k = ("chr6", 41000000, 41100000)
    assert k in d
    assert d[k]["DN"] == pytest.approx(0.30)
    assert d[k]["S3T3"] == pytest.approx(-0.20)
    assert d[k]["padj"] == pytest.approx(0.01)
    # sign call from the parsed PCs: DN=A, S3T3=B -> a real pairwise flip
    assert B.call(d[k]["DN"]) == "A" and B.call(d[k]["S3T3"]) == "B"
