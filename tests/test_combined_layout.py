"""Tests for the opt-in figure enhancements in bearing_hic_combined_plot.py:
loop-arc rows, slim RGB overview, gene-label curation, and the layout's
default-unchanged contract. These exercise the pure helpers/layout without
needing Hi-C data. ASCII only.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bearing_hic_combined_plot as bhc


GENES = [
    (100, 200, "Trbv13", "+"),
    (300, 400, "Trbd2", "-"),
    (500, 600, "Gm12345", "+"),
    (700, 800, "Sval2", "-"),
]


def test_filter_genes_empty_whitelist_is_identity():
    assert bhc._filter_genes_by_name(GENES, None) is GENES
    assert bhc._filter_genes_by_name(GENES, []) is GENES


def test_filter_genes_exact_and_substring_case_insensitive():
    out = bhc._filter_genes_by_name(GENES, ["trbv13", "trbd"])
    names = {g[2] for g in out}
    assert names == {"Trbv13", "Trbd2"}   # exact + substring, case-insensitive
    assert bhc._filter_genes_by_name(GENES, ["nope"]) == []


def _layout_keys(**kw):
    base = dict(show_rgb=False, show_hic_list=["A", "B", "diff"],
                hic_height_in=3.0, has_insul=True, num_beds=0, bed_styles=[],
                num_decomp=2, has_pval_fill=False)
    base.update(kw)
    fig, axes = bhc._combined_figure_layout(**base)
    keys = set(axes.keys())
    plt.close(fig)
    return keys


def test_default_layout_has_no_loop_or_rgb_rows():
    """Contract: without the new flags the panel set is unchanged."""
    keys = _layout_keys()
    assert "loops_a" not in keys and "loops_b" not in keys
    assert "rgb" not in keys
    # core panels still present
    for k in ("hic_A", "hic_B", "hic_diff", "insul", "delta_insul",
              "qcat_a", "qcat_b", "diff", "gene", "legend"):
        assert k in keys


def test_loop_rows_added_when_requested():
    keys = _layout_keys(has_loops_a=True, has_loops_b=True)
    assert "loops_a" in keys and "loops_b" in keys


def test_rgb_overview_row_added():
    keys = _layout_keys(show_rgb=True, rgb_slim=True)
    assert "rgb" in keys
