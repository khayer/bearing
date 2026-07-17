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


# --- regression tests for the loop-loading + RGB fixes ---------------------

def test_load_loops_chrom_and_delimiter_tolerant(tmp_path):
    from bearing_hic_plot import load_loops
    p = tmp_path / "loops.bedpe"
    p.write_text(
        "chrom1 start1 end1 chrom2 start2 end2 score\n"   # header, space-delimited
        "6 41010000 41012000 6 41060000 41062000 AnchorA\n"  # chrom '6' vs region 'chr6'
        "6\t41020000\t41022000\t6\t41090000\t41092000\t5.0\n"  # tab-delimited
        "7 100 200 7 300 400 1\n")                             # different chrom -> excluded
    loops = load_loops(str(p), "chr6", 41000000, 41100000)
    assert len(loops) == 2


def test_make_rgb_red_green_joint_match_is_yellow():
    import numpy as np
    from bearing_hic_plot import make_rgb_hic
    a = np.array([[10.0, 0.0], [0.0, 4.0]])
    rgb = make_rgb_hic(a, a.copy(), palette="red-green", joint_norm=True)
    px = rgb[0, 0]
    assert px[0] == px[1] and px[2] == 0          # equal contact -> yellow
    r = make_rgb_hic(np.array([[10.0]]), np.array([[0.0]]),
                     palette="red-green", joint_norm=True)[0, 0]
    assert r[0] > 0 and r[1] == 0 and r[2] == 0   # A-only -> red


def test_loops_with_offframe_partner_still_draw():
    """A loop with one anchor in-frame and its partner outside the plotted
    region must still draw a (partial) arc, not be dropped."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc
    from bearing_hic_plot import draw_loops_horizontal
    rs, re_ = 40_400_000, 42_300_000
    loops = [
        (41_000_000, 41_010_000, 41_500_000, 41_510_000, 1.0),  # both in-frame
        (41_200_000, 41_210_000, 43_000_000, 43_010_000, 1.0),  # partner off right
        (39_000_000, 39_010_000, 41_400_000, 41_410_000, 1.0),  # partner off left
    ]
    fig, ax = plt.subplots()
    draw_loops_horizontal(ax, loops, rs, re_, label="loops")
    arcs = [p for p in ax.patches if isinstance(p, Arc)]
    plt.close(fig)
    assert len(arcs) == 3
