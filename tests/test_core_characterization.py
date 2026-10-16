"""Characterization tests for the shared ``bearing`` package.

These pin the *current* behaviour of the core helpers (parsing, TSV/sheet
loading, command running, path validation, highlight drawing, BED loading,
contact-file resolution and Hi-C matrix loading) so refactors can be checked
for regressions. Inputs come from the small simulated dataset in
``tests/data`` (see tests/data/README.md).
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import pytest

import bearing
from bearing import diagnostics, hic_io, plot_loaders, runner, sheet, tsv, utils, validate
from bearing.track_primitives import add_horizontal_highlights, add_vertical_highlights

DATA = Path(__file__).parent / "data"


# --------------------------------------------------------------------- utils
class TestUtils:
    def test_key_value_strips_and_splits_on_first_equals(self):
        assert utils.parse_key_value_items([" A = x=y ", "B=2", "A=3"], "--m") == {"A": "3", "B": "2"}

    @pytest.mark.parametrize("item,msg", [("A", "expects CONDITION=VALUE"), ("=v", "Invalid --m item"), ("k= ", "Invalid --m item")])
    def test_key_value_errors(self, item, msg):
        with pytest.raises(ValueError, match=msg):
            utils.parse_key_value_items([item], "--m")

    def test_key_value_empty(self):
        assert utils.parse_key_value_items([], "--m") == {}

    def test_sanitize_token(self):
        assert utils.sanitize_token("a b/c.d-e_f:g") == "a_b_c.d-e_f_g"
        assert utils.sanitize_token("a.b", allowed="") == "a_b"
        assert utils.sanitize_token("é1") == "é1"  # unicode alnum kept

    def test_resolve_path(self, tmp_path):
        assert utils.resolve_path("/abs/x", tmp_path) == Path("/abs/x")
        assert utils.resolve_path("a/../b", tmp_path) == (tmp_path / "b").resolve()

    def test_iter_non_comment_lines_keeps_original_line(self):
        src = ["# c\n", "\n", "   \n", "  #x\n", " a\tb \n", "c\n"]
        assert list(utils.iter_non_comment_lines(src)) == [" a\tb \n", "c\n"]

    def test_parse_ucsc_region(self):
        assert utils.parse_ucsc_region("chr6:1,000-2,000") == ("chr6", 1000, 2000)
        with pytest.raises(ValueError):
            utils.parse_ucsc_region("chr6")

    def test_parse_named_region_item(self):
        assert utils.parse_named_region_item(" a = chr1:1-2 ") == {"name": "a", "region": "chr1:1-2"}
        with pytest.raises(ValueError, match="--r must be NAME"):
            utils.parse_named_region_item("chr1:1-2", "--r")
        with pytest.raises(ValueError, match="Invalid region entry"):
            utils.parse_named_region_item("a=")


# ----------------------------------------------------------------- tsv/sheet
class TestTsvAndSheet:
    def test_samples_sheet(self):
        fields, rows = tsv.read_tsv_table(DATA / "sim_samples.tsv")
        assert fields == ["sample", "condition", "bw1", "bw2", "bw3"]
        assert rows == [
            {"sample": "S1", "condition": "WT", "bw1": "s1_atac.bw", "bw2": "", "bw3": "s1_ctcf.bw"},
            {"sample": "S2", "condition": "KO", "bw1": "s2_atac.bw", "bw2": "s2_rna.bw", "bw3": ""},
        ]
        assert sheet.load_rows(DATA / "sim_samples.tsv") == rows
        assert [sheet.extract_bw_values(r) for r in rows] == [["s1_atac.bw", "s1_ctcf.bw"], ["s2_atac.bw", "s2_rna.bw"]]
        assert [sheet.sample_name_from_row(r) for r in rows] == ["S1", "S2"]

    def test_regions_and_required(self):
        rows = tsv.read_tsv_dict_rows(DATA / "sim_regions.tsv", required_columns=["name", "region"])
        assert [utils.parse_ucsc_region(r["region"]) for r in rows] == [("chrSim", 1000, 2000), ("chrSim", 5000, 9000)]
        with pytest.raises(ValueError, match="TSV missing required columns: x, y"):
            tsv.read_tsv_dict_rows(DATA / "sim_regions.tsv", required_columns=["x", "name", "y"])

    def test_extra_and_short_rows(self, tmp_path):
        p = tmp_path / "t.tsv"
        p.write_text("A\tB\n1\t2\t3\t4\n5\n")
        _, rows = tsv.read_tsv_table(p)
        assert rows == [{"a": "1", "b": "2", "": "3,4"}, {"a": "5", "b": ""}]

    def test_empty_tsv(self, tmp_path):
        p = tmp_path / "e.tsv"
        p.write_text("# only comments\n\n")
        with pytest.raises(ValueError, match="empty or missing header"):
            tsv.read_tsv_table(p)

    def test_require_columns(self):
        sheet.require_columns(["a", "b"], ["a"], "ctx")
        with pytest.raises(ValueError, match="ctx missing required columns: c, d"):
            sheet.require_columns(["a"], ["c", "a", "d"], "ctx")
        sheet.require_any_column(["a"], ["z", "a"], "ctx")
        with pytest.raises(ValueError, match="ctx must include at least one of: x, y"):
            sheet.require_any_column(["a"], ["x", "y"], "ctx")

    def test_sample_name_and_bw_edge_cases(self):
        assert sheet.sample_name_from_row({"sample": "", "name": "N"}) == "N"
        assert sheet.sample_name_from_row({}) == ""
        assert sheet.extract_bw_values({"bw": "a", "bw1": "ignored"}) == ["a"]
        assert sheet.extract_bw_values({"bw": "", "bw1": "x", "bw3": "skipped"}) == ["x"]
        assert sheet.extract_bw_values({}) == []


# -------------------------------------------------------- validate/diagnostics
class TestValidation:
    def test_paths(self, tmp_path):
        ok = tmp_path / "ok"
        ok.write_text("")
        validate.ensure_paths_exist([ok])
        with pytest.raises(FileNotFoundError, match=r"missing bigwig: .*nope1"):
            validate.ensure_paths_exist([ok, tmp_path / "nope1", tmp_path / "nope2"], context="bigwig")
        with pytest.raises(FileNotFoundError, match="missing input"):
            validate.ensure_paths_exist([tmp_path / "x"])
        diagnostics.require_files([ok], "f")
        with pytest.raises(FileNotFoundError, match="missing map"):
            diagnostics.require_mapping_paths({"a": "ok", "b": "gone"}, lambda v: tmp_path / v, "map")
        diagnostics.require_mapping_paths({"a": "ok"}, lambda v: tmp_path / v, "map")


# -------------------------------------------------------------------- runner
class TestRunner:
    def test_format_command(self):
        assert runner.format_command(["echo", "a b", 3]) == "echo 'a b' 3"

    def test_emit_and_dry_run(self, capsys):
        assert runner.execute_command(["echo", "x"], run=False, label="L") == 0
        assert capsys.readouterr().out == "[L]\n  $ echo x\n"
        assert runner.execute_commands([["a"], ["b c"]], run=False) == (2, 0)
        assert capsys.readouterr().out == "  $ a\n  $ 'b c'\n"

    def test_execute_command_real(self, capsys):
        assert runner.execute_command([sys.executable, "-c", "raise SystemExit(3)"], run=True) == 3
        assert capsys.readouterr().out.startswith("  $ ")

    @pytest.mark.parametrize("threads", [1, 3])
    def test_execute_commands_real(self, threads, capsys):
        ok = [sys.executable, "-c", "pass"]
        bad = [sys.executable, "-c", "raise SystemExit(1)"]
        assert runner.execute_commands([ok, bad, bad], run=True, threads=threads) == (3, 2)
        out = capsys.readouterr().out
        word = "running" if threads == 1 else "queued"
        assert f"[1/3] {word}\n" in out and f"[3/3] {word}\n" in out
        assert out.count("  $ ") == 3


# ---------------------------------------------------------- track primitives
class FakeAx:
    def __init__(self):
        self.calls = []

    def axvspan(self, a, b, **kw):
        self.calls.append(("v", a, b, kw))

    def axhspan(self, a, b, **kw):
        self.calls.append(("h", a, b, kw))


def _g2a(pos, s, e):
    return (pos - s) / (e - s)


class TestTrackPrimitives:
    def test_horizontal_clips_to_region(self):
        ax = FakeAx()
        add_horizontal_highlights(ax, [(-50, 50, "red"), (90, 500, "blue")], 0, 100, _g2a)
        assert ax.calls == [
            ("v", 0.0, 0.5, {"color": "red", "alpha": 0.18, "zorder": 0}),
            ("v", 0.9, 1.0, {"color": "blue", "alpha": 0.18, "zorder": 0}),
        ]

    def test_vertical_and_empty(self):
        ax = FakeAx()
        add_vertical_highlights(ax, [(10, 20, "g")], 0, 100, _g2a, alpha=0.5, zorder=2)
        add_vertical_highlights(ax, None, 0, 100, _g2a)
        add_horizontal_highlights(ax, [], 0, 100, _g2a)
        assert ax.calls == [("h", 0.1, 0.2, {"color": "g", "alpha": 0.5, "zorder": 2})]


# -------------------------------------------------------------- plot loaders
class TestPlotLoaders:
    @pytest.mark.parametrize("name,fmt", [("a.cool", "cool"), ("A.MCOOL", "cool"), ("x.hic", "hic")])
    def test_detect(self, name, fmt):
        assert plot_loaders.detect_hic_format(name) == fmt

    def test_detect_bad(self):
        with pytest.raises(ValueError, match="Cannot detect Hi-C format"):
            plot_loaders.detect_hic_format("x.txt")

    @pytest.mark.parametrize("name,expected", [
        ("d/s_bs_5000.cool", "d/s_bs_10000.cool"),
        ("s_RES_5000.cool", "s_RES_10000.cool"),
        ("s_resolution_1.cool", "s_resolution_10000.cool"),
        ("s_bin_1.cool", "s_bin_10000.cool"),
        ("s_binsize_1.cool", "s_binsize_10000.cool"),
        ("s_5000.cool", None),
        ("s_bs_5000.mcool", None),
    ])
    def test_cool_resolution_variant(self, name, expected):
        got = plot_loaders.cool_resolution_variant(name, 10000.0)
        assert got == (Path(expected) if expected else None)

    def test_resolve_contacts(self, tmp_path):
        r = plot_loaders.resolve_contacts_for_region_resolution
        a, b = tmp_path / "a_bs_5000.cool", tmp_path / "b_bs_5000.cool"
        assert r("a.hic", "b.cool", 1)[2].endswith("(non-cool format)")
        assert r("a.mcool", "b.cool", 1)[2].endswith("(.mcool uses internal resolutions)")
        assert r("a.hic", "b.hic", 1, hic_fmt="cool")[2].endswith("(unsupported extension)")
        assert r("a.cool", "b.cool", 1)[2].endswith("(no resolution token in filename)")
        assert r(a, b, 100) == (a, b, "no matching .cool files found; keeping original pair")
        (tmp_path / "a_bs_100.cool").write_text("")
        assert r(a, b, 100) == (a, b, "found only one matching .cool file; keeping original pair")
        (tmp_path / "b_bs_100.cool").write_text("")
        assert r(a, b, 100) == (tmp_path / "a_bs_100.cool", tmp_path / "b_bs_100.cool", "using resolution-matched .cool files")
        assert r(str(a), str(b), 100)[2] == "using resolution-matched .cool files"

    def test_load_bed_for_region(self):
        feats = plot_loaders.load_bed_for_region(DATA / "sim_features.bed", "chrSim", 200, 1300)
        assert [(f["name"], f["item_rgb"], f["bed_columns"]) for f in feats] == [
            ("", None, 3), ("geneA", None, 4), ("geneB", None, 6), ("enh1", (255, 0, 0), 9),
            ("enh2", (0, 255, 0), 9), ("badrgb", None, 9), ("junkrgb", None, 9),
        ]
        assert feats[0] == {"chrom": "chrSim", "start": 100, "end": 200, "name": "", "score": "",
                            "strand": None, "item_rgb": None, "bed_columns": 3}
        assert feats[2]["score"] == "7" and feats[2]["strand"] == "+"

    def test_load_bed_missing_and_other_chrom(self, tmp_path):
        assert plot_loaders.load_bed_for_region(tmp_path / "nope.bed", "chrSim", 0, 10) == []
        assert [f["name"] for f in plot_loaders.load_bed_for_region(DATA / "sim_features.bed", "chrOther", 0, 10**6)] == ["other"]


# -------------------------------------------------------------------- hic_io
def _sim_records():
    rows = (DATA / "sim_contacts.tsv").read_text().splitlines()[1:]
    return [tuple(float(x) for x in r.split("\t")) for r in rows]


def _expected(start, end, res):
    n = (end - start) // res + 1
    m = np.zeros((n, n))
    for x, y, c in _sim_records():
        i, j = int((x - start) // res), int((y - start) // res)
        if 0 <= i < n and 0 <= j < n:
            m[i, j] = m[j, i] = c
    return np.log1p(m)


class _Rec:
    def __init__(self, x, y, c):
        self.binX, self.binY, self.counts = int(x), int(y), c


def _install_fake_hicstraw(monkeypatch, resolutions, norms_ok):
    calls = []

    class ZD:
        def __init__(self, norm):
            self.norm = norm

        def getRecords(self, *a):
            if self.norm not in norms_ok:
                raise RuntimeError("no norm")
            return [_Rec(*r) for r in _sim_records()]

    class HiCFile:
        def __init__(self, path):
            pass

        def getResolutions(self):
            return resolutions

        def getMatrixZoomData(self, c1, c2, kind, norm, unit, res):
            calls.append((c1, norm, res))
            return ZD(norm)

    monkeypatch.setitem(sys.modules, "hicstraw", types.SimpleNamespace(HiCFile=HiCFile))
    return calls


class TestHicIO:
    def test_hic_kr(self, monkeypatch, capsys):
        calls = _install_fake_hicstraw(monkeypatch, [100, 1000], {"KR"})
        mat, res = hic_io.load_contact_matrix("x.hic", "chrSim", 0, 1000, 100)
        assert res == 100 and calls == [("Sim", "KR", 100)]
        np.testing.assert_allclose(mat, _expected(0, 1000, 100))
        assert capsys.readouterr().out == ""

    def test_hic_fallback_resolution_and_norm(self, monkeypatch, capsys):
        calls = _install_fake_hicstraw(monkeypatch, [100, 1000], {"VC_SQRT"})
        mat, res = hic_io.load_hic_matrix("x.hic", "chrSim", 0, 1000, 120)
        assert res == 100 and [c[1] for c in calls] == ["KR", "SCALE", "VC_SQRT"]
        out = capsys.readouterr().out
        assert "resolution 120 not in .hic file. Using 100" in out and "using VC_SQRT instead" in out
        assert mat.shape == (11, 11)

    def test_hic_no_norm_exits(self, monkeypatch):
        _install_fake_hicstraw(monkeypatch, [100], set())
        with pytest.raises(SystemExit, match="Could not load Hi-C matrix"):
            hic_io.load_hic_matrix("x.hic", "chrSim", 0, 1000, 100)

    def test_missing_libs_exit(self, monkeypatch):
        monkeypatch.setitem(sys.modules, "hicstraw", None)
        monkeypatch.setitem(sys.modules, "cooler", None)
        with pytest.raises(SystemExit, match="hicstraw is required"):
            hic_io.load_hic_matrix("x.hic", "chr1", 0, 1, 1)
        with pytest.raises(SystemExit, match="cooler is required"):
            hic_io.load_cool_matrix("x.cool", "chr1", 0, 1, 1)

    def _install_fake_cooler(self, monkeypatch, ok_uris, balanced=True, binsize=None):
        opened = []
        raw = np.array([[1.0, np.nan], [np.nan, 3.0]])

        class Matrix:
            def __init__(self, balance):
                self.balance = balance

            def fetch(self, region):
                opened.append(("fetch", region, self.balance))
                if self.balance and not balanced:
                    raise KeyError("weight")
                return raw.copy()

        class Cooler:
            def __init__(self, uri):
                if ok_uris is not None and uri not in ok_uris:
                    raise KeyError(uri)
                opened.append(("open", uri))
                self.binsize = binsize

            def matrix(self, balance):
                return Matrix(balance)

        monkeypatch.setitem(sys.modules, "cooler", types.SimpleNamespace(Cooler=Cooler))
        return opened, np.log1p(np.array([[1.0, 0.0], [0.0, 3.0]]))

    def test_cool_plain(self, monkeypatch):
        opened, exp = self._install_fake_cooler(monkeypatch, None, binsize=500)
        mat, res = hic_io.load_contact_matrix("x.cool", "chrSim", 0, 1000, 100)
        assert res == 500 and mat.dtype == np.float64
        np.testing.assert_allclose(mat, exp)
        assert opened == [("open", "x.cool"), ("fetch", "chrSim:0-1000", True)]

    def test_cool_unbalanced_no_binsize(self, monkeypatch, capsys):
        opened, _ = self._install_fake_cooler(monkeypatch, None, balanced=False)
        _, res = hic_io.load_cool_matrix("x.cool", "chrSim", 0, 1000, 100)
        assert res == 100 and opened[-1] == ("fetch", "chrSim:0-1000", False)
        assert "ICE balancing weights not found" in capsys.readouterr().out

    def test_mcool_direct(self, monkeypatch):
        opened, _ = self._install_fake_cooler(monkeypatch, {"x.mcool::/resolutions/100"}, binsize=999)
        _, res = hic_io.load_cool_matrix("x.mcool", "chrSim", 0, 1000, 100)
        assert res == 100 and opened[0] == ("open", "x.mcool::/resolutions/100")

    def test_mcool_closest_and_missing(self, monkeypatch, tmp_path, capsys):
        h5py = pytest.importorskip("h5py")
        p = tmp_path / "x.mcool"
        with h5py.File(p, "w") as f:
            for r in ("50", "200"):
                f.create_group(f"resolutions/{r}")
        opened, _ = self._install_fake_cooler(monkeypatch, {f"{p}::/resolutions/50"})
        _, res = hic_io.load_cool_matrix(p, "chrSim", 0, 1000, 70)
        assert res == 50 and "Using 50 instead" in capsys.readouterr().out
        with pytest.raises(SystemExit, match="Could not open"):
            hic_io.load_cool_matrix(tmp_path / "none.mcool", "chrSim", 0, 1000, 70)


def test_package_exports():
    for name in bearing.__all__:
        assert callable(getattr(bearing, name))


def test_null_min_bytes_default_and_env(monkeypatch):
    import bearing_pvalue
    monkeypatch.delenv("BEARING_NULL_MIN_BYTES", raising=False)
    assert bearing_pvalue._default_null_min_bytes() == 1_000_000
    monkeypatch.setenv("BEARING_NULL_MIN_BYTES", "0")
    assert bearing_pvalue._default_null_min_bytes() == 0
    monkeypatch.setenv("BEARING_NULL_MIN_BYTES", "abc")
    with pytest.raises(SystemExit, match="must be an integer"):
        bearing_pvalue._default_null_min_bytes()
