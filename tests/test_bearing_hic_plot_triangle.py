import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import importlib.util


HAS_MPL = importlib.util.find_spec("matplotlib") is not None
HAS_NUMPY = importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(HAS_MPL and HAS_NUMPY, "matplotlib/numpy not installed")
class TriangleFigureTests(unittest.TestCase):
    def test_make_triangle_figure_with_synthetic_inputs(self):
        import numpy as np

        from bearing_hic_plot_triangle import make_figure_triangle

        hic_a = np.array(
            [
                [10, 8, 5, 2],
                [8, 10, 6, 3],
                [5, 6, 10, 6],
                [2, 3, 6, 10],
            ],
            dtype=float,
        )
        hic_b = np.array(
            [
                [9, 7, 4, 1],
                [7, 9, 5, 2],
                [4, 5, 9, 5],
                [1, 2, 5, 9],
            ],
            dtype=float,
        )

        positions = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4000)]
        scores_a = np.array(
            [
                [0.6, 0.4],
                [0.2, 0.8],
                [0.7, 0.3],
                [0.5, 0.5],
            ],
            dtype=float,
        )
        scores_b = np.array(
            [
                [0.5, 0.5],
                [0.3, 0.7],
                [0.6, 0.4],
                [0.4, 0.6],
            ],
            dtype=float,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "triangle.png"

            with patch(
                "bearing_hic_plot_triangle.load_contact_matrix",
                side_effect=[(hic_a, 1000), (hic_b, 1000)],
            ), patch(
                "bearing_hic_plot_triangle.load_qcat_scores",
                side_effect=[
                    (positions, scores_a, 2),
                    (positions, scores_b, 2),
                ],
            ):
                make_figure_triangle(
                    hic_a_path="a.cool",
                    hic_b_path="b.cool",
                    qcat_a_path="a.qcat.bgz",
                    qcat_b_path="b.qcat.bgz",
                    region_str="chr1:0-4000",
                    resolution=1000,
                    out_path=str(out_path),
                    label_a="Reference",
                    label_b="Current",
                    categories=[("State1", "#00aa55"), ("State2", "#d33636")],
                    rgb_hic=False,
                )

            self.assertTrue(out_path.exists())
            self.assertGreater(out_path.stat().st_size, 0)

    def test_triangle_passes_diff_pval_overlay_inputs(self):
        import numpy as np

        from bearing_hic_plot_triangle import make_figure_triangle

        hic_a = np.array([[10, 8], [8, 10]], dtype=float)
        hic_b = np.array([[9, 7], [7, 9]], dtype=float)

        positions = [(0, 1000), (1000, 2000)]
        diff_positions = [(0, 1000), (1000, 2000)]
        scores_a = np.array([[0.6, 0.4], [0.2, 0.8]], dtype=float)
        scores_b = np.array([[0.5, 0.5], [0.3, 0.7]], dtype=float)
        diff_scores = np.array([[1.0, -2.0], [-3.0, 4.0]], dtype=float)
        pval_positions = [500, 1500]
        pval_values = np.array([3.0, -4.0], dtype=float)

        captured = []

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "triangle_overlay.png"

            with patch(
                "bearing_hic_plot_triangle.load_contact_matrix",
                side_effect=[(hic_a, 1000), (hic_b, 1000)],
            ), patch(
                "bearing_hic_plot_triangle.load_qcat_scores",
                side_effect=[
                    (positions, scores_a, 2),
                    (positions, scores_b, 2),
                    (diff_positions, diff_scores, 2),
                ],
            ), patch(
                "bearing_hic_plot_triangle.load_pval_track_values",
                return_value=(pval_positions, pval_values),
            ), patch(
                "bearing_hic_plot_triangle.draw_pval_diff_horizontal",
                side_effect=lambda *args, **kwargs: captured.append(kwargs),
            ):
                make_figure_triangle(
                    hic_a_path="a.cool",
                    hic_b_path="b.cool",
                    qcat_a_path="a.qcat.bgz",
                    qcat_b_path="b.qcat.bgz",
                    region_str="chr1:0-2000",
                    resolution=1000,
                    out_path=str(out_path),
                    label_a="Reference",
                    label_b="Current",
                    categories=[("State1", "#6495ed"), ("State2", "#d33636")],
                    diff_qcat_path="diff.qcat.bgz",
                    pval_diff_path="diff.bw",
                    rgb_hic=False,
                )

            self.assertTrue(out_path.exists())
            self.assertTrue(captured)
            self.assertEqual(captured[0]["diff_score_positions"], diff_positions)
            self.assertTrue(np.array_equal(captured[0]["diff_score_matrix"], -diff_scores))
            self.assertEqual(captured[0]["categories"], [("State1", "#6495ed"), ("State2", "#d33636")])

    def test_triangle_reads_pval_track_from_stats_tsv(self):
        import numpy as np

        from bearing_hic_plot_triangle import make_figure_triangle

        hic_a = np.array([[10, 8], [8, 10]], dtype=float)
        hic_b = np.array([[9, 7], [7, 9]], dtype=float)

        positions = [(0, 1000), (1000, 2000)]
        scores_a = np.array([[0.6, 0.4], [0.2, 0.8]], dtype=float)
        scores_b = np.array([[0.5, 0.5], [0.3, 0.7]], dtype=float)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            out_path = tmpdir / "triangle_stats_tsv.png"
            pval_tsv = tmpdir / "diff.stats.tsv"
            pval_tsv.write_text(
                "chrom\tstart\tend\tbearing_score\tnull_method\tscore_normalised\tpval\tpval_adj_bh\tsignificant_fdr0.05\tdirection\n"
                "chr1\t0\t1000\t1.0\tempirical\tfalse\t1.0e-03\t1.0e-03\t1\t+\n"
                "chr1\t1000\t2000\t-2.0\tempirical\tfalse\t1.0e-04\t1.0e-04\t1\t-\n"
            )

            captured = []

            with patch(
                "bearing_hic_plot_triangle.load_contact_matrix",
                side_effect=[(hic_a, 1000), (hic_b, 1000)],
            ), patch(
                "bearing_hic_plot_triangle.load_qcat_scores",
                side_effect=[
                    (positions, scores_a, 2),
                    (positions, scores_b, 2),
                ],
            ), patch(
                "bearing_hic_plot_triangle.draw_pval_diff_horizontal",
                side_effect=lambda *args, **kwargs: captured.append(kwargs),
            ):
                make_figure_triangle(
                    hic_a_path="a.cool",
                    hic_b_path="b.cool",
                    qcat_a_path="a.qcat.bgz",
                    qcat_b_path="b.qcat.bgz",
                    region_str="chr1:0-2000",
                    resolution=1000,
                    out_path=str(out_path),
                    label_a="Reference",
                    label_b="Current",
                    categories=[("State1", "#6495ed"), ("State2", "#d33636")],
                    pval_diff_path=str(pval_tsv),
                    rgb_hic=False,
                )

            self.assertTrue(out_path.exists())
            self.assertGreater(out_path.stat().st_size, 0)
            self.assertTrue(captured)

    def test_triangle_uses_overlay_helper_when_requested(self):
        import numpy as np

        from bearing_hic_plot_triangle import make_figure_triangle

        hic_a = np.array([[10, 8], [8, 10]], dtype=float)
        hic_b = np.array([[9, 7], [7, 9]], dtype=float)

        positions = [(0, 1000), (1000, 2000)]
        scores_a = np.array([[0.6, 0.4], [0.2, 0.8]], dtype=float)
        scores_b = np.array([[0.5, 0.5], [0.3, 0.7]], dtype=float)
        pval_positions = [500, 1500]
        pval_values = np.array([3.0, 4.0], dtype=float)

        overlay_calls = []

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "triangle_overlay_on.png"

            with patch(
                "bearing_hic_plot_triangle.load_contact_matrix",
                side_effect=[(hic_a, 1000), (hic_b, 1000)],
            ), patch(
                "bearing_hic_plot_triangle.load_qcat_scores",
                side_effect=[
                    (positions, scores_a, 2),
                    (positions, scores_b, 2),
                ],
            ), patch(
                "bearing_hic_plot_triangle.load_pval_track_values",
                side_effect=[
                    (pval_positions, pval_values),
                    (pval_positions, pval_values),
                ],
            ), patch(
                "bearing_hic_plot_triangle.draw_epilogos_with_pval_horizontal",
                side_effect=lambda *args, **kwargs: overlay_calls.append(kwargs),
            ):
                make_figure_triangle(
                    hic_a_path="a.cool",
                    hic_b_path="b.cool",
                    qcat_a_path="a.qcat.bgz",
                    qcat_b_path="b.qcat.bgz",
                    region_str="chr1:0-2000",
                    resolution=1000,
                    out_path=str(out_path),
                    label_a="Reference",
                    label_b="Current",
                    categories=[("State1", "#6495ed"), ("State2", "#d33636")],
                    pval_a_path="a.stats.tsv",
                    pval_b_path="b.stats.tsv",
                    pval_overlay=True,
                    rgb_hic=False,
                )

            self.assertTrue(out_path.exists())
            self.assertGreater(out_path.stat().st_size, 0)
            self.assertEqual(len(overlay_calls), 2)
            self.assertEqual(overlay_calls[0]["pval_positions"], pval_positions)
            self.assertEqual(overlay_calls[1]["pval_positions"], pval_positions)


class BatchTriangleRoutingTests(unittest.TestCase):
    def test_triangle_flag_routes_child_command(self):
        import batch_bearing_hic_plots as batch

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sheet = root / "sheet.tsv"
            regions = root / "regions.tsv"
            results = root / "results"
            results.mkdir()
            outdir = root / "out"
            contacts = root / "contacts"
            contacts.mkdir()

            (root / "ref.qcat.bgz").write_text("ref")
            (root / "cur.qcat.bgz").write_text("cur")
            (results / "ref1.neglog10p.bw").write_text("")
            (results / "cur1.neglog10p.bw").write_text("")
            (contacts / "ref.cool").write_text("")
            (contacts / "cur.cool").write_text("")
            sheet.write_text(
                "sample\tcondition\treplicate\tout\n"
                f"ref1\treference\t1\t{root / 'ref.qcat.bgz'}\n"
                f"cur1\tcurrent\t1\t{root / 'cur.qcat.bgz'}\n"
            )
            regions.write_text("name\tregion\nregion1\tchr1:0-4000\n")

            captured = []

            argv = [
                "batch_bearing_hic_plots.py",
                "--sheet", str(sheet),
                "--regions-file", str(regions),
                "--reference-condition", "reference",
                "--contact", f"reference={contacts / 'ref.cool'}",
                "--contact", f"current={contacts / 'cur.cool'}",
                "--results-dir", str(results),
                "--outdir", str(outdir),
                "--triangle",
            ]

            with patch.object(sys, "argv", argv), patch.object(batch, "execute_command", side_effect=lambda cmd, run, label: captured.append(cmd) or 0), patch.object(batch, "_load_regions", return_value=[]):
                rc = batch.main()

            self.assertEqual(rc, 0)
            self.assertTrue(captured)
            self.assertEqual(Path(captured[0][1]).name, "bearing_hic_plot_triangle.py")
            self.assertNotIn("--triangle", captured[0])


if __name__ == "__main__":
    unittest.main()
