import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


HAS_MPL = importlib.util.find_spec("matplotlib") is not None
HAS_NUMPY = importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(HAS_MPL and HAS_NUMPY, "matplotlib/numpy not installed")
class MakeFigureIntegrationTests(unittest.TestCase):
    def test_make_figure_with_synthetic_inputs(self):
        import numpy as np

        from bearing_hic_plot import make_figure

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
            out_path = Path(tmpdir) / "figure.png"

            with patch(
                "bearing_hic_plot.load_contact_matrix",
                side_effect=[(hic_a, 1000), (hic_b, 1000)],
            ), patch(
                "bearing_hic_plot.load_qcat_scores",
                side_effect=[
                    (positions, scores_a, 2),
                    (positions, scores_b, 2),
                ],
            ):
                make_figure(
                    hic_a_path="a.cool",
                    hic_b_path="b.cool",
                    qcat_a_path="a.qcat.bgz",
                    qcat_b_path="b.qcat.bgz",
                    region_str="chr1:0-4000",
                    resolution=1000,
                    out_path=str(out_path),
                    label_a="A",
                    label_b="B",
                    categories=[("State1", "#00aa55"), ("State2", "#d33636")],
                    rgb_hic=False,
                )

            self.assertTrue(out_path.exists())
            self.assertGreater(out_path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
