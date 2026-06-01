import importlib.util
import unittest

import numpy as np


HAS_MPL = importlib.util.find_spec("matplotlib") is not None
HAS_NUMPY = importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(HAS_MPL and HAS_NUMPY, "matplotlib/numpy not installed")
class PvalDiffLollipopOverlayTests(unittest.TestCase):
    def _path_collection_hexes(self, ax):
        from matplotlib.collections import PathCollection
        from matplotlib.colors import to_hex

        hexes = []
        for coll in ax.collections:
            if isinstance(coll, PathCollection) and len(coll.get_facecolors()) > 0:
                hexes.append(to_hex(coll.get_facecolors()[0]))
        return hexes

    def test_horizontal_lollipops_are_colored_by_dominant_track(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from bearing_hic_plot import draw_pval_diff_horizontal

        fig, ax = plt.subplots(figsize=(4, 2))
        positions = np.array([1000, 2000, 3000], dtype=float)
        values = np.array([5.0, -6.0, 0.5], dtype=float)
        diff_positions = positions + 100
        diff_scores = np.array([
            [10.0, 1.0],
            [-1.0, -10.0],
            [8.0, 2.0],
        ], dtype=float)
        categories = [("RNAseq+", "#6495ed"), ("CTCF", "#ff2200")]

        draw_pval_diff_horizontal(
            ax,
            positions,
            values,
            0,
            4000,
            cutoff_value=2.0,
            diff_score_positions=diff_positions,
            diff_score_matrix=diff_scores,
            categories=categories,
        )

        hexes = self._path_collection_hexes(ax)
        self.assertEqual(len(hexes), 2)
        self.assertEqual(set(hexes), {"#6495ed", "#ff2200"})
        plt.close(fig)

    def test_horizontal_lollipops_fall_back_to_grey_when_unaligned(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from bearing_hic_plot import draw_pval_diff_horizontal

        fig, ax = plt.subplots(figsize=(4, 2))
        positions = np.array([1000], dtype=float)
        values = np.array([5.0], dtype=float)
        diff_positions = np.array([15000], dtype=float)
        diff_scores = np.array([[10.0, 1.0]], dtype=float)
        categories = [("RNAseq+", "#6495ed"), ("CTCF", "#ff2200")]

        draw_pval_diff_horizontal(
            ax,
            positions,
            values,
            0,
            4000,
            cutoff_value=2.0,
            diff_score_positions=diff_positions,
            diff_score_matrix=diff_scores,
            categories=categories,
        )

        hexes = self._path_collection_hexes(ax)
        self.assertEqual(hexes, ["#9aa0a6"])
        plt.close(fig)

    def test_horizontal_without_diff_inputs_has_no_lollipop_overlay(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from bearing_hic_plot import draw_pval_diff_horizontal

        fig, ax = plt.subplots(figsize=(4, 2))
        positions = np.array([1000], dtype=float)
        values = np.array([5.0], dtype=float)

        draw_pval_diff_horizontal(
            ax,
            positions,
            values,
            0,
            4000,
            cutoff_value=2.0,
        )

        self.assertEqual(self._path_collection_hexes(ax), [])
        plt.close(fig)

    def test_vertical_lollipops_are_colored_by_dominant_track(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from bearing_hic_plot import draw_pval_diff_vertical

        fig, ax = plt.subplots(figsize=(4, 2))
        positions = np.array([1000, 2000], dtype=float)
        values = np.array([5.0, -6.0], dtype=float)
        diff_positions = positions + 100
        diff_scores = np.array([
            [10.0, 1.0],
            [-1.0, -10.0],
        ], dtype=float)
        categories = [("RNAseq+", "#6495ed"), ("CTCF", "#ff2200")]

        draw_pval_diff_vertical(
            ax,
            positions,
            values,
            0,
            4000,
            cutoff_value=2.0,
            diff_score_positions=diff_positions,
            diff_score_matrix=diff_scores,
            categories=categories,
        )

        hexes = self._path_collection_hexes(ax)
        self.assertEqual(set(hexes), {"#6495ed", "#ff2200"})
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
