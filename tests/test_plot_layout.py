import importlib.util
import unittest


HAS_MPL = importlib.util.find_spec("matplotlib") is not None


@unittest.skipUnless(HAS_MPL, "matplotlib not installed")
class PlotLayoutTests(unittest.TestCase):
    def test_create_main_figure_axes_keys(self):
        from bearing.plot_layout import create_main_figure

        fig, axes = create_main_figure(use_pval_overlay=False)
        try:
            expected = {
                "hic", "epi_a", "pval_a", "loop_a", "diff_a", "diff_pval_a",
                "gene_h", "axis_h",
                "epi_b", "pval_b", "loop_b", "diff_b", "diff_pval_b",
                "gene_v", "axis_v", "legend",
            }
            self.assertEqual(set(axes.keys()), expected)
            self.assertEqual(len(fig.axes), len(expected))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
