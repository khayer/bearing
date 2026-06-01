import importlib.util
import unittest


HAS_MPL = importlib.util.find_spec("matplotlib") is not None
HAS_NP = importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(HAS_MPL and HAS_NP, "matplotlib/numpy not installed")
class PlotLegendTests(unittest.TestCase):
    def test_rgb_triangle_points_downward(self):
        import matplotlib.pyplot as plt
        from bearing.plot_legend import draw_legend_panel

        fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
        try:
            categories = [("ATAC", "#00b050"), ("CTCF", "#ff2200")]
            draw_legend_panel(
                ax,
                categories,
                num_states=2,
                hic_cmap=None,
                rgb_mode=True,
                label_a="A",
                label_b="B",
                rgb_palette="magenta-green",
            )
            # In RGB mode, an inset axis with a triangular image is added.
            inset_axes = [a for a in fig.axes if a is not ax]
            self.assertTrue(inset_axes)
            grad_ax = inset_axes[0]
            img = grad_ax.images[0].get_array()
            top = img[int(img.shape[0] * 0.10), img.shape[1] // 2]
            bottom = img[int(img.shape[0] * 0.75), img.shape[1] // 2]
            self.assertLess(int(top[3]), 255)
            self.assertGreater(int(bottom[3]), 0)
            self.assertGreater(int(bottom[3]), int(top[3]))
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
