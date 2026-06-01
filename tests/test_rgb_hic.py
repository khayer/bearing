import importlib.util
import unittest


HAS_MPL = importlib.util.find_spec("matplotlib") is not None
HAS_NUMPY = importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(HAS_MPL and HAS_NUMPY, "matplotlib/numpy not installed")
class RgbHicTests(unittest.TestCase):
    def test_magenta_green_default_mapping(self):
        import numpy as np

        from bearing_hic_plot import make_rgb_hic

        mat_a = np.array([[10.0, 0.0], [10.0, 10.0]], dtype=np.float64)
        mat_b = np.array([[0.0, 10.0], [10.0, 0.0]], dtype=np.float64)

        rgb = make_rgb_hic(mat_a, mat_b)

        # A-only should be magenta.
        self.assertTrue(np.array_equal(rgb[0, 0], np.array([255, 0, 255], dtype=np.uint8)))
        # B-only should be green.
        self.assertTrue(np.array_equal(rgb[0, 1], np.array([0, 255, 0], dtype=np.uint8)))
        # Co-enriched should be white.
        self.assertTrue(np.array_equal(rgb[1, 0], np.array([255, 255, 255], dtype=np.uint8)))

    def test_red_green_palette_mapping(self):
        import numpy as np

        from bearing_hic_plot import make_rgb_hic

        mat_a = np.array([[10.0, 0.0], [10.0, 10.0]], dtype=np.float64)
        mat_b = np.array([[0.0, 10.0], [10.0, 0.0]], dtype=np.float64)

        rgb = make_rgb_hic(mat_a, mat_b, palette="red-green")

        # A-only should be red.
        self.assertTrue(np.array_equal(rgb[0, 0], np.array([255, 0, 0], dtype=np.uint8)))
        # B-only should be green.
        self.assertTrue(np.array_equal(rgb[0, 1], np.array([0, 255, 0], dtype=np.uint8)))
        # Co-enriched should be yellow.
        self.assertTrue(np.array_equal(rgb[1, 0], np.array([255, 255, 0], dtype=np.uint8)))


if __name__ == "__main__":
    unittest.main()
