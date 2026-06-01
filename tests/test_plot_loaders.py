import unittest

from bearing.plot_loaders import cool_resolution_variant, detect_hic_format


class PlotLoadersTests(unittest.TestCase):
    def test_detect_hic_format(self):
        self.assertEqual(detect_hic_format("a.cool"), "cool")
        self.assertEqual(detect_hic_format("a.mcool"), "cool")
        self.assertEqual(detect_hic_format("a.hic"), "hic")

    def test_cool_resolution_variant(self):
        p = cool_resolution_variant("sample_bs_10000.cool", 5000)
        self.assertTrue(str(p).endswith("sample_bs_5000.cool"))


if __name__ == "__main__":
    unittest.main()
