import unittest

from bearing.sheet import extract_bw_values, sample_name_from_row


class BearingSheetTests(unittest.TestCase):
    def test_sample_name_prefers_sample(self):
        self.assertEqual(sample_name_from_row({"sample": "S1", "name": "N1"}), "S1")

    def test_sample_name_falls_back_to_name(self):
        self.assertEqual(sample_name_from_row({"name": "N1"}), "N1")

    def test_extract_bw_values_from_bw_column(self):
        row = {"bw": "a.bw, b.bw ,,c.bw"}
        self.assertEqual(extract_bw_values(row), ["a.bw", "b.bw", "c.bw"])

    def test_extract_bw_values_from_indexed_columns(self):
        row = {"bw1": "a.bw", "bw2": "", "bw3": "c.bw"}
        self.assertEqual(extract_bw_values(row), ["a.bw", "c.bw"])


if __name__ == "__main__":
    unittest.main()
