import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bearing.tsv import read_tsv_dict_rows, read_tsv_table


class BearingTsvTests(unittest.TestCase):
    def test_read_tsv_dict_rows_skips_comments(self):
        content = "# comment\nname\tregion\nfoo\tchr1:1-2\n"
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "regions.tsv"
            p.write_text(content)
            rows = read_tsv_dict_rows(p, required_columns=["name", "region"])
        self.assertEqual(rows, [{"name": "foo", "region": "chr1:1-2"}])

    def test_read_tsv_table_returns_fields(self):
        content = "sample\tout\nS1\ta.qcat.bgz\n"
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "samples.tsv"
            p.write_text(content)
            fields, rows = read_tsv_table(p)
        self.assertEqual(fields, ["sample", "out"])
        self.assertEqual(rows[0]["sample"], "S1")

    def test_read_tsv_dict_rows_requires_columns(self):
        content = "sample\tout\nS1\ta.qcat.bgz\n"
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "samples.tsv"
            p.write_text(content)
            with self.assertRaises(ValueError):
                read_tsv_dict_rows(p, required_columns=["sample", "condition"])

    def test_read_tsv_table_normalizes_list_values(self):
        fake_row = {"name": ["wide", "ignored"], "region": ["chr1:1-2"], "label": None}

        class FakeReader:
            fieldnames = ["name", "region", "label"]

            def __iter__(self):
                return iter([fake_row])

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "regions.tsv"
            p.write_text("name\tregion\tlabel\nwide\tchr1:1-2\t\n")
            with patch("bearing.tsv.csv.DictReader", return_value=FakeReader()):
                fields, rows = read_tsv_table(p)

        self.assertEqual(fields, ["name", "region", "label"])
        self.assertEqual(rows[0]["name"], "wide,ignored")
        self.assertEqual(rows[0]["region"], "chr1:1-2")
        self.assertEqual(rows[0]["label"], "")


if __name__ == "__main__":
    unittest.main()
