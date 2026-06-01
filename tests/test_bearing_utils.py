import unittest
from pathlib import Path

from bearing.utils import (
    parse_key_value_items,
    parse_named_region_item,
    parse_ucsc_region,
    resolve_path,
    sanitize_token,
)


class BearingUtilsTests(unittest.TestCase):
    def test_parse_key_value_items(self):
        got = parse_key_value_items(["A=one", "B=two"], "--map")
        self.assertEqual(got, {"A": "one", "B": "two"})

    def test_parse_key_value_items_rejects_invalid(self):
        with self.assertRaises(ValueError):
            parse_key_value_items(["A"], "--map")

    def test_sanitize_token(self):
        self.assertEqual(sanitize_token("A B/C"), "A_B_C")

    def test_resolve_path(self):
        base = Path("/tmp/example")
        resolved = resolve_path("sub/file.txt", base)
        self.assertEqual(resolved, (base / "sub/file.txt").resolve())

    def test_parse_ucsc_region(self):
        self.assertEqual(parse_ucsc_region("chr6:41,000,000-41,100,000"), ("chr6", 41000000, 41100000))

    def test_parse_named_region_item(self):
        self.assertEqual(
            parse_named_region_item("vdj=chr6:40793981-41688054"),
            {"name": "vdj", "region": "chr6:40793981-41688054"},
        )


if __name__ == "__main__":
    unittest.main()
