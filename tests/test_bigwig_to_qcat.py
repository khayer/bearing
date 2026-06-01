import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np

import bigwig_to_qcat as bq


class BlacklistHelpersTests(unittest.TestCase):
    def test_load_blacklist_parses_and_sorts(self):
        with tempfile.TemporaryDirectory() as td:
            bed = Path(td) / "blacklist.bed"
            bed.write_text(
                "# comment\n"
                "track name=test\n"
                "chr1\t200\t300\tfoo\n"
                "chr1\t100\t150\n"
                "chr2\t50\t80\n"
            )
            out = bq.load_blacklist(str(bed))

        self.assertEqual(out["chr1"], [(100, 150), (200, 300)])
        self.assertEqual(out["chr2"], [(50, 80)])

    def test_bins_overlapping_blacklist_boundaries(self):
        bins = [(0, 200), (200, 400), (400, 600), (600, 800)]
        blacklist = {"chr1": [(150, 220), (590, 610)]}
        mask = bq.bins_overlapping_blacklist(bins, blacklist, "chr1")
        self.assertTrue(np.array_equal(mask, np.array([True, True, True, True], dtype=bool)))

    def test_write_blacklist_bed(self):
        flagged = {
            "chr1": np.array([False, True, True, False], dtype=bool),
            "chr2": np.array([True], dtype=bool),
        }
        chrom_sizes = {"chr1": 1000, "chr2": 150}

        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "detected.bed"
            bq.write_blacklist_bed(flagged, chrom_sizes, str(out), bin_size=200)
            lines = [ln.strip() for ln in out.read_text().splitlines() if ln and not ln.startswith("#")]

        self.assertEqual(lines, ["chr1\t200\t400", "chr1\t400\t600", "chr2\t0\t150"])


class DetectUnmappableTests(unittest.TestCase):
    def test_detect_unmappable_bins_threshold_and_min_samples(self):
        chrom_sizes = {"chr1": 400}
        sheet_jobs = [
            {"sample": "s1", "bw_paths": [Path("s1_a.bw"), Path("s1_b.bw")]},
            {"sample": "s2", "bw_paths": [Path("s2_a.bw"), Path("s2_b.bw")]},
        ]

        values = {
            "s1_a.bw": [0.0, 0.5],
            "s1_b.bw": [0.0, 0.0],
            "s2_a.bw": [0.0, 0.7],
            "s2_b.bw": [0.0, 0.0],
        }

        class FakeBW:
            def __init__(self, path):
                self.path = path

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def close(self):
                pass

            def stats(self, chrom, start, end, type, nBins):
                self.assert_equal = (chrom, start, end, type, nBins)
                return values[self.path]

        fake_module = types.SimpleNamespace(open=lambda path: FakeBW(path))
        old_pybigwig = sys.modules.get("pyBigWig")
        sys.modules["pyBigWig"] = fake_module
        try:
            flagged = bq.detect_unmappable_bins(
                sheet_jobs,
                chrom_sizes,
                zero_frac=0.9,
                min_samples=2,
                bin_size=200,
                signal_threshold=0.01,
            )
        finally:
            if old_pybigwig is None:
                del sys.modules["pyBigWig"]
            else:
                sys.modules["pyBigWig"] = old_pybigwig

        self.assertIn("chr1", flagged)
        # Bin0 is zero across all combos -> flagged
        # Bin1 has only 2/4 zero combos and only one sample fully zero -> not flagged
        self.assertTrue(np.array_equal(flagged["chr1"], np.array([True, False], dtype=bool)))

    def test_detect_unmappable_bins_aborts_when_bigwig_unreadable(self):
        chrom_sizes = {"chr1": 400}
        sheet_jobs = [{"sample": "s1", "bw_paths": [Path("missing.bw")]}]

        old_checker = bq.check_bigwig_paths
        bq.check_bigwig_paths = lambda paths: [("missing.bw", "open returned None")]
        try:
            with self.assertRaises(SystemExit) as exc:
                bq.detect_unmappable_bins(
                    sheet_jobs,
                    chrom_sizes,
                    zero_frac=0.9,
                    min_samples=1,
                    bin_size=200,
                )
        finally:
            bq.check_bigwig_paths = old_checker

        self.assertIn("missing/unreadable BigWig", str(exc.exception))


if __name__ == "__main__":
    unittest.main()
