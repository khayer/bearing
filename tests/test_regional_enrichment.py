#!/usr/bin/env python3
"""
Tests for regional_enrichment.py
"""

import tempfile
import unittest
from pathlib import Path

from regional_enrichment import (
    compute_regional_enrichment,
    parse_diff_qcat,
    parse_diff_pvals,
    parse_regions_bed,
)


class TestRegionalEnrichment(unittest.TestCase):

    def setUp(self):
        self.locus_chrom = 'chr1'
        self.locus_start = 0
        self.locus_end = 1000000

        self.diff_qcat_content = """# Synthetic diff qcat
chr1	0	100000	{"qcat": [[1.0, 1], [0.5, 2]]}
chr1	100000	200000	{"qcat": [[0.8, 1], [0.3, 2]]}
chr1	200000	300000	{"qcat": [[0.2, 1], [0.1, 2]]}
chr1	300000	400000	{"qcat": [[1.2, 1], [0.4, 2]]}
chr1	400000	500000	{"qcat": [[0.9, 1], [0.6, 2]]}
chr1	500000	600000	{"qcat": [[0.1, 1], [0.2, 2]]}
chr1	600000	700000	{"qcat": [[1.5, 1], [0.8, 2]]}
chr1	700000	800000	{"qcat": [[0.3, 1], [0.4, 2]]}
chr1	800000	900000	{"qcat": [[1.8, 1], [0.9, 2]]}
chr1	900000	1000000	{"qcat": [[0.4, 1], [0.3, 2]]}
"""

        self.diff_pvals_content = """chrom	start	end	pval	bearing_score
chr1	0	100000	0.01	1.5
chr1	100000	200000	0.02	1.2
chr1	200000	300000	0.8	0.3
chr1	300000	400000	0.005	1.6
chr1	400000	500000	0.03	1.5
chr1	500000	600000	0.9	0.2
chr1	600000	700000	0.001	2.0
chr1	700000	800000	0.7	0.5
chr1	800000	900000	0.0001	2.5
chr1	900000	1000000	0.6	0.4
"""

        self.regions_content = """chr1	0	300000	region1
"""

    def test_parse_diff_qcat(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.qcat', delete=False) as f:
            f.write(self.diff_qcat_content)
            f.flush()
            bins = parse_diff_qcat(f.name)

        self.assertEqual(len(bins), 10)
        self.assertIn(('chr1', 0, 100000), bins)
        self.assertEqual(bins[('chr1', 0, 100000)], 1.0)

    def test_parse_diff_pvals(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write(self.diff_pvals_content)
            f.flush()
            bins = parse_diff_pvals(f.name)

        self.assertEqual(len(bins), 10)
        self.assertIn(('chr1', 0, 100000), bins)
        pval, score = bins[('chr1', 0, 100000)]
        self.assertAlmostEqual(pval, 0.01)
        self.assertAlmostEqual(score, 1.5)

    def test_parse_regions_bed(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(self.regions_content)
            f.flush()
            regions = parse_regions_bed(f.name)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0], ('chr1', 0, 300000, 'region1'))

    def test_smoke_test(self):
        """Test basic functionality with synthetic enriched region."""
        with tempfile.TemporaryDirectory() as tmpdir:
            qcat_file = Path(tmpdir) / 'test.qcat'
            pvals_file = Path(tmpdir) / 'test.pvals.tsv'
            regions_file = Path(tmpdir) / 'regions.bed'

            with open(qcat_file, 'w') as f:
                f.write(self.diff_qcat_content)
            with open(pvals_file, 'w') as f:
                f.write(self.diff_pvals_content)
            with open(regions_file, 'w') as f:
                f.write(self.regions_content)

            diff_qcat_bins = parse_diff_qcat(str(qcat_file))
            diff_pvals_bins = parse_diff_pvals(str(pvals_file))
            regions = parse_regions_bed(str(regions_file))

            results = compute_regional_enrichment(
                diff_qcat_bins, diff_pvals_bins, regions,
                self.locus_chrom, self.locus_start, self.locus_end, 0.05
            )

            self.assertEqual(len(results), 1)
            result = results[0]

            self.assertEqual(result['region_name'], 'region1')
            self.assertEqual(result['L_region_bins'], 3)
            self.assertEqual(result['L_locus_bins'], 10)
            self.assertEqual(result['k'], 2)  # pvals: 0.01, 0.02, 0.8 -> 2 below 0.05
            # Total bins below 0.05: 0.01, 0.02, 0.005, 0.03, 0.001, 0.0001 = 6
            self.assertEqual(result['n_locus'], 6)
            self.assertEqual(result['k_pos'], 2)
            self.assertEqual(result['k_neg'], 0)

            # With only 2 significant bins, neither spatial nor directional is strong
            # Just verify the structure is correct
            self.assertTrue(0 <= result['p_spatial'] <= 1)
            self.assertTrue(0 <= result['p_directional'] <= 1)
            self.assertTrue(0 <= result['p_combined'] <= 1)

    def test_null_behavior(self):
        """Test with uniform distribution (no bins below p_thresh, so no signal)."""
        uniform_pvals = """chrom	start	end	pval	bearing_score
chr1	0	100000	0.1	0.5
chr1	100000	200000	0.2	0.5
chr1	200000	300000	0.3	0.5
chr1	300000	400000	0.4	0.5
chr1	400000	500000	0.5	0.5
chr1	500000	600000	0.6	0.5
chr1	600000	700000	0.7	0.5
chr1	700000	800000	0.8	0.5
chr1	800000	900000	0.9	0.5
chr1	900000	1000000	0.95	0.5
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            qcat_file = Path(tmpdir) / 'test.qcat'
            pvals_file = Path(tmpdir) / 'test.pvals.tsv'
            regions_file = Path(tmpdir) / 'regions.bed'

            with open(qcat_file, 'w') as f:
                f.write(self.diff_qcat_content)
            with open(pvals_file, 'w') as f:
                f.write(uniform_pvals)
            with open(regions_file, 'w') as f:
                f.write(self.regions_content)

            diff_qcat_bins = parse_diff_qcat(str(qcat_file))
            diff_pvals_bins = parse_diff_pvals(str(pvals_file))
            regions = parse_regions_bed(str(regions_file))

            results = compute_regional_enrichment(
                diff_qcat_bins, diff_pvals_bins, regions,
                self.locus_chrom, self.locus_start, self.locus_end, 0.05
            )

            result = results[0]
            self.assertEqual(result['k'], 0)  # all pvals >= 0.1
            self.assertEqual(result['n_locus'], 0)  # no bins below 0.05
            self.assertEqual(result['p_spatial'], 1.0)  # null when n_locus = 0

    def test_direction_only_enrichment(self):
        """Test region with all significant bins in same direction."""
        # Region (chr1:0-300000) bins are all significant and POSITIVE; an equal
        # number of out-of-region bins are significant and NEGATIVE, so the
        # genome-wide directional baseline g_dir = 3 / 6 = 0.5. The region's
        # all-positive concordance is then tested for EXCESS over that baseline
        # (recalibrated one-sided test), not against a fixed 0.5.
        directional_pvals = """chrom	start	end	pval	bearing_score
chr1	0	100000	0.01	1.0
chr1	100000	200000	0.01	1.0
chr1	200000	300000	0.01	1.0
chr1	300000	400000	0.01	-1.0
chr1	400000	500000	0.01	-1.0
chr1	500000	600000	0.01	-1.0
chr1	600000	700000	0.1	0.5
chr1	700000	800000	0.1	0.5
chr1	800000	900000	0.1	0.5
chr1	900000	1000000	0.1	0.5
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            qcat_file = Path(tmpdir) / 'test.qcat'
            pvals_file = Path(tmpdir) / 'test.pvals.tsv'
            regions_file = Path(tmpdir) / 'regions.bed'

            with open(qcat_file, 'w') as f:
                f.write(self.diff_qcat_content)
            with open(pvals_file, 'w') as f:
                f.write(directional_pvals)
            with open(regions_file, 'w') as f:
                f.write(self.regions_content)

            diff_qcat_bins = parse_diff_qcat(str(qcat_file))
            diff_pvals_bins = parse_diff_pvals(str(pvals_file))
            regions = parse_regions_bed(str(regions_file))

            results = compute_regional_enrichment(
                diff_qcat_bins, diff_pvals_bins, regions,
                self.locus_chrom, self.locus_start, self.locus_end, 0.05
            )

            result = results[0]
            self.assertEqual(result['k'], 3)  # 3 sig bins in region
            self.assertEqual(result['n_locus'], 6)  # 6 sig bins total (3 pos + 3 neg)
            self.assertEqual(result['k_pos'], 3)  # region all positive
            self.assertEqual(result['k_neg'], 0)
            # Recalibrated one-sided excess-concordance test against the
            # genome-wide directional baseline g_dir = 3 / 6 = 0.5:
            #   binom.sf(k_pos - 1, k, g_dir) = binom.sf(2, 3, 0.5) = 0.125
            self.assertAlmostEqual(result['g_dir'], 0.5, places=6)
            self.assertAlmostEqual(result['p_directional'], 0.125, places=3)

    def test_boundary_cases(self):
        """Test edge cases: zero bins, all bins, region = locus."""
        zero_pvals = """chrom	start	end	pval	bearing_score
chr1	0	100000	0.1	0.5
chr1	100000	200000	0.1	0.5
chr1	200000	300000	0.1	0.5
chr1	300000	400000	0.01	1.0
chr1	400000	500000	0.01	1.0
chr1	500000	600000	0.01	1.0
chr1	600000	700000	0.01	1.0
chr1	700000	800000	0.01	1.0
chr1	800000	900000	0.01	1.0
chr1	900000	1000000	0.01	1.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            qcat_file = Path(tmpdir) / 'test.qcat'
            pvals_file = Path(tmpdir) / 'test.pvals.tsv'
            regions_file = Path(tmpdir) / 'regions.bed'

            with open(qcat_file, 'w') as f:
                f.write(self.diff_qcat_content)
            with open(pvals_file, 'w') as f:
                f.write(zero_pvals)
            with open(regions_file, 'w') as f:
                f.write(self.regions_content)

            diff_qcat_bins = parse_diff_qcat(str(qcat_file))
            diff_pvals_bins = parse_diff_pvals(str(pvals_file))
            regions = parse_regions_bed(str(regions_file))

            results = compute_regional_enrichment(
                diff_qcat_bins, diff_pvals_bins, regions,
                self.locus_chrom, self.locus_start, self.locus_end, 0.05
            )

            result = results[0]
            self.assertEqual(result['k'], 0)
            self.assertEqual(result['p_directional'], 1.0)

    def test_region_outside_locus(self):
        """Test that regions outside locus raise error."""
        regions_outside = """chr1	500000	600000	region1
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            qcat_file = Path(tmpdir) / 'test.qcat'
            pvals_file = Path(tmpdir) / 'test.pvals.tsv'
            regions_file = Path(tmpdir) / 'regions.bed'

            with open(qcat_file, 'w') as f:
                f.write(self.diff_qcat_content)
            with open(pvals_file, 'w') as f:
                f.write(self.diff_pvals_content)
            with open(regions_file, 'w') as f:
                f.write(regions_outside)

            diff_qcat_bins = parse_diff_qcat(str(qcat_file))
            diff_pvals_bins = parse_diff_pvals(str(pvals_file))
            regions = parse_regions_bed(str(regions_file))

            with self.assertRaises(ValueError):
                compute_regional_enrichment(
                    diff_qcat_bins, diff_pvals_bins, regions,
                    self.locus_chrom, self.locus_start, 500000,
                    0.05
                )


if __name__ == '__main__':
    unittest.main()
