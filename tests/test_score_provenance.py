import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from score_provenance import score_provenance_signature
import assert_score_provenance as asp


def sig(**kw):
    defaults = dict(
        normalize_tracks=False,
        normalize_method="nonzero-quantile",
        score_method="kl",
        min_signal=0.1,
        categories_path="/repo/categories/mm10_6track_panel.yaml",
        bins_bed_path=None,
        cohort_reference_path=None,
    )
    defaults.update(kw)
    digest, _payload = score_provenance_signature(**defaults)
    return digest


class TestScoreProvenanceSignature(unittest.TestCase):

    def test_deterministic(self):
        self.assertEqual(sig(), sig())

    def test_normalize_on_off_differ(self):
        # The exact failure that bit qnorm: normalize on vs off must differ.
        self.assertNotEqual(sig(normalize_tracks=False),
                            sig(normalize_tracks=True))

    def test_kl_vs_jsd_differ(self):
        self.assertNotEqual(sig(score_method="kl"), sig(score_method="jsd"))

    def test_fixed_vs_adaptive_bins_differ(self):
        self.assertNotEqual(
            sig(bins_bed_path=None),
            sig(bins_bed_path="/repo/results_adaptive/adaptive_segmentation.bed"))

    def test_min_signal_differ(self):
        self.assertNotEqual(sig(min_signal=0.1), sig(min_signal=0.01))

    def test_cohort_reference_differs(self):
        self.assertNotEqual(
            sig(normalize_tracks=True, normalize_method="cohort-quantile",
                cohort_reference_path=None),
            sig(normalize_tracks=True, normalize_method="cohort-quantile",
                cohort_reference_path="/repo/results/cohort_quantile_reference.npz"))

    def test_method_ignored_when_normalize_off(self):
        # With normalization off, the method default must not leak into the sig.
        self.assertEqual(
            sig(normalize_tracks=False, normalize_method="nonzero-quantile"),
            sig(normalize_tracks=False, normalize_method="quantile"))

    def test_path_basename_robust(self):
        # Same basename via different parent dirs -> identical signature, so
        # abs-vs-rel path forms never cause a false mismatch.
        self.assertEqual(
            sig(categories_path="/abs/path/categories/panel.yaml"),
            sig(categories_path="../categories/panel.yaml"))


class TestParseNormFlags(unittest.TestCase):

    def test_empty_is_no_normalize(self):
        ntr, meth, cref = asp.parse_norm_flags("")
        self.assertFalse(ntr)
        self.assertIsNone(cref)

    def test_within_sample_qnorm(self):
        ntr, meth, cref = asp.parse_norm_flags(
            "--normalize-tracks --normalize-method nonzero-quantile")
        self.assertTrue(ntr)
        self.assertEqual(meth, "nonzero-quantile")
        self.assertIsNone(cref)

    def test_cohort(self):
        ntr, meth, cref = asp.parse_norm_flags(
            "--normalize-tracks --normalize-method cohort-quantile "
            "--cohort-reference /repo/results/cohort_quantile_reference.npz")
        self.assertTrue(ntr)
        self.assertEqual(meth, "cohort-quantile")
        self.assertTrue(cref.endswith("cohort_quantile_reference.npz"))


class TestScorerAndGateAgree(unittest.TestCase):
    """The digest the scorer stamps (from its argv) must equal the digest the
    gate expects (parsed from the resolved --norm-flags string). This is the
    contract that makes the gate trustworthy."""

    def _scorer_side(self, normalize_tracks, normalize_method, score_method,
                     min_signal, categories, bins_bed, cohort_reference):
        d, _ = score_provenance_signature(
            normalize_tracks, normalize_method, score_method, min_signal,
            categories, bins_bed, cohort_reference)
        return d

    def _gate_side(self, norm_flags, score_method, min_signal, categories, bins_bed):
        ntr, meth, cref = asp.parse_norm_flags(norm_flags)
        d, _ = score_provenance_signature(
            ntr, meth, score_method, min_signal, categories,
            (bins_bed or None), cref)
        return d

    def test_qnorm_agreement(self):
        cats = "/repo/categories/panel.yaml"
        scorer = self._scorer_side(
            normalize_tracks=True, normalize_method="nonzero-quantile",
            score_method="kl", min_signal=0.1, categories=cats,
            bins_bed=None, cohort_reference=None)
        gate = self._gate_side(
            norm_flags="--normalize-tracks --normalize-method nonzero-quantile",
            score_method="kl", min_signal=0.1, categories=cats, bins_bed="")
        self.assertEqual(scorer, gate)

    def test_production_agreement(self):
        cats = "/repo/categories/panel.yaml"
        # Production: argparse leaves normalize_method at its default even though
        # normalization is off; the gate sees an empty norm-flag string. Both
        # must still agree (method ignored when off).
        scorer = self._scorer_side(
            normalize_tracks=False, normalize_method="nonzero-quantile",
            score_method="kl", min_signal=0.1, categories=cats,
            bins_bed=None, cohort_reference=None)
        gate = self._gate_side(
            norm_flags="", score_method="kl", min_signal=0.1,
            categories=cats, bins_bed="")
        self.assertEqual(scorer, gate)

    def test_stale_reuse_is_caught(self):
        # qnorm config expects normalized; a reused production qcat stamped its
        # un-normalized sig -> the two digests differ -> gate fails. This is the
        # exact protection.
        cats = "/repo/categories/panel.yaml"
        expected_qnorm = self._gate_side(
            "--normalize-tracks --normalize-method nonzero-quantile",
            "kl", 0.1, cats, "")
        reused_production_sig = self._scorer_side(
            normalize_tracks=False, normalize_method="nonzero-quantile",
            score_method="kl", min_signal=0.1, categories=cats,
            bins_bed=None, cohort_reference=None)
        self.assertNotEqual(expected_qnorm, reused_production_sig)


if __name__ == "__main__":
    unittest.main()
