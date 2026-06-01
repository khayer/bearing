import tempfile
import unittest
import importlib.util
from pathlib import Path
from unittest.mock import patch

import numpy as np

import compare_qcat as cq


HAS_SCIPY = importlib.util.find_spec("scipy") is not None


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class CompareQcatNonzeroModeTests(unittest.TestCase):
    def test_primary_spearman_uses_global_nonzero_shared_set(self):
        # 4 informative bins + 1 jointly-zero bin (dropped by the union mask).
        bins = [
            ("chr1", 0, 100),
            ("chr1", 100, 200),
            ("chr1", 200, 300),
            ("chr1", 300, 400),
            ("chr1", 400, 500),
        ]
        mats = {
            "A_rep1.qcat.bgz": np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]]),
            "A_rep2.qcat.bgz": np.array([[1.1, 0.0], [1.9, 0.1], [0.0, 0.9], [0.9, 1.1], [0.0, 0.0]]),
            "B_rep1.qcat.bgz": np.array([[0.0, 1.0], [0.0, 2.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]]),
            "B_rep2.qcat.bgz": np.array([[0.1, 0.9], [0.0, 2.1], [0.9, 0.1], [1.0, 1.0], [0.0, 0.0]]),
        }
        samples = [
            {"sample": "A_rep1", "condition": "A", "replicate": 1, "qcat": Path("A_rep1.qcat.bgz")},
            {"sample": "A_rep2", "condition": "A", "replicate": 2, "qcat": Path("A_rep2.qcat.bgz")},
            {"sample": "B_rep1", "condition": "B", "replicate": 1, "qcat": Path("B_rep1.qcat.bgz")},
            {"sample": "B_rep2", "condition": "B", "replicate": 2, "qcat": Path("B_rep2.qcat.bgz")},
        ]

        spearman_calls = []
        spearman_plot_calls = []

        def fake_parse_qcat_bgz(path, chroms=None, include_raw=False):
            mat = mats[Path(path).name]
            return bins, mat, mat.shape[1], {}

        def fake_build_spearman_matrix(score_mats, sample_names, num_states,
                                       nonzero_mode="all", min_bins=2,
                                       return_counts=False):
            spearman_calls.append({
                "nonzero_mode": nonzero_mode,
                "return_counts": bool(return_counts),
                "n_rows": int(score_mats[0].shape[0]),
            })
            n = len(sample_names)
            rho_per_state = np.full((n, n, num_states), 0.25, dtype=np.float64)
            for si in range(num_states):
                np.fill_diagonal(rho_per_state[:, :, si], 1.0)
            rho_mean = np.full((n, n), 0.25, dtype=np.float64)
            np.fill_diagonal(rho_mean, 1.0)
            if return_counts:
                bins_total = np.full((n, n, num_states), score_mats[0].shape[0], dtype=np.int32)
                return rho_per_state, rho_mean, bins_total.copy(), bins_total
            return rho_per_state, rho_mean

        def fake_plot_spearman(rho_mean, rho_per_state, sample_names, conditions,
                               num_states, out_path, categories=None, title=None):
            spearman_plot_calls.append((Path(out_path).name, title))

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            with patch.object(cq, "load_sample_sheet", return_value=samples), \
                 patch.object(cq, "parse_qcat_bgz", side_effect=fake_parse_qcat_bgz), \
                 patch.object(cq, "plot_jsd_heatmap"), \
                 patch.object(cq, "plot_spearman", side_effect=fake_plot_spearman), \
                 patch.object(cq, "plot_total_saliency_spearman"), \
                 patch.object(cq, "plot_pca"), \
                 patch.object(cq, "write_compare_ini"), \
                 patch.object(cq, "plot_qcat_region"), \
                 patch.object(cq, "build_spearman_matrix", side_effect=fake_build_spearman_matrix), \
                 patch.object(cq, "build_jsd_matrix", return_value=np.zeros((4, 4), dtype=np.float64)), \
                 patch.object(cq, "write_spearman_nonzero_diagnostics"):
                cq.run(
                    sheet_path="unused.tsv",
                    out_dir=out_dir,
                    skip_diff=True,
                    skip_pca=True,
                    clip_for_ini=False,
                    spearman_nonzero_only=True,
                    spearman_min_bins=2,
                    global_nonzero_mode="any",
                    global_nonzero_min_bins=2,
                    skip_q_pair_jsd=True,
                )

            # Primary reported Spearman: all bins of the global-nonzero shared
            # set (the jointly-zero bin is dropped, so 4 of 5 bins).
            primary = [c for c in spearman_calls
                       if not c["return_counts"] and c["nonzero_mode"] == "all"
                       and c["n_rows"] == 4]
            # Supplementary: Spearman over ALL aligned bins (5 of 5).
            supp_all = [c for c in spearman_calls
                        if not c["return_counts"] and c["nonzero_mode"] == "all"
                        and c["n_rows"] == 5]
            self.assertTrue(primary,
                            "primary Spearman did not use the global-nonzero shared set")
            self.assertTrue(supp_all,
                            "supplementary all-aligned-bins Spearman missing")

            self.assertTrue((out_dir / "nonzero_bin_filter_summary.tsv").exists())
            self.assertTrue((out_dir / "total_saliency_spearman.tsv").exists())
            self.assertTrue((out_dir / "total_saliency_spearman_all_aligned_bins.tsv").exists())

            primary_titles = [t for name, t in spearman_plot_calls if name == "spearman.pdf"]
            self.assertTrue(any("global nonzero" in (t or "") for t in primary_titles))


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class TotalSaliencySpearmanTests(unittest.TestCase):
    def test_total_saliency_spearman_shape_and_group_separation(self):
        names = ["A_rep1", "A_rep2", "B_rep1", "B_rep2"]

        base_a = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], dtype=np.float64)
        base_b = base_a[::-1]

        jitter_a = np.array([0.0, 0.1, 0.0, 0.1, 0.0, -0.1, 0.0, -0.1, 0.0, 0.1, 0.0, 0.1])
        jitter_b = np.array([0.1, 0.0, 0.1, 0.0, -0.1, 0.0, -0.1, 0.0, 0.1, 0.0, 0.1, 0.0])

        totals = {
            "A_rep1": base_a,
            "A_rep2": np.clip(base_a + jitter_a, 0.0, None),
            "B_rep1": base_b,
            "B_rep2": np.clip(base_b + jitter_b, 0.0, None),
        }

        weights = {
            "A_rep1": np.array([0.55, 0.30, 0.15]),
            "A_rep2": np.array([0.52, 0.31, 0.17]),
            "B_rep1": np.array([0.20, 0.30, 0.50]),
            "B_rep2": np.array([0.18, 0.34, 0.48]),
        }

        mats = []
        for name in names:
            t = totals[name][:, None]
            w = weights[name][None, :]
            mats.append(t * w)

        rho = cq.build_total_saliency_spearman(
            mats,
            names,
            nonzero_mode="any",
            min_bins=2,
        )

        self.assertEqual(rho.shape, (4, 4))
        self.assertTrue(np.allclose(np.diag(rho), 1.0))

        within = np.array([rho[0, 1], rho[2, 3]], dtype=np.float64)
        cross = np.array([rho[0, 2], rho[0, 3], rho[1, 2], rho[1, 3]], dtype=np.float64)

        self.assertGreater(float(np.mean(within)), float(np.mean(cross)))


if __name__ == "__main__":
    unittest.main()
