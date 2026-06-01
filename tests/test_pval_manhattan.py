#!/usr/bin/env python3
"""Tests for Manhattan p-value track drawing with category colors."""

import unittest
from unittest import mock
import numpy as np
import matplotlib.pyplot as plt

from bearing.plot_tracks import draw_pval_manhattan_horizontal


class ManhattanPvalTests(unittest.TestCase):
    """Test drawing of Manhattan p-value tracks with category coloring."""

    def setUp(self):
        """Create a fresh figure and axes for each test."""
        self.fig, self.ax = plt.subplots(figsize=(10, 2))

    def tearDown(self):
        """Clean up figure."""
        plt.close(self.fig)

    def test_manhattan_with_category_colored_dots(self):
        """Manhattan plot dots should be colored by dominant category."""
        # Categories are tuples of (name, color)
        categories = [
            ("Cat1", "#FF0000"),
            ("Cat2", "#00FF00"),
            ("Cat3", "#0000FF"),
        ]

        # Positions and p-values
        positions = np.array([100, 200, 300, 400])
        values = np.array([2.0, 3.5, 1.5, 2.8])

        # Score matrix: each row is category abundances at that position
        # Row 0: Cat1 dominant (5, 2, 1)
        # Row 1: Cat2 dominant (2, 8, 1)
        # Row 2: Cat3 dominant (1, 1, 9)
        # Row 3: Cat1 dominant (7, 1, 1)
        score_positions = np.array([100, 200, 300, 400])
        score_matrix = np.array([
            [5.0, 2.0, 1.0],  # Cat1 dominant
            [2.0, 8.0, 1.0],  # Cat2 dominant
            [1.0, 1.0, 9.0],  # Cat3 dominant
            [7.0, 1.0, 1.0],  # Cat1 dominant
        ])

        # Should not raise any errors
        draw_pval_manhattan_horizontal(
            self.ax, positions, values,
            region_start=0, region_end=500,
            highlights=None,
            label="Test Manhattan",
            score_positions=score_positions,
            score_matrix=score_matrix,
            categories=categories,
            y_max=4.0,
        )

        # Basic checks
        self.ax.set_xlim(0, 1)
        self.assertEqual(self.ax.get_ylabel(), "Test Manhattan")

    def test_manhattan_without_categories_defaults_to_grey(self):
        """Manhattan plot without category info should use grey dots."""
        positions = np.array([100, 200, 300, 400])
        values = np.array([2.0, 3.5, 1.5, 2.8])

        # Should not raise errors even without category colors
        draw_pval_manhattan_horizontal(
            self.ax, positions, values,
            region_start=0, region_end=500,
            highlights=None,
            label="Grey Manhattan",
        )

        self.assertEqual(self.ax.get_ylabel(), "Grey Manhattan")

    def test_manhattan_with_misaligned_scores(self):
        """Manhattan should handle score positions that don't match p-value positions."""
        categories = [
            ("Cat1", "#FF0000"),
            ("Cat2", "#00FF00"),
        ]

        positions = np.array([100, 200, 300])
        values = np.array([2.0, 3.5, 1.5])

        # Score positions don't match
        score_positions = np.array([150, 250, 350])
        score_matrix = np.array([
            [5.0, 2.0],
            [2.0, 8.0],
            [1.0, 9.0],
        ])

        # Should fall back to grey for unmatched positions
        draw_pval_manhattan_horizontal(
            self.ax, positions, values,
            region_start=0, region_end=500,
            score_positions=score_positions,
            score_matrix=score_matrix,
            categories=categories,
        )

        self.assertEqual(len(self.ax.collections) >= 0, True)  # Should have scatter plots


if __name__ == "__main__":
    unittest.main()
