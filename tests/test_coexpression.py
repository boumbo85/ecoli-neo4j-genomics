"""Unit tests for the co-expression cohesion functions."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import coexpression  # noqa: E402


# A tiny symmetric network:
#   A-B (score 100), B-C (score 50), C-D (score 10), D-E (score 5)
# Genes {A,B,C} form a tight cluster; D,E are loosely attached.
ADJ = {
    "A": {"B": 100},
    "B": {"A": 100, "C": 50},
    "C": {"B": 50, "D": 10},
    "D": {"C": 10, "E": 5},
    "E": {"D": 5},
}
POP = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]  # F-J have no edges


class TestIntraSetStrength(unittest.TestCase):
    def test_counts_internal_edges_once(self):
        # A,B,C: edges A-B (100) + B-C (50) = 150
        self.assertEqual(coexpression.intra_set_strength(ADJ, {"A", "B", "C"}), 150)

    def test_ignores_edges_leaving_the_set(self):
        # Just {A,B}: only A-B = 100 (B-C leaves the set)
        self.assertEqual(coexpression.intra_set_strength(ADJ, {"A", "B"}), 100)

    def test_disconnected_genes_score_zero(self):
        self.assertEqual(coexpression.intra_set_strength(ADJ, {"F", "G"}), 0)


class TestIntraSetEdges(unittest.TestCase):
    def test_edge_count(self):
        self.assertEqual(coexpression.intra_set_edges(ADJ, {"A", "B", "C"}), 2)

    def test_each_edge_once(self):
        self.assertEqual(coexpression.intra_set_edges(ADJ, {"A", "B"}), 1)


class TestPermutationTest(unittest.TestCase):
    def test_tight_cluster_beats_random(self):
        res = coexpression.permutation_test(ADJ, {"A", "B", "C"}, POP,
                                            n_perm=500, seed=1)
        # the real cluster should be well above the random mean
        self.assertGreater(res["observed"], res["mean_random"])
        self.assertLessEqual(res["p_value"], 0.05)

    def test_pvalue_never_zero(self):
        res = coexpression.permutation_test(ADJ, {"A", "B", "C"}, POP,
                                            n_perm=100, seed=1)
        self.assertGreater(res["p_value"], 0.0)
        self.assertLessEqual(res["p_value"], 1.0)

    def test_reproducible_with_seed(self):
        r1 = coexpression.permutation_test(ADJ, {"A", "C"}, POP, n_perm=200, seed=42)
        r2 = coexpression.permutation_test(ADJ, {"A", "C"}, POP, n_perm=200, seed=42)
        self.assertEqual(r1["p_value"], r2["p_value"])

    def test_requires_at_least_two_genes(self):
        with self.assertRaises(ValueError):
            coexpression.permutation_test(ADJ, {"A"}, POP, n_perm=10)


if __name__ == "__main__":
    unittest.main()
