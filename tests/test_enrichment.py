"""Unit tests for the statistical enrichment measures.

Run with:  python -m pytest    (or)    python -m unittest
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import enrichment  # noqa: E402


class TestCoverage(unittest.TestCase):
    def test_perfect_overlap_scores_zero(self):
        # query fully inside target and target fully covered → best score
        self.assertAlmostEqual(enrichment.coverage(k=5, t=5, q=5, N=100), 0.0)

    def test_no_overlap_scores_one(self):
        self.assertAlmostEqual(enrichment.coverage(k=0, t=5, q=5, N=100), 1.0)

    def test_empty_sets_are_safe(self):
        self.assertEqual(enrichment.coverage(k=0, t=0, q=0, N=100), 1.0)


class TestBinomial(unittest.TestCase):
    def test_returns_probability(self):
        p = enrichment.binomial(k=3, t=10, q=5, N=100)
        self.assertTrue(0.0 <= p <= 1.0)

    def test_strong_enrichment_is_significant(self):
        # 10/10 query genes in a tiny target set → very small p-value
        p = enrichment.binomial(k=10, t=12, q=10, N=4315)
        self.assertLess(p, 1e-10)

    def test_more_overlap_lowers_pvalue(self):
        weak = enrichment.binomial(k=2, t=50, q=10, N=1000)
        strong = enrichment.binomial(k=8, t=50, q=10, N=1000)
        self.assertLess(strong, weak)


class TestHypergeometric(unittest.TestCase):
    def test_returns_probability(self):
        p = enrichment.hypergeometric(k=3, t=10, q=5, N=100)
        self.assertTrue(0.0 <= p <= 1.0)

    def test_full_overlap_is_significant(self):
        p = enrichment.hypergeometric(k=10, t=10, q=10, N=4315)
        self.assertLess(p, 1e-10)

    def test_single_success_equals_prob_at_least_one(self):
        # P(X>=1) = 1 - P(X=0); here with k=1 the right tail is essentially 1
        p = enrichment.hypergeometric(k=1, t=50, q=10, N=100)
        self.assertGreater(p, 0.9)


class TestChi2(unittest.TestCase):
    def test_returns_probability(self):
        p = enrichment.chi2(k=5, t=20, q=20, N=200)
        self.assertTrue(0.0 <= p <= 1.0)

    def test_strong_association_is_significant(self):
        p = enrichment.chi2(k=18, t=20, q=20, N=200)
        self.assertLess(p, 0.05)


class TestDispatch(unittest.TestCase):
    def test_dispatch_matches_direct_call(self):
        direct = enrichment.binomial(k=3, t=10, q=5, N=100)
        vianame = enrichment.score("binomial", k=3, t=10, q=5, N=100)
        self.assertEqual(direct, vianame)

    def test_unknown_measure_raises(self):
        with self.assertRaises(ValueError):
            enrichment.score("nope", k=1, t=1, q=1, N=1)


class TestFdrFilter(unittest.TestCase):
    def test_keeps_significant_sorted_results(self):
        results = [
            {"p-value": 1e-10},
            {"p-value": 1e-5},
            {"p-value": 0.9},
        ]
        kept = enrichment.fdr_filter(results, alpha=0.05, query_size=10)
        self.assertEqual(len(kept), 2)

    def test_empty_input(self):
        self.assertEqual(enrichment.fdr_filter([], alpha=0.05, query_size=10), [])


if __name__ == "__main__":
    unittest.main()
