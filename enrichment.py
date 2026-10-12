"""Statistical enrichment measures.

Pure functions, independent of Neo4j, so they can be unit-tested in isolation.

Notation used throughout:
    k  = number of query genes that also belong to the target set (successes/overlap)
    t  = size of the target set in the whole genome
    q  = size of the query set
    N  = population size (total number of genes in the database)
"""

from scipy.stats import binom, hypergeom, chi2_contingency

MEASURES = ("binomial", "hypergeometric", "chi2", "coverage")


def coverage(k, t, q, N):
    """Simple overlap score (not a p-value). Lower is better/more enriched."""
    if q == 0 or t == 0:
        return 1.0
    return 1 - (k / q) * (k / t)


def binomial(k, t, q, N):
    """Right-tail binomial probability P(X >= k).

    Computed on the complementary event (failures) for numerical stability
    when the enrichment is strong and the naive 1 - cdf(k-1) underflows.
    """
    p = t / N  # probability a random draw hits the target set
    return binom.cdf(q - k, q, 1 - p)


def hypergeometric(k, t, q, N):
    """Right-tail hypergeometric probability P(X >= k) (sampling without replacement)."""
    return 1 - hypergeom.cdf(k - 1, N, t, q)


def chi2(k, t, q, N):
    """Chi-squared test on the 2x2 contingency table, returns the p-value."""
    observed = [
        [k, q - k],
        [t - k, N - (t + q - k)],
    ]
    _, pvalue, _, _ = chi2_contingency(observed)
    return pvalue


_DISPATCH = {
    "binomial": binomial,
    "hypergeometric": hypergeometric,
    "chi2": chi2,
    "coverage": coverage,
}


def score(measure, k, t, q, N):
    """Dispatch to the requested measure. Raises ValueError on unknown measure."""
    try:
        fn = _DISPATCH[measure]
    except KeyError:
        raise ValueError(f"unknown measure {measure!r}, expected one of {MEASURES}")
    return fn(k, t, q, N)


def fdr_filter(results, alpha, query_size):
    """Benjamini-Hochberg style filter over a list of result dicts.

    `results` must be sorted by ascending p-value. Each item needs a 'p-value' key.
    Returns the subset whose adjusted p-value passes `alpha`.
    """
    kept = []
    for i, r in enumerate(results):
        adjusted = (r["p-value"] * query_size) / (i + 1)
        if adjusted <= alpha:
            kept.append(r)
    return kept
