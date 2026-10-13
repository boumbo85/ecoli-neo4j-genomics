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


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR adjustment.

    Returns q-values in the SAME order as the input p-values.

    The adjustment is p_(i) * m / rank_i, where m is the number of tests and
    rank_i the ascending rank of p_(i). The classic monotonicity constraint is
    enforced (a q-value is never larger than a higher-ranked one) and values
    are capped at 1.

    Uses statsmodels when available (reference implementation); otherwise falls
    back to an equivalent pure-Python version.

    Note on the earlier version of this code: the correction read
        adjusted = p * query_size / rank
    i.e. it multiplied by the *size of the query gene set* instead of by *m, the
    number of terms tested*. That is wrong — BH scales by the number of
    hypotheses, not by the sample size — and it explains the project report's
    observation that adjusted and raw result counts came out almost identical.
    """
    try:
        from statsmodels.stats.multitest import multipletests
        _, qvals, _, _ = multipletests(pvalues, method="fdr_bh")
        return list(qvals)
    except Exception:
        pass

    m = len(pvalues)
    if m == 0:
        return []
    # rank p-values ascending, apply p * m / rank, then enforce monotonicity
    order = sorted(range(m), key=lambda i: pvalues[i])
    qvals = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):          # from largest p-value down to smallest
        idx = order[rank - 1]
        q = pvalues[idx] * m / rank
        prev = min(prev, q)               # monotone non-decreasing with p-value
        qvals[idx] = min(prev, 1.0)
    return qvals


def fdr_filter(results, alpha):
    """Keep results whose Benjamini-Hochberg q-value passes `alpha`.

    Each item is a dict with a 'p-value' key; the surviving items are returned
    with an added 'q-value' key. Input order is preserved.
    """
    if not results:
        return []
    qvals = benjamini_hochberg([r["p-value"] for r in results])
    kept = []
    for r, q in zip(results, qvals):
        if q <= alpha:
            r = dict(r, **{"q-value": q})
            kept.append(r)
    return kept
