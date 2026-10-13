"""Co-expression cohesion of a gene set, from STRING edges.

Answers: "are the genes in my query more co-expressed *with each other* than a
random set of the same size?" — a sanity check that a query is a real functional
module rather than an arbitrary list.

Pure functions (no I/O, no Neo4j) so they can be unit-tested. The co-expression
network is passed as an adjacency dict:  {gene: {neighbour: score, ...}, ...}
assumed symmetric (undirected STRING edges).
"""

import random


def intra_set_strength(adjacency, gene_set):
    """Sum of STRING scores of edges whose *both* endpoints are in gene_set.

    Each undirected edge is counted once.
    """
    genes = set(gene_set)
    total = 0
    seen = set()
    for g in genes:
        for nbr, score in adjacency.get(g, {}).items():
            if nbr in genes:
                key = (g, nbr) if g < nbr else (nbr, g)
                if key not in seen:
                    seen.add(key)
                    total += score
    return total


def intra_set_edges(adjacency, gene_set):
    """Number of edges with both endpoints inside gene_set (each counted once)."""
    genes = set(gene_set)
    seen = set()
    for g in genes:
        for nbr in adjacency.get(g, {}):
            if nbr in genes:
                key = (g, nbr) if g < nbr else (nbr, g)
                seen.add(key)
    return len(seen)


def permutation_test(adjacency, gene_set, population, n_perm=1000,
                     statistic=intra_set_strength, seed=0):
    """Empirical p-value for the cohesion of `gene_set`.

    Draws `n_perm` random gene sets of the same size from `population`, computes
    the statistic on each, and returns how often a random set reaches or beats
    the observed value.

    Returns a dict: observed, mean_random, p_value, fold (obs / mean_random).
    The p-value uses the standard (b + 1) / (n_perm + 1) estimator so it is
    never exactly zero.
    """
    genes = set(gene_set)
    population = list(population)
    size = len(genes)
    if size < 2:
        raise ValueError("gene_set must contain at least 2 genes")

    observed = statistic(adjacency, genes)

    rng = random.Random(seed)
    at_least = 0
    rand_total = 0.0
    for _ in range(n_perm):
        sample = rng.sample(population, size)
        val = statistic(adjacency, sample)
        rand_total += val
        if val >= observed:
            at_least += 1

    mean_random = rand_total / n_perm if n_perm else 0.0
    p_value = (at_least + 1) / (n_perm + 1)
    fold = (observed / mean_random) if mean_random > 0 else float("inf")
    return {
        "observed": observed,
        "mean_random": mean_random,
        "p_value": p_value,
        "fold": fold,
        "n_perm": n_perm,
    }
