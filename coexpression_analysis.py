#!/usr/bin/env python
"""Is a gene set more co-expressed than chance? (STRING, permutation test)

Uses the STRING co-expression edges that are loaded into the graph but not used
by the enrichment test — this is the "unused data" the README flags.

Runs WITHOUT Neo4j: reads csv/string.coexpression.csv directly, so it is
reproducible in CI. Mirrors the (Gene)-[:coexpression]->(Gene) relation.

Usage:
    python coexpression_analysis.py query/query6.txt
    python coexpression_analysis.py query/positive_control_trp.txt -n 5000
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coexpression

HERE = os.path.dirname(os.path.abspath(__file__))
STRING_CSV = os.path.join(HERE, "csv", "string.coexpression.csv")
GENOME = os.path.join(HERE, "csv", "mapping.bed.tsv")


def load_network(min_score=0):
    """Build a symmetric adjacency dict from the STRING co-expression file."""
    adjacency = defaultdict(dict)
    with open(STRING_CSV) as f:
        for row in csv.DictReader(f):
            score = int(row["coexpression"])
            if score < min_score:
                continue
            a, b = row["id1"], row["id2"]
            adjacency[a][b] = score
            adjacency[b][a] = score
    return adjacency


def load_population():
    genes = set()
    with open(GENOME) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["bnumber"])
    return genes


def load_query(path):
    if os.path.isfile(path):
        with open(path) as f:
            return set(f.read().split())
    return set(path.split())


def main():
    ap = argparse.ArgumentParser(description="Co-expression cohesion of a gene set (STRING permutation test).")
    ap.add_argument("query", help="Gene set: a file of bnumbers or a quoted string.")
    ap.add_argument("-n", "--permutations", type=int, default=2000, help="Number of random draws.")
    ap.add_argument("-s", "--min-score", type=int, default=0, help="Minimum STRING co-expression score to keep an edge.")
    ap.add_argument("--seed", type=int, default=0, help="Random seed (reproducibility).")
    args = ap.parse_args()

    adjacency = load_network(args.min_score)
    population = load_population()
    query = load_query(args.query) & population  # only genes present in the genome

    res = coexpression.permutation_test(
        adjacency, query, population,
        n_perm=args.permutations, seed=args.seed,
    )

    print(f"Co-expression cohesion — {len(query)} genes, {args.permutations} permutations\n")
    print(f"  observed intra-set strength : {res['observed']:.0f}")
    print(f"  mean of random sets         : {res['mean_random']:.0f}")
    print(f"  fold vs random              : {res['fold']:.1f}x")
    print(f"  empirical p-value           : {res['p_value']:.2e}")

    if res["p_value"] <= 0.05:
        print(f"\n=> The set is significantly more co-expressed than chance "
              f"({res['fold']:.1f}x, p={res['p_value']:.2e}): it behaves like a functional module.")
    else:
        print("\n=> No significant co-expression cohesion beyond chance.")


if __name__ == "__main__":
    main()
