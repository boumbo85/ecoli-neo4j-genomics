#!/usr/bin/env python
"""Negative control: a random gene set should NOT be significantly enriched.

The positive control (validate.py) shows the method finds known biology. The
negative control shows the symmetric property: fed random genes, the method
should return (almost) nothing significant. If random input produced strong
"enrichment", the test would be broken (inflated false positives).

Runs without Neo4j: reads csv/ edge files, scores EcoCyc pathways with the
binomial test + Benjamini-Hochberg, and checks how many random sets of size n
yield any significant pathway.

Exit code 0 if the false-positive rate stays near the expected level.
"""

import csv
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import enrichment

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "csv")

SET_SIZE = 113       # same size as the study set (query6)
N_TRIALS = 200       # number of random sets to draw
ALPHA = 0.05
SEED = 0


def load_pathways():
    members = defaultdict(set)
    with open(os.path.join(CSV, "ecocyc.pathways.genes.csv")) as f:
        for row in csv.DictReader(f):
            members[row["ecocyID"]].add(row["bnumber"])
    genes = set()
    with open(os.path.join(CSV, "mapping.bed.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["bnumber"])
    return members, sorted(genes)


def significant_pathways(members, query, N):
    q = len(query)
    results = []
    for pid, genes in members.items():
        common = genes & query
        if not common:
            continue
        p = enrichment.binomial(k=len(common), t=len(genes), q=q, N=N)
        results.append({"id": pid, "p-value": p})
    results.sort(key=lambda r: r["p-value"])
    return enrichment.fdr_filter(results, alpha=ALPHA)


def main():
    members, population = load_pathways()
    N = len(population)
    rng = random.Random(SEED)

    hits = 0           # random sets with >=1 significant pathway
    total_terms = 0
    for _ in range(N_TRIALS):
        query = set(rng.sample(population, SET_SIZE))
        sig = significant_pathways(members, query, N)
        if sig:
            hits += 1
            total_terms += len(sig)

    rate = hits / N_TRIALS
    print(f"Negative control — {N_TRIALS} random sets of {SET_SIZE} genes, population N={N}\n")
    print(f"  random sets with >=1 significant pathway : {hits}/{N_TRIALS} ({rate:.1%})")
    print(f"  mean significant terms per random set    : {total_terms / N_TRIALS:.2f}")
    print(f"  (study set query6, for comparison, has 11 significant pathways)\n")

    # A sound test keeps the false-positive rate modest. We allow a generous
    # ceiling because random E. coli genes still occasionally overlap a pathway,
    # but it must be far below the study set's 11 hits.
    if rate <= 0.25 and total_terms / N_TRIALS < 1.0:
        print("PASS: random input yields little to no enrichment, as expected.")
        return 0
    print("FAIL: random input produces too much 'enrichment' — the test is inflated.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
