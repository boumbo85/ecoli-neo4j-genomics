#!/usr/bin/env python
"""Positive control: does the enrichment recover a known operon?

The trp operon (trpA-E, b1260-b1264) is a textbook example — five genes that
build L-tryptophan. If the method works, this gene set must come back enriched
for the "L-tryptophan biosynthesis" pathway (TRPSYN-PWY) with a tiny p-value.

This check runs WITHOUT Neo4j: it reads the same csv/ edge files the graph is
built from, so it is fully reproducible in CI. It mirrors what app.py computes
via Cypher (gene -> Pathway membership), then scores with enrichment.py.

Exit code 0 if the expected pathway is the top hit, 1 otherwise.
"""

import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import enrichment

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "csv")

QUERY = {"b1260", "b1261", "b1262", "b1263", "b1264"}  # trp operon
EXPECTED_PATHWAY = "TRPSYN-PWY"                          # L-tryptophan biosynthesis


def load_pathway_membership():
    """pathway_id -> set(genes), and the total gene population size."""
    members = defaultdict(set)
    with open(os.path.join(CSV, "ecocyc.pathways.genes.csv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            members[row["ecocyID"]].add(row["bnumber"])
    # population = all genes in the genome mapping
    genes = set()
    with open(os.path.join(CSV, "mapping.bed.tsv")) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genes.add(row["bnumber"])
    return members, len(genes)


def descriptions():
    desc = {}
    with open(os.path.join(CSV, "ecocyc.pathways.csv")) as f:
        for row in csv.DictReader(f):
            desc[row["ecocyID"]] = row["desc"]
    return desc


def main():
    members, N = load_pathway_membership()
    desc = descriptions()
    q = len(QUERY)

    results = []
    for pid, genes in members.items():
        common = genes & QUERY
        if not common:
            continue
        p = enrichment.binomial(k=len(common), t=len(genes), q=q, N=N)
        results.append({"id": pid, "p-value": p,
                        "common.n": len(common), "target.n": len(genes),
                        "desc": desc.get(pid, "")})

    results.sort(key=lambda r: r["p-value"])
    results = enrichment.fdr_filter(results, alpha=0.05)

    print(f"Positive control — trp operon ({q} genes), population N={N}\n")
    print(f"{'pathway':<16}{'p-value':<12}{'q-value':<12}{'overlap':<9}description")
    for r in results[:5]:
        print(f"{r['id']:<16}{r['p-value']:<12.2e}{r.get('q-value', 0):<12.2e}"
              f"{r['common.n']}/{r['target.n']:<6} {r['desc']}")

    if not results:
        print("\nFAIL: no significant pathway found.")
        return 1
    top = results[0]
    if top["id"] == EXPECTED_PATHWAY:
        print(f"\nPASS: top hit is {EXPECTED_PATHWAY} "
              f"({top['desc']}), p = {top['p-value']:.2e}")
        return 0
    print(f"\nFAIL: expected {EXPECTED_PATHWAY} as top hit, got {top['id']}.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
