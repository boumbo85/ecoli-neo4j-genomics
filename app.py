#!/usr/bin/env python
"""Search enriched terms/categories in a provided gene set.

Queries a Neo4j knowledge graph of E. coli K-12 MG1655 and reports which
annotation categories (GO terms, pathways, transcription units, protein
families, keywords, literature) are over-represented in the query gene set.

Credentials are read from environment variables:
    export NEO4J_URI="bolt://localhost:7687"
    export NEO4J_USER="neo4j"
    export NEO4J_PASSWORD="your-password"
"""

import argparse
import csv
import os
import sys
from os.path import isfile

from py2neo import Graph
from py2neo.matching import NodeMatcher

import enrichment


# ── Cypher templates per target type ──────────────────────────────────────
# Each entry describes how to reach a target node from a Gene, and whether the
# target carries a human-readable description. Parameters ($genes, $set_id) are
# passed to Neo4j rather than interpolated into the query string.
TARGETS = {
    "Keyword":  {"rel": "-->",                        "desc": False},
    "GOTerm":   {"rel": "-[:is_a|part_of|annotates*]->", "desc": True},
    "TU":       {"rel": "-[:harbors]->",              "desc": True},
    "Pathway":  {"rel": "-[:requires]->",             "desc": True},
    "InterPro": {"rel": "-[:encoded_in]->",           "desc": False},
    "PubMed":   {"rel": "-[:cites]->",                "desc": False},
}


def connect():
    return Graph(
        os.environ.get("NEO4J_URI", "bolt://localhost:7687"),
        auth=(
            os.environ.get("NEO4J_USER", "neo4j"),
            os.environ.get("NEO4J_PASSWORD", "neo4j"),
        ),
    )


def load_query(text):
    """Load a gene set from a file (whitespace/newline separated) or a string."""
    genes = set()
    if isfile(text):
        with open(text) as f:
            genes |= set(f.read().split())
    else:
        genes |= set(text.split())
    return genes


def candidate_terms(neo, target, genes):
    """Return [(id, desc)] of target nodes linked to any gene in the query."""
    spec = TARGETS[target]
    ret = "k.id, k.desc" if spec["desc"] else "k.id, NULL"
    cypher = (
        f"MATCH (k:{target}){spec['rel']}(g:Gene) "
        f"WHERE g.id IN $genes RETURN DISTINCT {ret}"
    )
    return [(row[0], row[1]) for row in neo.run(cypher, genes=list(genes))]


def genes_of_term(neo, target, set_id):
    """Return the set of gene ids associated with a given target node."""
    spec = TARGETS[target]
    cypher = (
        f"MATCH (t:{target}){spec['rel']}(g:Gene) "
        f"WHERE t.id = $set_id RETURN DISTINCT g.id"
    )
    return {row[0] for row in neo.run(cypher, set_id=set_id)}


def main():
    parser = argparse.ArgumentParser(
        description="Search enriched terms/categories in the provided (gene) set"
    )
    parser.add_argument("-q", "--query", required=True,
                        help="Query set: a file of gene ids, or a quoted string.")
    parser.add_argument("-t", "--type", required=True, choices=list(TARGETS),
                        help="Type of target sets.")
    parser.add_argument("-a", "--alpha", type=float, default=0.05,
                        help="Significance threshold.")
    parser.add_argument("-c", "--adjust", action="store_true",
                        help="Adjust for multiple testing (FDR).")
    parser.add_argument("-m", "--measure", default="binomial",
                        choices=enrichment.MEASURES,
                        help="Enrichment measure (default: binomial).")
    parser.add_argument("-l", "--limit", type=int, default=0,
                        help="Maximum number of results to report (0 = all).")
    parser.add_argument("-o", "--output",
                        help="Write results as TSV to this file (default: stdout).")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Talk a lot.")
    param = parser.parse_args()

    neo = connect()
    nodes = NodeMatcher(neo)

    query = load_query(param.query)
    query_size = len(query)
    N = nodes.match("Gene").count()

    if param.verbose:
        print(f"query set ({query_size} genes): {sorted(query)}", file=sys.stderr)
        print(f"population size: {N}", file=sys.stderr)

    # Score every candidate term
    results = []
    for set_id, desc in candidate_terms(neo, param.type, query):
        elements = genes_of_term(neo, param.type, set_id)
        common = elements & query
        k, t = len(common), len(elements)
        pvalue = enrichment.score(param.measure, k, t, query_size, N)
        results.append({
            "id": set_id,
            "p-value": pvalue,
            "common.n": k,
            "target.n": t,
            "desc": desc if desc is not None else "",
            "elements.common": sorted(common),
        })

    results.sort(key=lambda x: x["p-value"])

    if param.adjust:
        results = enrichment.fdr_filter(results, param.alpha)

    if param.limit and param.limit > 0:
        results = results[: param.limit]

    if param.verbose:
        print(f"reporting {len(results)} term(s)", file=sys.stderr)

    # Output ───────────────────────────────────────────────────────────────
    out = open(param.output, "w", newline="") if param.output else sys.stdout
    try:
        writer = csv.writer(out, delimiter="\t")
        header = ["id", "p-value", "common", "target", "description", "genes"]
        if param.adjust:
            header.insert(2, "q-value")  # BH-adjusted
        writer.writerow(header)
        for r in results:
            row = [r["id"], r["p-value"], r["common.n"], r["target.n"],
                   r["desc"], ", ".join(r["elements.common"])]
            if param.adjust:
                row.insert(2, r.get("q-value", ""))
            writer.writerow(row)
    finally:
        if param.output:
            out.close()


if __name__ == "__main__":
    main()
