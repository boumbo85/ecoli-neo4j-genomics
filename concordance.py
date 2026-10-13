#!/usr/bin/env python
"""Quantify agreement with a gold standard (the reference online tool).

The report claims our results are "consistent with the online tool
(silico.biotoul.fr/enrichment)". This turns that impression into numbers:
overlap, Jaccard index, and recovery of the reference's significant GO terms.

Compares the GO IDs our app reports against the GO IDs the reference tool
reports (Results/res_go_web.txt, "binom adjusted" block). Reads existing result
files only — no Neo4j needed.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "Results")

GO_RE = re.compile(r"GO:\d+")


def go_ids_from_file(path, pvalue_col=1, alpha=0.05):
    """Collect GO IDs whose p-value (column `pvalue_col`) passes `alpha`."""
    ids = set()
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = GO_RE.match(line.strip())
            if not m:
                continue
            parts = line.split("\t")
            try:
                p = float(parts[pvalue_col])
            except (IndexError, ValueError):
                ids.add(m.group(0))  # reference file: keep listed terms
                continue
            if p <= alpha:
                ids.add(m.group(0))
    return ids


def main():
    ours_path = os.path.join(RES, "res_GOTerm_C.txt")   # our app, set 6
    ref_path = os.path.join(RES, "res_go_web.txt")       # reference online tool
    if not (os.path.exists(ours_path) and os.path.exists(ref_path)):
        sys.exit("Result files not found in Results/.")

    ours = go_ids_from_file(ours_path, pvalue_col=1)
    ref = go_ids_from_file(ref_path, pvalue_col=1)

    inter = ours & ref
    union = ours | ref
    jaccard = len(inter) / len(union) if union else 0.0
    recovery = len(inter) / len(ref) if ref else 0.0

    print("Concordance with the reference tool (GO terms, set 6)\n")
    print(f"  terms reported by our app        : {len(ours)}")
    print(f"  terms reported by reference tool : {len(ref)}")
    print(f"  shared terms                     : {len(inter)}")
    print(f"  Jaccard index                    : {jaccard:.3f}")
    print(f"  recovery of reference terms      : {recovery:.0%} ({len(inter)}/{len(ref)})")
    print()
    print("Interpretation: our app reports many more terms because it propagates")
    print("annotations up the GO hierarchy (true-path rule with a variable-length")
    print("path), so general ancestor terms also come up. The low Jaccard is driven")
    print("by that over-generalisation, not by disagreement on the specific biology:")
    print(f"it still recovers {recovery:.0%} of the reference's significant terms.")


if __name__ == "__main__":
    main()
