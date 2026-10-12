#!/usr/bin/env python
"""Plot the top enriched terms from an enrichment result file.

Reads a TSV produced by app.py (or the Results/res_* files) and draws a
horizontal bar chart of the most enriched terms by -log10(p-value).

Usage:
    python app.py -q query/set1.txt -t Pathway -c -o out.tsv
    python plot_enrichment.py out.tsv -o docs/figures/pathway.png --title "Enriched pathways"

The input is tab-separated with columns:
    id   p-value   common   target   description   genes
Lines that are not valid data rows (headers, py2neo log lines) are skipped.
"""

import argparse
import csv
import math
import re
import sys

import matplotlib
matplotlib.use("Agg")  # headless-safe (CI, servers)
import matplotlib.pyplot as plt


def _clean(text):
    """Strip simple HTML tags that EcoCyc descriptions sometimes carry."""
    return re.sub(r"</?[a-zA-Z]+>", "", text).strip()


def parse_results(path, top=10, alpha=0.05):
    """Return [(label, neg_log10_p, p)] for the top significant rows.

    Columns (TSV): id  p-value  common/target  description  genes
    """
    rows = []
    with open(path, newline="") as f:
        for parts in csv.reader(f, delimiter="\t"):
            if len(parts) < 3:
                continue
            term_id, pval, ratio = parts[0], parts[1], parts[2]
            try:
                p = float(pval)
            except ValueError:
                continue  # skip header / py2neo log lines
            if p > alpha:
                continue  # keep only significant terms
            desc = _clean(parts[3]) if len(parts) > 3 and parts[3] else term_id
            # floor tiny p-values so -log10 stays finite on the plot
            neg_log10 = -math.log10(max(p, 1e-300))
            label = desc if len(desc) <= 45 else desc[:42] + "..."
            rows.append((f"{label}  ({ratio})", neg_log10, p))
    rows.sort(key=lambda x: x[1], reverse=True)
    return rows[:top]


def plot(rows, out, title):
    if not rows:
        sys.exit("No plottable rows found in the input file.")
    labels = [r[0] for r in rows][::-1]
    values = [r[1] for r in rows][::-1]

    fig, ax = plt.subplots(figsize=(9, 0.5 * len(rows) + 1.5))
    bars = ax.barh(labels, values, color="#2a7ab9")
    ax.set_xlabel(r"$-\log_{10}(p\mathrm{-value})$")
    ax.set_title(title)
    ax.margins(x=0.12)
    for bar, (_, _, p) in zip(bars, rows[::-1]):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f"  {p:.1e}", va="center", fontsize=8, color="#333")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


def main():
    ap = argparse.ArgumentParser(description="Plot top enriched terms from an app.py result file.")
    ap.add_argument("results", help="TSV result file (from app.py -o).")
    ap.add_argument("-o", "--output", default="enrichment.png", help="Output image path.")
    ap.add_argument("-n", "--top", type=int, default=10, help="Number of terms to show.")
    ap.add_argument("-t", "--title", default="Top enriched terms", help="Figure title.")
    ap.add_argument("-a", "--alpha", type=float, default=0.05, help="Significance cutoff.")
    args = ap.parse_args()
    plot(parse_results(args.results, args.top, args.alpha), args.output, args.title)


if __name__ == "__main__":
    main()
