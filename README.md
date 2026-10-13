# E. coli Functional Enrichment on a Neo4j Knowledge Graph

[![tests](https://github.com/boumbo85/ecoli-neo4j-genomics/actions/workflows/tests.yml/badge.svg)](https://github.com/boumbo85/ecoli-neo4j-genomics/actions/workflows/tests.yml)

> 🇬🇧 English · [🇫🇷 Version française](./README.fr.md)

Give this tool a list of *E. coli* genes, and it tells you what they do together — which
metabolic pathways, GO functions, protein families and even which papers describe them.
Everything runs as Cypher queries against a Neo4j knowledge graph, with the statistics done in Python.

**Context:** M2 Bioinformatics (Université Paul Sabatier), NoSQL / graph-databases module.
**Organism:** *Escherichia coli* K-12 MG1655 (`GCF_000005845.2`).

---

## The problem

A raw list of 100–200 gene identifiers (`b0001 b0002 b0003 …`) carries no biological meaning on
its own. This project turns that list into a ranked, p-value-sorted answer to the question
*"what connects these genes?"* — and then lets you follow the thread all the way to the literature.

```
gene set  →  Neo4j graph  →  enrichment test  →  ranked terms  →  PubMed articles  →  biology
```

## The knowledge graph

A single `Gene` node sits at the centre, linked to six annotation layers by typed relationships:

```
                         ┌─────────────┐
            annotates    │   GOTerm    │  is_a / part_of
          ┌──────────────│  (GO:xxxx)  │──────────────┐
          │              └─────────────┘              │ (hierarchy)
          │                                           ▼
   ┌──────────┐  requires   ┌─────────┐        ┌─────────────┐
   │  Pathway │────────────▶│         │        │   GOTerm    │
   │ (EcoCyc) │             │         │        └─────────────┘
   └──────────┘             │         │
                            │  Gene   │  harbors   ┌─────────┐
   ┌──────────┐  encoded_in │ (bNNNN) │◀───────────│   TU    │
   │ InterPro │────────────▶│         │            │(EcoCyc) │
   └──────────┘             │         │            └─────────┘
                            │         │  cites     ┌─────────┐
   ┌──────────┐  annotates  │         │◀───────────│ PubMed  │
   │ Keyword  │────────────▶│         │            │ (NCBI)  │
   │(UniProt) │             └────┬────┘            └─────────┘
   └──────────┘                  │ coexpression (STRING)
                                 ▼
                              (Gene)
```

| Target node | Relationship to `Gene`           | Source      |
|-------------|----------------------------------|-------------|
| `GOTerm`    | `is_a`, `part_of`, `annotates`   | Gene Ontology (`go-basic.obo`) |
| `Pathway`   | `requires`                       | EcoCyc      |
| `TU` (transcription unit) | `harbors`          | EcoCyc      |
| `InterPro`  | `encoded_in`                     | UniProt     |
| `Keyword`   | annotates                        | UniProt     |
| `PubMed`    | `cites`                          | NCBI        |
| `Gene`→`Gene` | `coexpression`                 | STRING v12  |

## The statistics

Each candidate term is tested for over-representation in the query set versus the whole genome
(4315 genes). Four measures, plus optional FDR correction:

- **Binomial** (default) — computed on failures for numerical stability
- **Hypergeometric** — sampling without replacement
- **Chi²** — contingency-table test
- **Coverage** — simple overlap score

---

## Worked example — from gene list to biology

Running the tool on query set **C** (113 genes):

```bash
python app.py -q query/set1.txt -t Pathway -c
```

**Top enriched pathways (EcoCyc):**

| Pathway ID | p-value   | genes | Description |
|------------|-----------|-------|-------------|
| PWY-8243   | 4.8e-83   | 11/11 | colanic acid (E. coli K12) biosynthesis |
| PWY-66     | 1.4e-10   | 2/2   | GDP-L-fucose biosynthesis I |
| COLANSYN-PWY | 1.3e-09 | 4/11  | colanic acid building blocks |

![Top enriched pathways for query set C](docs/figures/top_pathways.png)

*Generated with `python plot_enrichment.py <results.tsv> -o fig.png`.*

**Top enriched GO terms** point the same way — colanic acid biosynthetic process
(GO:0046377, p = 3.3e-96, 17/23 genes), plus a second cluster around **c-di-GMP signalling**
(diguanylate cyclase GO:0052621, phosphodiesterase GO:0071111) and **negative regulation of
flagellar motility** (GO:1902201).

**Following the thread to the literature** (`-t PubMed`) surfaces PMID **16672475**
(p = 5.2e-279, 35 genes) — the reference describing capsular colanic-acid synthesis.

**So what does it mean?** The gene block **b2044–b2059** encodes the biosynthesis of
**colanic acid**, a capsular exopolysaccharide. Its appearance alongside c-di-GMP enzymes and
motility repressors tells a coherent story: these genes drive the switch from a motile,
planktonic lifestyle toward a **sessile, biofilm-forming state** under stress — the enrichment
analysis recovers that biology straight from an anonymous list of identifiers.

---

## Validation (positive control)

A tool that "runs" is not a tool that "works". To show the method recovers known biology,
`validate.py` feeds it a textbook case: the **trp operon** (`trpA–E`, b1260–b1264), five genes
that build L-tryptophan.

```bash
python validate.py
```

```
pathway         p-value     q-value     overlap  description
TRPSYN-PWY      2.11e-15    1.27e-14    5/5      L-tryptophan biosynthesis
COMPLETE-ARO-PWY 2.16e-12   6.49e-12    5/20     superpathway of aromatic amino acid biosynthesis
ALL-CHORISMATE-PWY 4.44e-10 8.88e-10    5/58     superpathway of chorismate metabolism
PASS: top hit is TRPSYN-PWY (L-tryptophan biosynthesis), p = 2.11e-15
```

The operon comes back as the top hit on L-tryptophan biosynthesis, followed by the aromatic
amino acid and chorismate superpathways — exactly the expected hierarchy, since tryptophan
derives from chorismate. The check runs without Neo4j (it reads the same `csv/` edge files) so
it also runs in CI.

## Co-expression cohesion (STRING)

Enrichment tells you *which* functions a gene set is annotated for. A different, complementary
question: do these genes actually **behave as a module** — are they more co-expressed with each
other than a random set of the same size? `coexpression_analysis.py` answers it with a
permutation test on the STRING co-expression edges.

```bash
python coexpression_analysis.py query/query6.txt -n 2000
```

| Gene set | genes | observed strength | random mean | fold | p-value |
|----------|-------|-------------------|-------------|------|---------|
| study set (query6) | 113 | 305,975 | 23,140 | **13.2×** | ≤ 5e-4 |
| trp operon (control) | 5 | 9,934 | 39 | **252×** | ≤ 5e-4 |

The study set is 13× more internally co-expressed than chance — strong evidence it is a genuine
functional module, not an arbitrary list, which backs up the biofilm interpretation. The trp
operon, co-transcribed by definition, scores 252× as a positive control. This uses the STRING
edges that are otherwise only stored in the graph.

## Methods & design choices

These are the questions a reviewer asks, answered up front.

- **Where do the gene sets come from?** `query/set1.txt` is a *control* set with a known answer,
  used to build and sanity-check the tool. `query/query6.txt` (= "set 6", 113 genes) is the
  study set that was actually analysed. `positive_control_trp.txt` is the trp operon used above.
- **Why binomial by default, not hypergeometric?** Gene-set enrichment is sampling *without*
  replacement, so the exact model is hypergeometric. The binomial approximates it well when the
  query is small relative to the population — the usual rule is `n/N < 0.1`. Here
  `n/N = 113/4315 = 0.026`, so the approximation is legitimate, and binomial/hypergeometric give
  the same hits in practice. Both are available via `--measure`; use `hypergeometric` if you
  prefer the exact test.
- **What is the background?** The whole genome (4315 genes). If your query comes from an assay
  that only probed a subset of genes (e.g. a microarray), the correct background is that subset,
  not the genome — this changes the p-values. Supporting a custom background is on the roadmap.
- **Multiple testing.** p-values are adjusted with Benjamini-Hochberg (`-c`), reported as
  q-values. See the note below.
- **STRING coexpression** powers a second, complementary check (see below): is the query set more
  co-expressed with itself than chance?

> **A fixed bug, worth telling.** The project report already flagged that FDR-adjusted and raw
> result counts came out almost identical "which should not happen — probably a bug". It was: the
> original correction scaled p-values by the *query size* instead of by the *number of tests*.
> `enrichment.py` now implements Benjamini-Hochberg correctly (via `statsmodels`, with an
> equivalent pure-Python fallback), and the unit tests pin the expected q-values.

---

## Project layout

```
.
├── app.py                     # main enrichment tool (CLI + Cypher queries)
├── enrichment.py              # statistical measures + BH FDR (pure, unit-tested)
├── coexpression.py            # STRING co-expression cohesion (pure, unit-tested)
├── coexpression_analysis.py   # permutation test on a gene set, no Neo4j needed
├── validate.py                # positive control (trp operon), no Neo4j needed
├── plot_enrichment.py         # bar chart of top enriched terms
├── load_graph.cypher          # builds the graph from csv/ (LOAD CSV)
├── load_graph.sh              # copies CSVs + runs load_graph.cypher
├── docker-compose.yml         # reproducible Neo4j instance
├── tests/                     # unit tests for the enrichment measures
├── scripts/
│   ├── GeneOntology.py        # OBO parser → GO nodes & is_a/part_of edges
│   ├── gbff.to.bed.py         # genome annotation → BED (gene coordinates)
│   ├── gbff.to.bnumber.GeneID.py  # bnumber ↔ NCBI GeneID mapping
│   ├── ncbi.py                # gene ↔ PubMed citation mapping
│   └── proba_binom1.py        # binomial enrichment prototype
├── csv/                       # ready-to-load CSV/TSV (graph nodes & edges)
├── query/                     # example gene sets (set1-3, query6)
├── Results/                   # enrichment outputs + Revigo tables
└── course-material/           # UniProt / similarity practicals (not part of the pipeline)
```

### Results index

Files are named `res_<Type>_<set>.txt`. Each line is tab-separated:
`term_id   p-value   common/target   description   genes`.

| File pattern            | Content |
|-------------------------|---------|
| `res_GOTerm_*.txt`      | Enriched GO terms (BP/MF/CC) per gene set |
| `res_Pathway_*.txt`     | Enriched EcoCyc pathways |
| `res_TU_*.txt`          | Enriched transcription units |
| `res_InterPro_*.txt`    | Enriched InterPro protein families |
| `res_keywords*`         | Enriched UniProt keywords |
| `res_PubMed_*.txt`      | Over-represented literature (PMIDs) |
| `Revigo_BP/MF_*.tsv`    | Revigo semantic clustering of GO terms (for plotting) |

> **Note on `course-material/`** — `TP01/02/03_*.Rmd` are UniProt / sequence-similarity
> practicals from the course, kept for reference. They are not part of the enrichment pipeline.

---

## Quick start (Docker)

Clone, start Neo4j, build the graph, run — four commands:

```bash
docker compose up -d                 # Neo4j at http://localhost:7474
./load_graph.sh                      # build the graph from csv/
pip install -r requirements.txt
python app.py -q query/set1.txt -t Pathway -c
```

Default local credentials are `neo4j / testpassword` (see `docker-compose.yml`). Override
them for anything shared.

## Manual setup (existing Neo4j)

Point your Neo4j import directory at `csv/`, set the connection, and load the graph:

```bash
export NEO4J_URI="bolt://localhost:7687"
export NEO4J_USER="neo4j"
export NEO4J_PASSWORD="your-password"

NEO4J_IMPORT_DIR=/path/to/neo4j/import ./load_graph.sh
python app.py -q query/set1.txt -t GOTerm -c -v
```

Credentials are read from the environment variables above — nothing is hard-coded.

## Tests

The statistical measures live in `enrichment.py` and are covered by unit tests:

```bash
python -m unittest discover -s tests
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-q, --query`   | Gene set: a file of bnumbers, or a quoted string | required |
| `-t, --type`    | Target: `GOTerm` `Pathway` `TU` `InterPro` `Keyword` `PubMed` | required |
| `-a, --alpha`   | Significance threshold | `0.05` |
| `-c, --adjust`  | Apply FDR correction | off |
| `-m, --measure` | `binomial` `hypergeometric` `chi2` `coverage` | `binomial` |
| `-l, --limit`   | Max results reported | `0` (all) |
| `-o, --output`  | Write results as TSV to a file | stdout |
| `-v, --verbose` | Print queries and sets | off |

### Plotting

Turn any result file into a bar chart of the top terms:

```bash
python app.py -q query/set1.txt -t Pathway -c -o pathways.tsv
python plot_enrichment.py pathways.tsv -o fig.png -t "Top enriched pathways"
```

## Data sources

NCBI (`GCF_000005845.2`), EcoCyc, UniProt (proteome `UP000000625`),
Gene Ontology (`go-basic.obo`), STRING v12.

## Report

The full project report (corrected version) is in
[`report/`](report/CAMPOS_JULIE_M2BBS_GNDS_corrige.pdf) — it covers the database design,
the enrichment method, the biological interpretation, the FDR fix and the co-expression
validation.

## Roadmap / next steps

Things that would push this from a solid academic project to a portfolio piece:

- [x] **Correct FDR** — Benjamini-Hochberg via `statsmodels`, q-values reported, unit-tested.
- [x] **Positive control** — trp operon validation, runs in CI.
- [x] **Use the STRING edges** — co-expression cohesion permutation test (`coexpression_analysis.py`).
- [ ] **Custom background** — accept a background gene list (assay-specific) instead of the whole genome.
- [ ] **End-to-end CI** — spin up Neo4j as a service in GitHub Actions and run one real enrichment query.
- [ ] **Package the CLI** — a `pyproject.toml` + console entry point (`ecoli-enrich ...`).
- [ ] **Generalise beyond E. coli** — parameterise the organism so the same schema works for any annotated genome.

> **Lessons learned (worth stating).** Three habits this project taught me, now applied from commit 1:
> secrets and local paths never go in a repo (env vars + `.gitignore` from the start); any project
> touching an external database ships with its `docker-compose` so it stays reproducible; and
> mathematical functions get unit tests — it is the cheapest, most convincing quality signal.

## Author

Julie Campos — M2 Bioinformatics, Université Paul Sabatier.
