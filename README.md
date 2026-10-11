# E. coli Functional Enrichment on a Neo4j Knowledge Graph

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

## Project layout

```
.
├── app.py                     # main enrichment tool (Cypher + stats)
├── scripts/
│   ├── GeneOntology.py        # OBO parser → GO nodes & is_a/part_of edges
│   ├── gbff.to.bed.py         # genome annotation → BED (gene coordinates)
│   ├── gbff.to.bnumber.GeneID.py  # bnumber ↔ NCBI GeneID mapping
│   ├── ncbi.py                # gene ↔ PubMed citation mapping
│   └── proba_binom1.py        # binomial enrichment prototype
├── *.Rmd                      # UniProt / similarity practicals — course material, see note below
├── neo4j.import/              # ready-to-load CSV/TSV (graph nodes & edges)
├── query/                     # example gene sets (set1-3, query6)
├── ref.sets/                  # reference annotation sets (JSON)
├── downloads/                 # raw sources (NCBI, EcoCyc, UniProt, GO, STRING)
└── Results/                   # enrichment outputs + Revigo tables
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

> **Note on the `.Rmd` files** — `TP01/02/03_*.Rmd` are UniProt / sequence-similarity practicals
> from the course and are kept here only for reference. They are not part of the enrichment pipeline.

---

## Requirements

- Python 3.9+
- A running Neo4j instance (bolt on `localhost:7687`)
- Python packages: `py2neo`, `scipy`

```bash
pip install -r requirements.txt
```

## Setup & run

**1. Load the graph.** Point Neo4j's import directory at `neo4j.import/` and load the nodes
and edges from the provided CSV/TSV files.

**2. Set the connection.** Credentials are read from environment variables (never hard-code them):

```bash
export NEO4J_URI="bolt://localhost:7687"
export NEO4J_USER="neo4j"
export NEO4J_PASSWORD="your-password"
```

**3. Run an analysis:**

```bash
python app.py -q query/set1.txt -t GOTerm -c -v
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
| `-v, --verbose` | Print queries and sets | off |

## Data sources

NCBI (`GCF_000005845.2`), EcoCyc, UniProt (proteome `UP000000625`),
Gene Ontology (`go-basic.obo`), STRING v12.

## Author

Julie Campos — M2 Bioinformatics, Université Paul Sabatier.
