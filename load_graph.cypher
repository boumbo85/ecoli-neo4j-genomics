// ════════════════════════════════════════════════════════════════════════
// Load the E. coli K-12 MG1655 knowledge graph into Neo4j.
//
// Expects the CSV/TSV files from csv/ to be reachable via file:/// — i.e.
// copied into Neo4j's import directory. Use load_graph.sh, which does the
// copy and then runs this script with cypher-shell.
//
// Relationship directions match what app.py queries:
//   (GOTerm)-[:is_a|part_of|annotates]->(Gene)   note: is_a/part_of are GO→GO
//   (Pathway)-[:requires]->(Gene)
//   (TU)-[:harbors]->(Gene)
//   (InterPro)-[:encoded_in]->(Gene)
//   (Keyword)-->(Gene)        (untyped match in app.py; loaded as :describes)
//   (PubMed)-[:cites]->(Gene)
// Expected counts (from the original build) are noted inline.
// ════════════════════════════════════════════════════════════════════════

// Clean slate ─────────────────────────────────────────────────────────────
MATCH (n) DETACH DELETE n;

// ── Gene nodes (4315) ──────────────────────────────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///mapping.bed.tsv' AS row FIELDTERMINATOR '\t'
CREATE (n:Gene)
SET n.id = row.bnumber,
    n.organism = 511145,
    n.rank = toInteger(row.rank),
    n.strand = row.strand,
    n.begin = toInteger(row.begin),
    n.end = toInteger(row.end);
CREATE INDEX gene_id IF NOT EXISTS FOR (n:Gene) ON (n.id);

// ── GO term nodes (42887) + hierarchy ───────────────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///go.nodes.tsv' AS row FIELDTERMINATOR '\t'
CREATE (n:GOTerm)
SET n.id = row.id,
    n.name = row.desc,
    n.desc = row.def,
    n.namespace = row.namespace;
CREATE INDEX goterm_id IF NOT EXISTS FOR (n:GOTerm) ON (n.id);

LOAD CSV WITH HEADERS FROM 'file:///go.is_a.edges.tsv' AS line FIELDTERMINATOR '\t'
MATCH (t1:GOTerm {id: line.term1}), (t2:GOTerm {id: line.term2})
MERGE (t2)-[:is_a]->(t1);

LOAD CSV WITH HEADERS FROM 'file:///go.part_of.edges.tsv' AS line FIELDTERMINATOR '\t'
MATCH (t1:GOTerm {id: line.term1}), (t2:GOTerm {id: line.term2})
MERGE (t2)-[:part_of]->(t1);

// GOTerm -> Gene annotations (22373)
LOAD CSV WITH HEADERS FROM 'file:///uniprot.GOTerm.bnumber.csv' AS line
MATCH (t:GOTerm {id: line.GOTerm}), (g:Gene {id: line.bnumber})
MERGE (t)-[:annotates]->(g);

// ── Keyword nodes (385) + relations (27323) ─────────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///uniprot.keywords.csv' AS row
CREATE (n:Keyword) SET n.id = row.keyword;
CREATE INDEX keyword_id IF NOT EXISTS FOR (n:Keyword) ON (n.id);

LOAD CSV WITH HEADERS FROM 'file:///uniprot.keywords.genes.csv' AS line
MATCH (k:Keyword {id: line.keyword}), (g:Gene {id: line.bnumber})
MERGE (k)-[:describes]->(g);

// ── InterPro nodes + relations (encoded_in) ─────────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///uniprot.interpro.csv' AS row
CREATE (n:InterPro) SET n.id = row.keyword;
CREATE INDEX interpro_id IF NOT EXISTS FOR (n:InterPro) ON (n.id);

LOAD CSV WITH HEADERS FROM 'file:///uniprot.interpro.genes.csv' AS line
MATCH (k:InterPro {id: line.keyword}), (g:Gene {id: line.bnumber})
MERGE (k)-[:encoded_in]->(g);

// ── PubMed nodes (35320) + citations (101841) ───────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///ncbi.pmid.csv' AS row
CREATE (n:PubMed) SET n.id = row.PMID;
CREATE INDEX pubmed_id IF NOT EXISTS FOR (n:PubMed) ON (n.id);

LOAD CSV WITH HEADERS FROM 'file:///ncbi.pmid.genes.csv' AS line
MATCH (p:PubMed {id: line.PMID}), (g:Gene {id: line.gene_id})
MERGE (p)-[:cites]->(g);

// ── Transcription Unit nodes (3635) + harbors relations ──────────────────────
LOAD CSV WITH HEADERS FROM 'file:///ecocyc.tu.csv' AS row
CREATE (n:TU)
SET n.id = row.ecocyID,
    n.strand = row.Strand,
    n.begin = toInteger(row.Left),
    n.end = toInteger(row.Right),
    n.desc = row.desc,
    n.organism = 511145;
CREATE INDEX tu_id IF NOT EXISTS FOR (n:TU) ON (n.id);

LOAD CSV WITH HEADERS FROM 'file:///ecocyc.tu.genes.csv' AS line
MATCH (p:TU {id: line.ecocyID}), (g:Gene {id: line.bnumber})
MERGE (p)-[:harbors]->(g);

// ── Pathway nodes (438) + requires relations ─────────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///ecocyc.pathways.csv' AS row
CREATE (n:Pathway)
SET n.id = row.ecocyID,
    n.organism = 511145,
    n.desc = row.desc;
CREATE INDEX pathway_id IF NOT EXISTS FOR (n:Pathway) ON (n.id);

LOAD CSV WITH HEADERS FROM 'file:///ecocyc.pathways.genes.csv' AS line
MATCH (p:Pathway {id: line.ecocyID}), (g:Gene {id: line.bnumber})
MERGE (p)-[:requires]->(g);

// ── STRING coexpression (Gene <-> Gene) ──────────────────────────────────────
LOAD CSV WITH HEADERS FROM 'file:///string.coexpression.csv' AS line
MATCH (g1:Gene {id: line.id1}), (g2:Gene {id: line.id2})
MERGE (g1)-[r:coexpression]->(g2)
SET r.score = toInteger(line.coexpression);

// ── Sanity check ─────────────────────────────────────────────────────────────
MATCH (g:Gene) RETURN count(g) AS genes;
