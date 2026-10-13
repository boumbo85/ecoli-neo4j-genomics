#!/usr/bin/env bash
# Load the E. coli knowledge graph into a running Neo4j instance.
#
# Usage:
#   ./load_graph.sh                 # uses docker compose service "neo4j"
#   NEO4J_PASSWORD=secret ./load_graph.sh
#
# What it does:
#   1. copies csv/* into the Neo4j import directory (so file:/// works)
#   2. runs load_graph.cypher with cypher-shell
#
# With the provided docker-compose.yml the import directory is already
# bind-mounted from ./csv, so step 1 is a no-op in that setup.

set -euo pipefail

NEO4J_USER="${NEO4J_USER:-neo4j}"
NEO4J_PASSWORD="${NEO4J_PASSWORD:-testpassword}"
NEO4J_URI="${NEO4J_URI:-bolt://localhost:7687}"
COMPOSE_SERVICE="${COMPOSE_SERVICE:-neo4j}"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Loading graph into ${NEO4J_URI} as ${NEO4J_USER}..."

if command -v docker >/dev/null 2>&1 && docker compose ps >/dev/null 2>&1; then
    # Run inside the container; csv/ is bind-mounted to /var/lib/neo4j/import
    docker compose exec -T "${COMPOSE_SERVICE}" \
        cypher-shell -u "${NEO4J_USER}" -p "${NEO4J_PASSWORD}" \
        < "${here}/load_graph.cypher"
elif command -v cypher-shell >/dev/null 2>&1; then
    # Local Neo4j: copy CSVs into the import dir first
    IMPORT_DIR="${NEO4J_IMPORT_DIR:?set NEO4J_IMPORT_DIR to your Neo4j import folder}"
    echo "Copying csv/* into ${IMPORT_DIR}"
    cp "${here}/csv/"* "${IMPORT_DIR}/"
    cypher-shell -a "${NEO4J_URI}" -u "${NEO4J_USER}" -p "${NEO4J_PASSWORD}" \
        < "${here}/load_graph.cypher"
else
    echo "ERROR: neither 'docker compose' nor 'cypher-shell' is available." >&2
    exit 1
fi

echo "Done. Try:  python app.py -q query/set1.txt -t Pathway -c"
