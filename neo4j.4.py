#!/bin/env python
# conda install python3 biopython requests
# ./gbff.to.bnumber.GeneID.py > mapping.bnumber_ncbi.tsv

# UPLOAD PACKAGES
from py2neo import * # conda install py2neo, monotonic
from py2neo.matching import *
import argparse
from pprint import pprint
from os.path import isfile
import json
from scipy.stats import binom, hypergeom, chi2_contingency  # pip install scipy


# CONNECTION TO NEO4J
neo = Graph("bolt://localhost:7687",
            auth=("neo4j","omybioinfo"))
nodes = NodeMatcher(neo)
 

# SCRIPT PARAMETERS #######################################################
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set (gene nodes in Neo4j).')
parser.add_argument('-t', '--type', required=True, choices=['Keyword', 'Pathway','TU','GOTerm','InterPro','PubMed'], help='Type of target sets: keyword, pathway,TU,GOterm,Interpro,pubmed')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')



# Path Query
path = '../query/set1.txt'
param = parser.parse_args(['-q', path, 
                           '-t', 'PubMed',
                           '-c', '-v', '-l', '10'])

# LOAD QUERY #######################################################################
text = param.query
query = set()
if isfile(text):
    with open(text) as f:
        content = ' '.join(f.read().split('\n')).split()
        query |= set(content)
else: # parse string
    query |= set(text.split())

if param.verbose:
  print(f'query set: {query}')

# LOAD TARGET/REFERENCE SETS #######################################################
population = nodes.match("Gene")
pprint(population)

# POPULATION SIZE = total number of genes in the database
n_gene_total = population.count() # len(population)
pprint(f"population size: {n_gene_total}")

# QUERY SIZE = number of genes in the query
query_size = len(query)
pprint(f"Number of genes in the query: {query_size}")

# EXTRACT query genes from the database:
target = param.type
if target == 'Keyword': # if we want to search for keywords
  cypher = f"MATCH (k:{target})-->(g:Gene) WHERE g.id IN {list(query)} RETURN DISTINCT k"
  res = neo.run(cypher)#.to_table()

elif target == 'GOTerm': # EXTRACT GO TERMS FROM query genes from the database:
  cypher = f"MATCH (k:{target})-[:is_a|part_of|annotates*]->(g:Gene) WHERE g.id IN {list(query)} RETURN DISTINCT k.id, k.desc" #:is_a|part_of|annotates a la place on peut mettre *
  res = neo.run(cypher).to_table()

elif target == 'TU': # EXTRACT TU FROM query genes from the database:
  cypher = f"MATCH (k:{target})-[:harbors]->(g:Gene) WHERE g.id IN {list(query)} RETURN DISTINCT k.id, k.desc" 
  res = neo.run(cypher).to_table()

elif target == 'Pathway': # EXTRACT PATHWAY FROM query genes from the database:
  cypher = f"MATCH (k:{target})-[:requires]->(g:Gene) WHERE g.id IN {list(query)} RETURN DISTINCT k.id, k.desc" 
  res = neo.run(cypher).to_table()

elif target == 'InterPro': # EXTRACT Interpro FROM query genes from the database:
  cypher = f"MATCH (k:{target})-[:encoded_in]->(g:Gene) WHERE g.id IN {list(query)} RETURN DISTINCT k.id" 
  res = neo.run(cypher).to_table()

elif target == 'PubMed': # EXTRACT PubMed FROM query genes from the database:
  cypher = f"MATCH (k:{target})-[:cites]->(g:Gene) WHERE g.id IN {list(query)} RETURN DISTINCT k.id" 
  res = neo.run(cypher).to_table()

if param.verbose:
  pprint(f'Cyper request: {cypher}')
  pprint(f'target set: {res}')


#  STATISTICAL TESTS #######################################################################
results = [] 
liste_proba = [] # list of dictionaries to store results (probabilities)
dico_proba = {}
for r in res:
    # EXTRACT GENES ASSOCIATED (bnumber) WITH CYPHER REQUESTS
    if target == 'Keyword':
        set_id = r['k']['id']
        desc = []
        cypher2 = f"MATCH (t:{target})-->(g:Gene) WHERE t.id = '{set_id}' RETURN DISTINCT g.id" 
    elif target == 'GOTerm':   
        set_id = r[0] 
        desc = r[1]
        cypher2 = f"MATCH (t:{target})-[:is_a|part_of|annotates*]->(g:Gene) WHERE t.id = '{set_id}' RETURN DISTINCT g.id"
    elif target == 'TU':   
        set_id = r[0] 
        cypher2 = f"MATCH (t:{target})-[:harbors]->(g:Gene) WHERE t.id = '{set_id}' RETURN DISTINCT g.id"
    elif target == 'Pathway':   
        set_id = r[0] 
        cypher2 = f"MATCH (t:{target})-[:requires]->(g:Gene) WHERE t.id = '{set_id}' RETURN DISTINCT g.id"    
    elif target == 'InterPro':   
        set_id = r[0] 
        cypher2 = f"MATCH (t:{target})-[:encoded_in]->(g:Gene) WHERE t.id = '{set_id}' RETURN DISTINCT g.id"    
    elif target == 'PubMed':   
        set_id = r[0] 
        cypher2 = f"MATCH (t:{target})-[:cites]->(g:Gene) WHERE t.id = '{set_id}' RETURN DISTINCT g.id"    
    
    # EXECUTE CYPHER REQUESTS ACCORDING TO THE TARGET TYPE
    result2 = neo.run(cypher2)
    
    # EXTRACT UNIQUE GENES ASSOCIATED WITH CYPHER REQUESTS
    elements = set(i[0] for i in result2.to_table()) 
    common_elements = elements.intersection(query)  # common elements between query and target
    k = len(common_elements) # c/k:number of successes in the population = common between query and elements
    t = len(elements)
    p = len(elements)/n_gene_total # proba of success each essay

    # COVERAGE ####################################################
    if param.measure == 'coverage':
        proba = 1 - (k/query_size) * (k/t)

    # BINOMIAL TEST ####################################################
    if param.measure == 'binomial': # by default binom.cdf(<=success, attempts, proba). we need calculate p(X>=x)
        # based on success
        # proba = 1 - binom.cdf(len(common_elements)-1, query_size, p) #query_size - len(common_elements)
        # based on failure : best
        proba = binom.cdf(query_size - k, query_size, 1 - p)

    # HYPERGEOMETRIC TEST ##################################################
    if param.measure == 'hypergeometric':
        proba = 1 - hypergeom.cdf(k - 1, n_gene_total, t, query_size)
        # x = k, M = n_gene_total, t = len(elements), n = query_size pdf = hypergeom.pmf(x, M, n, N)
        # k success parmi n du query/target, N = population size, t = number of target in the population, N = number of query

    # CHI2 ########################################################
    if param.measure == 'chi2': # adjusted chi square cdf(x, df, loc=0, scale=1)
        # Creation of the contingency table as observed = [[c, q - c], [t - c, g - q - t + c]]
        observed = [[k, query_size - k], [ t - k, n_gene_total - ( t + query_size - k)]]
        # chi2_contingency test
        chi2_stat, proba, _, _ = chi2_contingency(observed)


    # if proba <= param.alpha:
    dico_proba = {
                  "id": set_id, 
                  "p-value": proba,
                  'common.n':len(common_elements), 
                  'target.n': len(elements),
                  'elements.target': elements, 
                  'elements.common': common_elements ,
                  'desc': desc
                }      
    liste_proba.append(dico_proba)
    
liste_triee = sorted(liste_proba, key=lambda x: x["p-value"])

if param.adjust:
    liste_ajust = []
    
    for i, r in enumerate(liste_triee):
        p_val_ajust = (r["p-value"] * query_size) / (i + 1) 
        if p_val_ajust <= param.alpha:
            liste_ajust.append(r)
            
    pprint(f'Number of significative adjusted p-value(0.05): {len(liste_ajust)}')
    for r in liste_ajust:
        print("{}\t{}\t{}/{}\t{}\t{}".format(r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))

else:
    pprint(f'Number of significative adjusted p-value(0.05): {len(liste_triee)}')
    for r in liste_triee:
        print("{}\t{}\t{}/{}\t{}\t{}".format(r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))
