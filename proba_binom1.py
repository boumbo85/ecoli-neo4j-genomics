#!/usr/bin/env python

# UPLOAD PACKAGES
import argparse
from os.path import isfile
import json
import pprint
from scipy.stats import binom, hypergeom # pip install scipy
  # conda install python3 biopython requests


# Path Query
path = '../query/set1.txt'


# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets filename.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')


# param = parser.parse_args(['-q', '/chemin/vers/votre/fichier/query.txt', '-t', '/chemin/vers/votre/fichier/ref.sets/uniprot.keywords.sets.json'])
param = parser.parse_args()


# LOAD DATAS #######################################################
## LOAD QUERY
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


## LOAD TARGET/REFERENCE SETS
sets = json.loads(open(param.sets).read())
# set is a list of dictionaries
if param.verbose:
    print('First target sets: ', sets[0:2])
###########################################################################




# STATISTICAL TESTS #######################################################
###### EXTRACT LENGTH OF QUERY AND REFERENCE SETS #########################
query_size = len(query)
len_sets = len(sets) # number of sets in ref.sets
print(f"Number of genes in the query : {query_size}")

###### POPULATION SIZE ####################################################
population = set() # make sure that elements are uniques
for s in sets:
  valeur_cle3 = s.get('elements')
  population.update(valeur_cle3)

# EXTRACT LENGTH OF POPULATION
n_gene_total = len(population)
print(f"Total number of genes: {n_gene_total}") # Expected 4400 genes for uniprot.keywords.sets.json



###### COMPARE QUERY VS TARGET SETS #######################################################
results = [] # list of dictionaries

for s in sets:
  elements = set(s['elements']) # unique elements in the set
  common_elements = elements.intersection(query) # common elements between query and elements
  k = len(common_elements) # k:number of successes in the population = common between query and elements
  
  # BINOMIAL TEST ####################################################
  if param.measure == 'binomial': # by default binom.cdf(<=success, attempts, proba). we need calculate p(X>=x)
    p = 1 - binom.cdf(k - 1, len(elements), query_size / n_gene_total)
  
  # COVERAGE ####################################################
  if param.measure == 'coverage':
    p = len(sets[5]['elements']) / n_gene_total

  # HYPERGEOMETRIQUE ############################################
  if param.measure == 'hypergeometric':
    p = 1 - hypergeom.cdf(k - 1, n_gene_total, len(elements), query_size)
     
  # CHI2 ########################################################
  #if param.measure == 'chi2':
   # p = 1 - chi2.cdf(k - 1, len(elements), query_size / n_gene_total)
  
  
  #else:
     # print(f'sorry, {param.measure} not (yet) implemented')
      #exit(1)
  r = { 
        'id': s['id'], 
        'desc': s['desc'], 
        'common.n':len(common_elements), 
        'target.n': len(elements), 
        'p-value': p, 
        'elements.target': elements, 
        'elements.common': common_elements 
        }
  results.append( r )

if param.verbose:
  pprint(results)
  pprint(f"Number of results: {len(results)}")

## PRINT SIGNIFICANT RESULTS
significant_results = [] # create a list of significant results
results.sort(key=lambda an_item: an_item['p-value']) # sort results by p-value
for r in results:
  if r['p-value'] <= param.alpha: 
    significant_results.append((r['p-value']))
    # OUTPUT
    print("{}\t{}\t{}/{}\t{}\t{}".format(r['id'], 
                                        r['p-value'], 
                                        r['common.n'], 
                                        r['target.n'], 
                                        r['desc'], ', '.join(r['elements.common'])))
  
  
  # FDR ADJUSTMENT with p value adjusted  
  if param.adjust:
    p_ajusted = []
    for i, proba in enumerate(significant_results):
      p = (proba*len_sets) / (i+1)
      if p <= param.alpha :
        p_ajusted.append(p) 

# OUTPUT p
if param.adjust:
  print(f"Number of significant results: {len(p_ajusted)}")  
else:
  print(f"Number of significant results: {len(significant_results)}")


