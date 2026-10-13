#!/usr/bin/env python

import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom #pip install scipy

# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets filename.')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')


# param = parser.parse_args(['-q', '/chemin/vers/votre/fichier/query.txt', '-t', '/chemin/vers/votre/fichier/ref.sets/uniprot.keywords.sets.json'])
param = parser.parse_args(['-q','../query/set1.txt','-t','../ref.sets/uniprot.keywords.sets.json','-c','-v','-l','10'])

# LOAD QUERY
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

  # LOAD REFERENCE SETS
sets = json.loads(open(param.sets).read())
# set est une liste de dictionnaires
if param.verbose:
    print('first target sets: ', sets[0:2])

len_query = len(query)
len_sets = len(sets) # nombre d ensembles de reference

# COMPUTE POPULATION SIZE, deduire la taille de la population dans ref set
population = set() # faire en sorte que population soit un ensemble
for s in sets:
  valeur_cle3 = s.get('elements')
  population.update(valeur_cle3)

# Maintenant, "valeurs_uniques" contient toutes les valeurs uniques de la troisième clé ('cle3')
n_gene_total = len(population)

# Affichez le nombre de valeurs uniques d elements dans la population
print(n_gene_total) # 4400

# EVALUATE SETS
results = []
query_size = len(query)
for s in sets:
    elements = set(s['elements']) # liste d elements unique dans elements
    common_elements = elements.intersection(query) #k, nombre de succes dans la population = commun entre query et elements
    k = len(common_elements) # nombre de succes dans la population = commun entre query et elements
    if param.measure=='binomial': # par défaut binom.cdf(<=success, attempts, proba). il nous faut calculer p(X>=x)
        p = 1 - binom.cdf(k - 1, len(elements), query_size / n_gene_total)
    else:
        print(f'sorry, {param.measure} not (yet) implemented')
        exit(1)
    r = { 'id': s['id'], 
         'desc': s['desc'], 
         'common.n':len(common_elements), 
         'target.n': len(elements), 
         'p-value': p, 
         'elements.target': elements, 
         'elements.common': common_elements }
    results.append( r )

if param.verbose:
  print(results)

significant_results = []
#count_significant = 0
# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item['p-value'])
for r in results:
  if r['p-value'] < param.alpha: 
    significant_results.append((r['p-value']))
    print("{}\t{}\t{}/{}\t{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], r['desc'], ', '.join(r['elements.common'])))

  if param.adjust:
    p_ajusted = []
    for i, proba in enumerate(significant_results):
      p = (proba*len_sets) / (i+1)
      if p < param.alpha :
        p_ajusted.append(p) 
        print(len(p_ajusted))