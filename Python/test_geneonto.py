#!/bin/env python

import geneontology as gom
import gm
from pprint import pprint

print('Loading OBO GeneOntology virion_component subgraph')
vir = gom.load_OBO('Python/data/go-virion_component.obo')

print('loading GOA annotation file for SARS-Cov2')
test = gom.load_GOA(vir, 'Python/data/uniprot_sars-cov-2.gaf', warnings=False)

print(vir)
print(len(vir.nodes))
pprint(vir.edges)

ia = 0
po = 0
for u, attr in vir.edges.items():
  for v, i in attr.items() :
    if i['relationship'] == 'is_a':
      ia += 1
    elif i['relationship'] == 'part_of':
      po +=1
print(ia, po) #45 is a et 20 part of apres est ce que c'est bon tel est la question !

pprint(vir.nodes.values())

cc, bp, mf = 0, 0, 0
for attr in vir.nodes.values():
  namespace = attr.get('namespace')
  if namespace == 'cellular_component':
    cc += 1
  elif namespace == 'biological_process':
    bp += 1
  elif namespace == 'molecular_function':
    mf += 1

print(cc, bp, mf)
