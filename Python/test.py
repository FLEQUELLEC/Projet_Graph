#!/bin/env python
import gm
from pprint import pprint
#test = gm.read_delim('data/dressing.tsv')

#pprint(test)

#test BFS sur le noeud 'chemise' dans le graphe test
#pprint(gm.BFS(test, 'sous-vetements'))

prot_link = gm.graph.read_delim('data/511145.protein.links.experimental.txt', column_separator=' ')
#print(prot_link)

pprint(prot_link.BFS('511145.b0014', cible='511145.b3738', chemin=True))
