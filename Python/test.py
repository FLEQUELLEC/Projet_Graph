#!/bin/env python
import gm
from pprint import pprint
#test = gm.read_delim('data/dressing.tsv')

#pprint(test)

#test BFS sur le noeud 'chemise' dans le graphe test
#pprint(gm.BFS(test, 'sous-vetements'))

prot_link = gm.read_delim('data/511145.protein.links.experimental.txt', column_separator=' ')
pprint(gm.BFS(prot_link, '511145.b0014'))
