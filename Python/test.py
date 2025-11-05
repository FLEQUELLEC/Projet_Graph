#!/bin/env python
import gm
from pprint import pprint
test = gm.read_delim('data/dressing.tsv')

pprint(test)

#test BFS sur le noeud 'chemise' dans le graphe test
pprint(gm.BFS(test, 'sous-vetements'))
