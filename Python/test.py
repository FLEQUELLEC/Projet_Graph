#!/bin/env python
import gm
from pprint import pprint
test = gm.graph.read_delim('Python/data/dressing.tsv', column_separator="\t")

#print(test)

#test BFS sur le noeud 'chemise' dans le graphe test
#pprint(gm.BFS(test, 'sous-vetements'))

#prot_link = gm.graph.read_delim('data/511145.protein.links.experimental.txt', column_separator=' ')
#print(prot_link)

#pprint(test.BFS('sous-vetements'))
#pprint(test.sousgraphe_induit('sous-vetement'))

bfs_result = test.BFS('sous-vetements')

visited_nodes = list(bfs_result['Distance'].keys())

subgraph = test.sousgraphe_induit(visited_nodes)
print(subgraph)