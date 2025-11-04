#!/bin/env python
import gm
from pprint import pprint
test = gm.read_delim('data/dressing.tsv')

pprint(gm.nb_nodes(test))

for u in test['edges']:
    for v in test['edges'][u]:
        print(f"edges ({u},{v})")

for u in gm.nodes(test):
    for v in gm.neighbors(test, u):
        print(f"edges ({u},{v})")
