#!/bin/env python

from pprint import pprint

def create_graph(directed = True, weighted = False, weight_attribute = None):
    """
    create a dictionnary representing a graph and returns it.
    """
    g = { 'nodes': {}, 'edges': {}, 'directed': directed, 'weighted': weighted, 'weight_attribute': weight_attribute }
    return g

def node_exists(g,n): return n in g['nodes']

def add_node(g, node_id, attributes = None): 
    """
    add a node with node_id (node id provided as a string or int) to the graph g.
    attributes on the node can be provided by a dict.
    returns the node n attributes.
    """
    if not node_exists(g, node_id): # ensure node does not already exist
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['nodes'][node_id] = attributes
        g['edges'][node_id] = {} # init outgoing edges
    return g['nodes'][node_id] # return node attributes

def edge_exists(g, n1, n2): return node_exists(g, n1) and n2 in g['edges'].get(n1, {})


def add_edge(g, node_id1, node_id2, attributes = None): 
    # create nodes if they do not exist
    if not node_exists(g, node_id1): add_node(g, node_id1)
    if not node_exists(g, node_id2): add_node(g, node_id2)
    # add edge(s) only if they do not exist
    if not edge_exists(g,node_id1,node_id2):
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['edges'][node_id1][node_id2] = attributes
        if not g['directed']:
            g['edges'][node_id2][node_id1] = g['edges'][node_id1][node_id2] # share the same attributes as n1->n2
    return g['edges'][node_id1][node_id2] # return edge attributes

def read_delim(filename, column_separator='\t', directed=True, weighted=False, weight_attribute=None): 
    """
    Parse a text file which columns are separated by the specified column_separator and 
    returns a graph.
    
    line syntax: node_id1   node_id2    att1    att2    att3    ...
    """
    g = create_graph(directed, weighted, weight_attribute)
    with open(filename) as f: 
        # GET COLUMNS NAMES
        tmp = f.readline().rstrip()
        attNames= tmp.split(column_separator)
        # REMOVES FIRST TWO COLUMNS WHICH CORRESPONDS TO THE LABELS OF THE CONNECTED VERTICES
        attNames.pop(0)  # remove first column name (source node not to be in attribute names)
        attNames.pop(0)  # remove second column (target node ...)
        # PROCESS THE REMAINING LINES
        row = f.readline().rstrip()
        while row:
            vals = row.split(column_separator)
            u = vals.pop(0)
            v = vals.pop(0)
            att = {}
            for i in range(len(attNames)):
                att[ attNames[i] ] = vals[i]
            add_edge(g, u, v, att)
            row = f.readline().rstrip() # NEXT LINE
        return g
    
def nodes(g) : return g['nodes']

def nb_nodes(g): return len(g.get('nodes'))

def nb_edges(g): return len(g.get('edges'))

def neighbors(g, node_id): return list(g['edges'][node_id].keys())

def directed(g) : return g('directed')

def edges_tuples (g) : return [(u,v) for u in nodes(g) for v in neighbors(g,u)]

##### main → tests #####
if __name__ == "__main__":
  print("# Graph lib tests")
  print("## create_graph")
  g = create_graph()
  pprint(g)
  
  print("## add nodes and edges")
  g = create_graph()
  add_node(g, 'A')
  add_node(g, 'B')
  add_edge(g, 'A', 'B', { 'weight': 5 } )
  pprint(g)