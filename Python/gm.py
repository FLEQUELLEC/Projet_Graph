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

def nodes(g) : return sorted(g['nodes'].keys())

def nb_nodes(g): return len(g.get('nodes'))

def nb_edges(g): return len(g.get('edges'))

def neighbors(g, node_id): return list(g['edges'][node_id].keys())

def directed(g) : return g('directed')

def edges_tuples (g) : return [(u,v) for u in nodes(g) for v in neighbors(g,u)]

def BFS(g, s, non_visites=True) :
    """
    l'algorithme de parcours en largeur (Breadth-First Search) a pour objectif de visiter tous les sommets d'un graphe G, afin de determine le chemin le plus court entre un sommet de depart s et tous les autres sommets du graphe.

    G est votre graphe sous forme de dictionnaire
    s est le sommet de depart
    non_visites est un booléen qui indique si on veut inclure les sommets non traités dans le résultat final

    La fonction retourne un dictionnaire contenant pour chaque sommet son état (blanc, gris, noir), sa distance par rapport au sommet de départ, et son parent dans le parcours.
    """
    # création variable des états des sommets : couleurs pour chaque sommet non visité (blanc), en cours de visite (gris), visité (noir), distances pour la distance entre le sommet de départ et chaque sommet, parents pour le parent de chaque sommet dans le parcours
    etat = dict()
    distances = dict()
    parents = dict()

    for u in g['nodes']: # creation des données intiales pour chaque sommet avant parcours
        etat[u] = 'blanc'
        distances[u] = float('inf')
        parents[u] = None

    etat[s] = 'gris' # initialisation du sommet de départ, etat gris, c'est à dire en cours de visite , distance 0, pas de parent
    distances[s] = 0
    attente = [s] # création de la file d'attente pour le parcours
    while len(attente) != 0: # tant que la file n'est pas vide, la boucle continue
        u = attente.pop(0) # extraction du premier sommet de la file, pour signifier qu'on le visite
        for voisin in g['edges'].get(u, {}):# pour chaque voisin du sommet u
            if etat[voisin] == 'blanc': # si le voisin n'a pas été visité (couleur blanc)
                etat[voisin] = 'gris' # on le marque comme en cours de visite (etat gris)
                distances[voisin] = distances[u] + 1 # on met à jour la distance du voisin en fonction de la distance du sommet u
                parents[voisin] = u # on met à jour le parent du voisin comme étant u
                attente.append(voisin) # et on les rajoutes dans la file d'attente
        etat[u] = 'noir' #sommet visité, on le marque en noir
    result = {"état" : etat, "Distance" : distances, "parents" : parents, "source" : s} #sauvegarde des résultats dans un dictionnaire
    if not non_visites:
        visited = {node for node, state in etat.items() if state != 'blanc'}
        result["état"] = {node: etat[node] for node in visited}
        result["Distance"] = {node: distances[node] for node in visited}
        result["parents"] = {node: parents[node] for node in visited}
    return result



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
