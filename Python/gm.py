#!/bin/env python
# -*- coding: utf-8 -*-
"""
Graph Manipulation Library (gm.py)
==================================
Auteurs : Florent LE QUELLEC

Cette bibliothèque fournit un ensemble de fonctions pour manipuler des graphes
représentés sous forme de dictionnaires Python. Compatible avec des graphes
dirigés ou non dirigés, pondérés ou non, et intégrable dans des workflows
bioinformatiques ou analytiques légers.
"""

import polars as pl
import pandas as pd
from pprint import pprint

class graph :
    def __init__(self, directed=True, weighted=False, weight_attribute=None):
        """
        Initialise un graphe vide.

        Parameters
        ----------
        directed : bool, optional
            Indique si le graphe est dirigé (True) ou non dirigé (False).
        weighted : bool, optional
            Indique si le graphe possède des poids d’arêtes.
        weight_attribute : str, optional
            Nom de l’attribut de poids si applicable.
        """
        self.nodes = {}
        self.edges = {}
        self.directed = directed
        self.weighted = weighted
        self.weight_attribute = weight_attribute


    def node_exists(self, n):
        """
        Vérifie si un nœud existe dans le graphe.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.
        n : str or int
            Identifiant du nœud à tester.

        Returns
        -------
        bool
            True si le nœud existe, False sinon.
        """
        return n in self.nodes


    def add_node(self, node_id, attributes=None):
        """
        Ajoute un nœud au graphe s’il n’existe pas déjà.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.
        node_id : str or int
            Identifiant du nœud à ajouter.
        attributes : dict, optional
            Dictionnaire d’attributs associés au nœud.

        Returns
        -------
        dict
            Dictionnaire des attributs du nœud ajouté (ou existant).
        """
        if node_id is not self.nodes:
            self.nodes[node_id] = attributes or {}
            self.edges[node_id] = {}  # initialise les arêtes sortantes
        return self.nodes[node_id]


    def edge_exists(self, n1, n2):
        """
        Vérifie si une arête entre deux nœuds existe dans le graphe.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.
        n1 : str or int
            Nœud source.
        n2 : str or int
            Nœud cible.

        Returns
        -------
        bool
            True si l’arête (n1, n2) existe, False sinon.
        """
        return n1 in self.edges and n2 in self.edges[n1]


    def add_edge(self, node_id1, node_id2, attributes=None):
        """
        Ajoute une arête entre deux nœuds. Crée les nœuds si nécessaire.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.
        node_id1 : str or int
            Identifiant du nœud source.
        node_id2 : str or int
            Identifiant du nœud cible.
        attributes : dict, optional
            Dictionnaire d’attributs de l’arête.

        Returns
        -------
        dict
            Dictionnaire des attributs de l’arête ajoutée.
        """
        if node_id1 is not self.nodes:
            self.add_node(node_id1)
        if node_id2 is not self.nodes:
            self.add_node(node_id2)

        if not  self.edge_exists(node_id1, node_id2):
            if self.edge_exists(node_id1,node_id2) :
                self.edges[node_id1][node_id2] = attributes or {}
            if not self.directed:
                self.edges[node_id2][node_id1] = self.edges[node_id1][node_id2]
        return self.edges[node_id1][node_id2]

    def nodes(self):
        """
        Renvoie la liste triée des identifiants de nœuds du graphe.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.

        Returns
        -------
        list
            Liste triée des nœuds du graphe.
        """
        return sorted(self.nodes.keys())


    def nb_nodes(self):
        """
        Renvoie le nombre de nœuds du graphe.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.

        Returns
        -------
        int
            Nombre de nœuds.
        """
        return len(self.nodes)


    def nb_edges(self):
        """
        Renvoie le nombre de sommets ayant au moins une arête sortante.

        (Attention : pour les graphes non dirigés, chaque arête apparaît deux fois.)

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.

        Returns
        -------
        int
            Nombre d’entrées dans la table des arêtes.
        """
        return sum(len(v) for v in self.edges.values())


    def neighbors(self, node_id):
        """
        Renvoie la liste des voisins d’un nœud donné.

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.
        node_id : str or int
            Identifiant du nœud.

        Returns
        -------
        list
            Liste des identifiants de nœuds voisins.
        """
        return list(self.edges[node_id].keys())


    def edges_tuples(self):
        """
        Renvoie la liste de toutes les arêtes sous forme de tuples (source, cible).

        Parameters
        ----------
        g : dict
            Graphe sous forme de dictionnaire.

        Returns
        -------
        list of tuple
            Liste de tuples (u, v) représentant les arêtes.
        """
        return [(u, v) for u in self.nodes for v in self.neighbors(u)]

    @classmethod
    def read_delim(cls, filename, column_separator='\t', directed=True, weighted=False, weight_attribute=None):
        """
        Lit un fichier délimité (ex: TSV, CSV) et construit un graphe.

        Les deux premières colonnes représentent les nœuds connectés,
        les suivantes contiennent les attributs d’arêtes.

        Parameters
        ----------
        filename : str
            Chemin vers le fichier à lire.
        column_separator : str, optional
            Caractère de séparation des colonnes (par défaut : tabulation).
        directed : bool, optional
            Indique si le graphe doit être dirigé.
        weighted : bool, optional
            Indique si le graphe doit être pondéré.
        weight_attribute : str, optional
            Nom de l’attribut de poids.

        Returns
        -------
        dict
            Graphe construit à partir du fichier.
        """
        df = pl.read_csv(
            filename,
            separator=column_separator,
            has_header=True,
            infer_schema_length=1000,
            quote_char=None
        )

        cols = df.columns
        if len(cols) < 2:
            raise ValueError("Le fichier doit contenir au moins deux colonnes (source, target).")

        src_col, tgt_col = cols[0], cols[1]
        att_cols = cols[2:]

        g = cls(directed, directed=True, weighted=False, weight_attribute=None)
        pdf = df.to_pandas()

        for _, row in pdf.iterrows():
            u = row[src_col]
            v = row[tgt_col]
            att = {col: row[col] for col in att_cols}
            g.add_edge(g, u, v, att)

        return g

    def BFS(self, s) :
        """
        l'algorithme de parcours en largeur (Breadth-First Search) a pour objectif de visiter tous les sommets d'un graphe G, afin de determine le chemin le plus court entre un sommet de depart s et tous les autres sommets du graphe.

        G est votre graphe sous forme de dictionnaire
        s est le sommet de depart
        non_visites est un booléen qui indique si on veut inclure les sommets non traités dans le résultat final

        La fonction retourne un dictionnaire contenant pour chaque sommet son état (blanc, gris, noir), sa distance par rapport au sommet de départ, et son parent dans le parcours.
        """
    # création variable des états des sommets : couleurs pour chaque sommet non visité (blanc), en cours de visite (gris), visité (noir), distances pour la distance entre le sommet de départ et chaque sommet, parents pour le parent de chaque sommet dans le parcours
        etat, distances, parents = {}, {}, {}

        etat[s] = 'gris' # initialisation du sommet de départ, etat gris, c'est à dire en cours de visite , distance 0, pas de parent
        distances[s] = 0
        attente = [s] # création de la file d'attente pour le parcours

        while len(attente) != 0: # tant que la file n'est pas vide, la boucle continue
            u = attente.pop(0) # extraction du premier sommet de la file, pour signifier qu'on le visite
            for voisin in self.edges[u]:# pour chaque voisin du sommet u
                if voisin not in etat: # si le voisin n'a pas encore été visité
                    etat[voisin] = 'gris' # on le marque comme en cours de visite (etat gris)
                    distances[voisin] = distances[u] + 1 # on met à jour la distance du voisin en fonction de la distance du sommet u
                    parents[voisin] = u # on met à jour le parent du voisin comme étant u
                    attente.append(voisin) # et on les rajoutes dans la file d'attente
            etat[u] = 'noir' #sommet visité, on le marque en noir
        return {"état" : etat, "Distance" : distances, "parents" : parents, "source" : s}



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

# Créer un graphe d’exemple
g = create_graph(directed=False)
add_edge(g, 'A', 'B')
add_edge(g, 'A', 'C')
add_edge(g, 'B', 'D')
add_edge(g, 'C', 'E')

print("Graphe :")
pprint(g)

# Lancer BFS
visited = BFS(g, 'A')
print("\nOrdre de visite BFS depuis 'A' :", visited)
