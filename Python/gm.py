#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Graph Manipulation Library (gm.py)
==================================
Auteurs : Florent LE QUELLEC

Cette bibliothèque fournit un ensemble de fonctions pour manipuler des graphes
représentés sous forme de dictionnaires Python. Compatible avec des graphes
dirigés ou non dirigés, pondérés ou non, et intégrable dans des workflows
bioinformatiques ou analytiques légers.

Dépendances : polars, pandas
"""

import polars as pl
import pandas as pd
# pprint est utile pour le debug, on le garde
from pprint import pprint

class graph:
    """
    Classe représentant un graphe (orienté ou non, pondéré ou non).
    Structure de données : Liste d'adjacence via des dictionnaires imbriqués.
    """

    def __init__(self, directed=True, weighted=False, weight_attribute=None):
        """
        Initialise un graphe vide.

        Parameters
        ----------
        directed : bool, optional
            Indique si le graphe est dirigé (True) ou non dirigé (False). Par défaut True.
        weighted : bool, optional
            Indique si le graphe possède des poids d’arêtes. Par défaut False.
        weight_attribute : str, optional
            Nom de l’attribut de poids si applicable (ex: 'weight', 'distance').
        """
        self.nodes = {}  # Stocke les nœuds et leurs attributs : {id: {attr: val}}
        self.edges = {}  # Stocke la structure : {source: {cible: {attr: val}}}
        self.directed = directed
        self.weighted = weighted
        self.weight_attribute = weight_attribute

    def __str__(self):
        """Représentation textuelle du graphe."""
        lines = [
            f"Graphe {'dirigé' if self.directed else 'non dirigé'}",
            f"Nombre de nœuds : {len(self.nodes)}",
            f"Nombre d'arêtes : {self.nb_edges()}",
            "--- Aperçu des arêtes ---"
        ]
        # On limite l'affichage pour éviter de saturer la console sur les gros graphes
        count = 0
        for u, targets in self.edges.items():
            for v, attrs in targets.items():
                if count < 10:
                    lines.append(f"  {u} -> {v}  {attrs}")
                count += 1
        if count >= 10:
            lines.append(f"  ... (+ {count - 10} autres arêtes)")
        return "\n".join(lines)

    def node_exists(self, n):
        """
        Vérifie si un nœud existe dans le graphe.

        Parameters
        ----------
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
        node_id : str or int
            Identifiant du nœud à ajouter.
        attributes : dict, optional
            Dictionnaire d’attributs associés au nœud.

        Returns
        -------
        dict
            Dictionnaire des attributs du nœud ajouté (ou existant).
        """
        if node_id not in self.nodes:
            self.nodes[node_id] = attributes or {}
            self.edges[node_id] = {}  # initialise le dictionnaire des voisins
        return self.nodes[node_id]

    def edge_exists(self, n1, n2):
        """
        Vérifie si une arête entre deux nœuds existe.

        Parameters
        ----------
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
        Ajoute une arête entre deux nœuds. Crée les nœuds automatiquement s'ils manquent.

        Si le graphe est non dirigé, l'arête est ajoutée dans les deux sens
        (node_id1 -> node_id2 et node_id2 -> node_id1) partageant les mêmes attributs.

        Parameters
        ----------
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
        self.add_node(node_id1)
        self.add_node(node_id2)

        # On évite d'écraser une arête existante sauf si nécessaire
        if not self.edge_exists(node_id1, node_id2):
            attr = attributes or {}
            self.edges[node_id1][node_id2] = attr

            if not self.directed:
                # Partage de la référence du dictionnaire d'attributs
                self.edges[node_id2][node_id1] = self.edges[node_id1][node_id2]

        return self.edges[node_id1][node_id2]

    def get_nodes(self):
        """
        Renvoie la liste triée des identifiants de nœuds du graphe.
        (Renommé de 'nodes' pour éviter la confusion avec l'attribut self.nodes).

        Returns
        -------
        list
            Liste triée des clés des nœuds.
        """
        return sorted(self.nodes.keys())

    def nb_nodes(self):
        """
        Renvoie le nombre de nœuds du graphe.

        Returns
        -------
        int
            Nombre de nœuds.
        """
        return len(self.nodes)

    def nb_edges(self):
        """
        Renvoie le nombre d'arêtes du graphe.

        Pour les graphes non dirigés, divise le total par 2 car chaque arête
        est stockée deux fois (u->v et v->u).

        Returns
        -------
        int
            Nombre d'arêtes.
        """
        count = sum(len(v) for v in self.edges.values())
        return count // 2 if not self.directed else count

    def neighbors(self, node_id):
        """
        Renvoie la liste des voisins (successeurs) d’un nœud donné.

        Parameters
        ----------
        node_id : str or int
            Identifiant du nœud.

        Returns
        -------
        list
            Liste des identifiants de nœuds voisins.
        """
        if node_id in self.edges:
            return list(self.edges[node_id].keys())
        return []

    def predecessors(self, node_id):
        """
        Renvoie la liste des nœuds qui pointent vers node_id.
        (Parents dans un arbre, ou prédécesseurs dans un graphe orienté).

        Parameters
        ----------
        node_id : str or int
            Identifiant du nœud cible.

        Returns
        -------
        list
            Liste des prédécesseurs.
        """
        # Si le graphe n'est pas orienté, voisins = prédécesseurs
        if not self.directed:
            return self.neighbors(node_id)

        preds = []
        # On parcourt tous les nœuds sources possibles
        for u in self.edges:
            if node_id in self.edges[u]:
                preds.append(u)
        return preds

    def edges_tuples(self):
        """
        Renvoie la liste de toutes les arêtes sous forme de tuples (source, cible).

        Returns
        -------
        list of tuple
            Liste de tuples (u, v).
        """
        return [(u, v) for u in self.nodes for v in self.neighbors(u)]

    @classmethod
    def read_delim(cls, filename, column_separator='\t', directed=True, weighted=False, weight_attribute=None):
        """
        Lit un fichier délimité (ex: TSV, CSV) et construit un objet graph.

        Les deux premières colonnes représentent les nœuds connectés (Source, Target),
        les suivantes contiennent les attributs d’arêtes.

        Parameters
        ----------
        filename : str
            Chemin vers le fichier à lire.
        column_separator : str, optional
            Caractère de séparation (par défaut : tabulation).
        directed : bool, optional
            Indique si le graphe doit être dirigé.
        weighted : bool, optional
            Indique si le graphe doit être pondéré.
        weight_attribute : str, optional
            Nom de l’attribut de poids.

        Returns
        -------
        graph
            Instance de la classe graph remplie.
        """
        # Lecture optimisée avec Polars
        try:
            df = pl.read_csv(
                filename,
                separator=column_separator,
                has_header=True,
                infer_schema_length=1000,
                quote_char=None,
                truncate_ragged_lines=True
            )
        except Exception as e:
            print(f"Erreur lors de la lecture du fichier : {e}")
            return cls()

        cols = df.columns
        if len(cols) < 2:
            raise ValueError("Le fichier doit contenir au moins deux colonnes (source, target).")

        src_col, tgt_col = cols[0], cols[1]
        att_cols = cols[2:]

        # Création de l'instance
        g = cls(directed=directed, weighted=weighted, weight_attribute=weight_attribute)

        # Conversion en pandas pour l'itération (souvent plus simple pour les lignes mixtes)
        pdf = df.to_pandas()

        for _, row in pdf.iterrows():
            u = row[src_col]
            v = row[tgt_col]
            # Création du dictionnaire d'attributs pour l'arête
            att = {col: row[col] for col in att_cols}
            g.add_edge(u, v, att)

        return g

    def BFS(self, s, cible=None):
        """
        Algorithme de parcours en largeur (Breadth-First Search).
        Calcule les plus courts chemins (en nombre d'arêtes) depuis une source.

        Parameters
        ----------
        s : str or int
            Sommet de départ.
        cible : str or int, optional
            Sommet d'arrivée. Si spécifié, l'algorithme s'arrête dès qu'il est trouvé.

        Returns
        -------
        dict
            Si cible est trouvée : {'Distance': int, 'chemin': list, 'source': s}
            Sinon (parcours complet) : {'etats': dict, 'distances': dict, 'parents': dict, 'source': s}
        """
        if s not in self.nodes:
            print(f"Erreur : Le nœud de départ {s} n'existe pas.")
            return {}

        # Initialisation
        etats = {u: 'blanc' for u in self.nodes}
        distances = {u: float('inf') for u in self.nodes}
        parents = {u: None for u in self.nodes}

        etats[s] = 'gris'  # En cours de visite
        distances[s] = 0
        attente = [s]      # File (FIFO)

        while attente:
            u = attente.pop(0)

            if cible and u == cible:
                break

            for voisin in self.neighbors(u):
                if etats[voisin] == 'blanc':
                    etats[voisin] = 'gris'
                    distances[voisin] = distances[u] + 1
                    parents[voisin] = u
                    attente.append(voisin)

            etats[u] = 'noir'  # Visité

        # Reconstruction du chemin si une cible était demandée
        if cible:
            if parents[cible] is not None or cible == s:
                chemin = [cible]
                curr = cible
                while curr != s:
                    curr = parents[curr]
                    chemin.append(curr)
                chemin.reverse()
                return {"Distance": distances[cible], "chemin": chemin, "source": s}
            else:
                return {"Distance": float('inf'), "chemin": [], "source": s, "error": "Cible inatteignable"}

        return {"etats": etats, "distances": distances, "parents": parents, "source": s}

    def connected_components(self):
        """
        Identifie les composantes connexes d’un graphe NON orienté via BFS.

        Returns
        -------
        dict
            Dictionnaire {node_id: component_id}.
        """
        if self.directed:
            print("Attention : connected_components est conçu pour les graphes non orientés.")
            # On continue quand même, mais le résultat représente des composantes faiblement connexes
            # si on ignore la direction, ou juste l'accessibilité si on la garde.

        n_CC = 0
        CC = {u: None for u in self.nodes}

        for u in self.nodes:
            if CC[u] is None:
                # On lance un BFS complet depuis ce nœud non visité
                parcours = self.BFS(u)

                # Tous les nœuds atteints appartiennent à la composante n_CC
                # Le BFS retourne 'distances' pour les nœuds visités
                visited_nodes = [n for n, d in parcours['distances'].items() if d != float('inf')]

                for v in visited_nodes:
                    CC[v] = n_CC
                n_CC += 1
        return CC

    def sousgraphe_induit(self, nodes_subset):
        """
        Crée un sous-graphe induit par une liste de nœuds donnée.

        Parameters
        ----------
        nodes_subset : list
            Liste des identifiants de nœuds à conserver.

        Returns
        -------
        graph
            Nouveau graphe contenant uniquement les nœuds spécifiés et les arêtes entre eux.
        """
        sg = graph(
            directed=self.directed,
            weighted=self.weighted,
            weight_attribute=self.weight_attribute
        )

        # Ajout des nœuds
        for u in nodes_subset:
            if u in self.nodes:
                sg.add_node(u, self.nodes[u])

        # Ajout des arêtes existantes entre ces nœuds
        for u in nodes_subset:
            if u in self.edges:
                for v, attrs in self.edges[u].items():
                    if v in nodes_subset:
                        sg.add_edge(u, v, attrs)

        return sg

    def edges_filter(self, attribut, seuil):
        """
        Filtre le graphe pour ne garder que les arêtes respectant un critère numérique.

        Parameters
        ----------
        attribut : str
            Nom de l'attribut de l'arête à tester (ex: 'weight').
        seuil : float
            Valeur minimale pour conserver l'arête.

        Returns
        -------
        graph
            Nouveau graphe filtré.
        """
        G_filtre = graph(directed=self.directed, weighted=self.weighted)

        # On copie tous les nœuds
        for node_id, attrs in self.nodes.items():
            G_filtre.add_node(node_id, attrs)

        # On ne copie que les arêtes valides
        for u in self.edges:
            for v, attrs in self.edges[u].items():
                if attribut in attrs and attrs[attribut] >= seuil:
                    G_filtre.add_edge(u, v, attrs)
        return G_filtre

    def clustering_coefficient(self):
        """
        Calcule le coefficient de clustering local pour chaque nœud.
        C = (nombre de liens entre voisins) / (nombre de liens possibles entre voisins).

        Returns
        -------
        dict
            {node_id: coefficient (float entre 0 et 1)}
        """
        coeffs = {}
        for node_id in self.nodes:
            voisins = self.neighbors(node_id)
            k = len(voisins)

            # Si moins de 2 voisins, pas de connexions possibles entre eux -> coeff = 0
            if k < 2:
                coeffs[node_id] = 0.0
                continue

            liens_reels = 0
            # On regarde toutes les paires de voisins uniques
            for v1 in voisins:
                for v2 in voisins:
                    if v1 != v2:
                        if self.edge_exists(v1, v2):
                            liens_reels += 1

            # Nombre max de liens possibles entre k voisins
            # k * (k-1) pour dirigé, k * (k-1) / 2 pour non dirigé.
            # Ici l'algo compte v1->v2 et v2->v1 séparément dans la boucle,
            # donc on divise par k*(k-1) dans tous les cas pour normaliser.
            coeffs[node_id] = liens_reels / (k * (k - 1))

        return coeffs

    def DFS(self):
        """
        Algorithme de parcours en profondeur (Depth-First Search).
        Utilisé pour la classification des arêtes, la détection de cycles et le tri topologique.

        Returns
        -------
        dict
            Structure contenant l'état du parcours, les temps de découverte/fin,
            la classification des arêtes et l'ordre de fin.
        """
        a = {
            'etat': {},
            'parents': {},
            'decouvert': {}, # Timestamp de début
            'fin': {},       # Timestamp de fin
            'classification': {}, # Types d'arêtes
            'ordre_fin': [], # Liste pour tri topologique
            'temps': 0
        }

        # Initialisation
        for u in self.nodes:
            a['etat'][u] = 'inexplore'
            a['parents'][u] = None

        # Boucle principale pour traiter les composantes disjointes
        for u in self.nodes:
            if a['etat'][u] == 'inexplore':
                self._DFSvisite(u, a)
        return a

    def _DFSvisite(self, u, a):
        """Fonction interne récursive pour le DFS."""
        a['etat'][u] = 'decouvert'
        a['temps'] += 1
        a['decouvert'][u] = a['temps']

        for v in self.neighbors(u):
            if a['etat'][v] == 'inexplore':
                a['parents'][v] = u
                a['classification'][(u, v)] = 'branche'
                self._DFSvisite(v, a)
            elif a['etat'][v] == 'decouvert':
                # Arête vers un ancêtre non terminé -> Cycle détecté
                a['classification'][(u, v)] = 'retour'
            elif a['decouvert'][u] < a['decouvert'][v]:
                a['classification'][(u, v)] = 'arete avant'
            else:
                a['classification'][(u, v)] = 'transversale' # Cross edge

        a['etat'][u] = 'traite'
        a['temps'] += 1
        a['fin'][u] = a['temps']

        # Ajout à la liste pour le tri topologique (post-ordre)
        a['ordre_fin'].append(u)
        return a

    def is_acyclic(self):
        """
        Vérifie si le graphe est acyclique (DAG).

        Returns
        -------
        bool
            True si acyclique, False si contient un cycle.
        """
        # Un cycle existe ssi le DFS trouve une arête 'retour'
        res = self.DFS()
        for type_arete in res['classification'].values():
            if type_arete == 'retour':
                return False
        return True

    def topological_sort(self):
        """
        Effectue un tri topologique des nœuds (uniquement pour les DAG).

        Returns
        -------
        list or None
            Liste des nœuds ordonnée linéairement, ou None si cycle détecté.
        """
        if not self.is_acyclic():
            print("Erreur : Le graphe contient un cycle, tri topologique impossible.")
            return None

        res = self.DFS()
        # L'ordre topologique est l'inverse de l'ordre de fin de traitement
        return res['ordre_fin'][::-1]


##### Tests unitaires basiques #####
if __name__ == "__main__":
    print("# Graph lib tests")

    print("\n## 1. Création et affichage")
    g = graph(directed=True)
    g.add_node('A')
    g.add_node('B')
    g.add_edge('A', 'B', {'weight': 5})
    print(g)

    print("\n## 2. Graphe non dirigé")
    g2 = graph(directed=False)
    g2.add_edge('A', 'B')
    g2.add_edge('A', 'C')
    g2.add_edge('B', 'D')
    g2.add_edge('C', 'E')
    # Création d'un cycle pour tester is_acyclic plus tard (A-B-D-A impossible direct mais A-B-A oui en non dirigé)
    # Pour le test clustering
    g2.add_edge('B', 'C')

    print(g2)
    print(f"Nombre de noeuds: {g2.nb_nodes()}")
    print(f"Nombre d'arêtes: {g2.nb_edges()}")

    print("\n## 3. BFS (Plus court chemin)")
    # Chemin A -> E
    res_bfs = g2.BFS('A', 'E')
    print("Chemin A -> E :", res_bfs.get('chemin'))
    print("Distance :", res_bfs.get('Distance'))

    print("\n## 4. Clustering Coefficient")
    cc = g2.clustering_coefficient()
    print("Coefficients de clustering :", cc)

    print("\n## 5. DFS et Tri Topologique (Sur un DAG)")
    dag = graph(directed=True)
    dag.add_edge('Chemise', 'Ceinture')
    dag.add_edge('Chemise', 'Cravate')
    dag.add_edge('Cravate', 'Veste')
    dag.add_edge('Ceinture', 'Veste')
    dag.add_edge('Pantalon', 'Ceinture')
    dag.add_edge('Pantalon', 'Chaussures')
    dag.add_edge('Chaussettes', 'Chaussures')

    print("Est acyclique ?", dag.is_acyclic())
    print("Tri topologique (Ordre d'habillage) :", dag.topological_sort())

    print("\n## 6. Détection de cycle")
    cycle_g = graph(directed=True)
    cycle_g.add_edge('A', 'B')
    cycle_g.add_edge('B', 'C')
    cycle_g.add_edge('C', 'A') # Cycle
    print("Graphe cyclique est acyclique ?", cycle_g.is_acyclic())
