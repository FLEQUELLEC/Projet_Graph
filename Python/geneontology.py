#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gene Ontology Analysis Module (geneontology.py)
===============================================
Auteurs : Florent LE QUELLEC

Ce module permet de charger et de manipuler les données de la Gene Ontology (GO).
Il construit un graphe orienté acyclique (DAG) à partir des fichiers de définition
(.obo) et y associe les annotations géniques (.goa/.gaf).

Fonctionnalités principales :
- Parsing robuste des fichiers OBO (gestion des forward references).
- Parsing des fichiers d'annotations (GOA).
- Exploration de la hiérarchie (Ancêtres / Descendants).
- Calcul optimisé de la profondeur des ontologies.

Dépendances :
    - gm (Module de gestion de graphe local)
    - re (Expressions régulières)
"""

import re
import sys
import os
import gm

# --- Visualisation conceptuelle ---
# Le graphe GO est un DAG (Directed Acyclic Graph).
# Les arcs 'is_a' vont de l'Enfant vers le Parent (Spécifique -> Générique).
#

def load_OBO(filename='go-basic.obo'):
    """
    Parse un fichier OBO et construit un graphe de termes GO.

    Cette fonction gère les références anticipées ("forward references") :
    si un terme enfant fait référence à un parent qui n'a pas encore été
    lu dans le fichier, le parent est créé à la volée (noeud squelette)
    pour garantir la connectivité du graphe.

    Parameters
    ----------
    filename : str, optional
        Chemin vers le fichier .obo. Par défaut 'go-basic.obo'.

    Returns
    -------
    gm.graph
        Le graphe orienté contenant les termes GO.
        Structure des arêtes : Enfant -> Parent (relations 'is_a', 'part_of').
    """
    go_graph = gm.graph(directed=True, weighted=False)
    go_graph.alt_id = {} # Stockage des IDs alternatifs pour redirection

    def parseTerm(lines):
        """Fonction interne pour parser un bloc [Term] ligne par ligne."""
        go_id = None
        is_obsolete = False

        # 1. Extraction ID et Obsolescence
        for line in lines:
            if line.startswith('id:'):
                # Exemple : "id: GO:0000001" -> on garde "GO:0000001"
                go_id = line.split()[1]
            elif line.startswith('is_obsolete: true'):
                is_obsolete = True

        # On ignore les termes sans ID ou obsolètes
        if not go_id or is_obsolete:
            return

        # 2. Création ou Mise à jour du nœud
        # Note : Le nœud peut déjà exister (créé comme parent par un enfant lu avant)
        if go_id not in go_graph.nodes:
            go_graph.add_node(go_id, {'type': 'GOTerm'})
        else:
            # Si le nœud existait déjà, on confirme son type
            go_graph.nodes[go_id]['type'] = 'GOTerm'

        go_attr = go_graph.nodes[go_id]

        # 3. Parsing des attributs et relations
        for line in lines:
            # Nettoyage des commentaires inline (ex: " ... ! description")
            if ' ! ' in line:
                line = line.split(' ! ')[0].strip()

            if line.startswith('name:'):
                go_attr['name'] = line.replace('name:', '').strip()
            elif line.startswith('namespace:'):
                go_attr['namespace'] = line.replace('namespace:', '').strip()
            elif line.startswith('def:'):
                # Nettoyage des guillemets autour de la définition
                go_attr['def'] = line.replace('def:', '').strip().strip('"')
            elif line.startswith('alt_id:'):
                parts = line.split()
                if len(parts) > 1:
                    go_graph.alt_id[parts[1]] = go_id

            # --- GESTION DES RELATIONS ---
            parent_id = None
            relation = None

            if line.startswith('is_a:'):
                parts = line.split()
                if len(parts) >= 2:
                    parent_id = parts[1]
                    relation = 'is_a'

            elif line.startswith('relationship: part_of'):
                parts = line.split()
                # On cherche l'élément qui ressemble à un ID GO
                for p in parts:
                    if p.startswith('GO:'):
                        parent_id = p
                        relation = 'part_of'
                        break

            # Ajout de l'arête avec sécurité "Forward Reference"
            if parent_id:
                # Si le parent n'existe pas encore, on le crée (squelette)
                if parent_id not in go_graph.nodes:
                    go_graph.add_node(parent_id, {'type': 'GOTerm'})

                # Ajout de l'arête : Terme Courant (Enfant) -> Parent
                go_graph.add_edge(go_id, parent_id, {'relationship': relation})

    # Lecture principale du fichier
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            buff = []
            in_term = False
            for line in f:
                line = line.strip()
                if not line: continue

                if line == '[Term]':
                    # Si on a fini de lire un terme précédent, on le parse
                    if buff: parseTerm(buff)
                    buff = []
                    in_term = True
                elif line == '[Typedef]':
                    # On arrête le parsing si on tombe sur les définitions de types
                    if buff: parseTerm(buff)
                    break
                elif in_term:
                    buff.append(line)
            # Ne pas oublier le dernier buffer à la fin du fichier
            if buff and in_term: parseTerm(buff)

    except FileNotFoundError:
        print(f"Erreur critique: Le fichier {filename} est introuvable.")
        return gm.graph()

    return go_graph


def load_GOA(go, filename, warnings=True):
    """
    Parse un fichier GOA (Gene Ontology Annotation) et ajoute les produits géniques au graphe.

    Le format attendu est le format GAF (Gene Association File).
    Les arcs créés vont du Gène vers le Terme GO.

    Parameters
    ----------
    go : gm.graph
        Le graphe GO chargé précédemment (objet mutable).
    filename : str
        Chemin vers le fichier .gaf/.goa.
    warnings : bool, optional
        Afficher les avertissements si un terme annoté n'existe pas dans l'OBO.

    Returns
    -------
    None
        Le graphe est modifié en place.
    """
    # Optimisation : Cache local pour éviter les lookups répétés sur self.nodes
    nodes = go.nodes
    alt_ids = go.alt_id

    #

    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('!'): continue # Ignore les commentaires

                cols = line.rstrip().split('\t')
                # Sécurité format GAF (min 15 colonnes, ici on vérifie 11 pour le nécessaire)
                if len(cols) < 11: continue

                gp_id = cols[1]  # DB_Object_ID (Identifiant unique du gène)
                gt_id = cols[4]  # GO ID (Terme annoté)

                # 1. Résolution d'identifiant alternatif
                # (Si le fichier d'annotation utilise un vieil ID redirigé dans l'OBO)
                if gt_id not in nodes and gt_id in alt_ids:
                    gt_id = alt_ids[gt_id]

                # 2. Vérification existence du terme
                if gt_id not in nodes:
                    if warnings:
                        # Ce cas arrive si l'OBO est plus vieux que le GOA ou incomplet
                        # print(f"⚠️ Impossible de rattacher {gp_id} à {gt_id} (Terme inconnu)")
                        pass
                    continue

                # 3. Création du Produit Génique (s'il n'existe pas encore)
                if gp_id not in nodes:
                    go.add_node(gp_id, {
                        'id': gp_id,
                        'type': 'GeneProduct',
                        'name': cols[2],     # Symbol
                        'desc': cols[9],     # Full Name
                        'aliases': cols[10].split('|')
                    })

                # 4. Ajout de l'annotation (Arête GeneProduct -> GOTerm)
                e_attr = go.add_edge(gp_id, gt_id, {'relationship': 'annotation'})

                # Ajout du code de preuve (ex: IEA, EXP...) à la liste des preuves pour cette arête
                e_attr.setdefault('evidence-codes', []).append(cols[6])

    except FileNotFoundError:
        print(f"Erreur: Le fichier d'annotation {filename} est introuvable.")


def GOTerms(go, gp_id, recursive=False):
    """
    Retourne les termes GO associés à un produit génique.

    Parameters
    ----------
    go : gm.graph
        Le graphe contenant l'ontologie et les annotations.
    gp_id : str
        L'identifiant du gène (GeneProduct).
    recursive : bool, optional
        Si True, remonte l'arbre pour inclure tous les termes parents (Ancêtres).
        Si False, retourne uniquement les annotations directes.

    Returns
    -------
    list
        Liste des identifiants GO trouvés.
    """
    if gp_id not in go.nodes:
        return []

    # Les gènes pointent vers les termes (neighbors directs)
    termes_initiaux = go.neighbors(gp_id)

    if not recursive:
        return list(termes_initiaux)
    else:
        # Parcours en largeur (BFS) pour remonter vers les parents (relation is_a)
        # Note : Dans notre graphe, Edge = Enfant -> Parent, donc neighbors = Parents.
        res = set(termes_initiaux)
        file_attente = list(termes_initiaux)

        while len(file_attente) != 0:
            courant = file_attente.pop(0)
            parents = go.neighbors(courant)

            for parent in parents:
                if parent not in res:
                    res.add(parent)
                    file_attente.append(parent)
        return list(res)


def GeneProducts(go, go_id, recursive=False):
    """
    Retourne les produits géniques associés à un terme GO.

    Parameters
    ----------
    go : gm.graph
        Le graphe.
    go_id : str
        L'identifiant du terme GO.
    recursive : bool, optional
        Si True, descend dans l'arbre pour inclure les gènes annotés
        sur des termes enfants (Descendants / Plus spécifiques).

    Returns
    -------
    list
        Liste des identifiants de gènes trouvés.
    """
    if go_id not in go.nodes:
        return []

    termes_cible = {go_id}

    # 1. Si récursif : On descend dans l'arbre pour trouver tous les termes spécifiques
    if recursive:
        file = [go_id]
        ensemble_visites = {go_id}
        while len(file) != 0:
            term = file.pop()
            # On cherche qui pointe vers nous (Predecessors = Enfants)
            # Rappel : Edge = Enfant -> Parent
            entrants = go.predecessors(term)

            for entrant in entrants:
                # On ne descend que dans les termes GO (on ne saute pas directement aux gènes)
                # Utilisation de .get() pour sécurité (au cas où un nœud n'a pas de type)
                if go.nodes[entrant].get('type') == 'GOTerm' and entrant not in ensemble_visites:
                    termes_cible.add(entrant)
                    ensemble_visites.add(entrant)
                    file.append(entrant)

    # 2. Pour tous les termes identifiés (racine + descendants), on récupère les gènes
    # Les gènes sont des prédécesseurs des termes (Gene -> Terme)
    genes_trouves = set()
    for term in termes_cible:
        sources = go.predecessors(term)
        for source in sources:
            if go.nodes[source].get('type') == 'GeneProduct':
                genes_trouves.add(source)

    return list(genes_trouves)


def max_depth(go):
    """
    Calcule la profondeur maximale des 3 sous-ontologies (BP, MF, CC).

    L'algorithme calcule la "Hauteur" de la racine, définie comme la longueur
    du chemin le plus long de la racine vers une feuille (terme le plus spécifique).

    Optimisation :
    Utilise un index inversé pré-calculé pour transformer la recherche des enfants
    (qui nécessite un scan complet O(N) dans gm.py) en accès direct O(1).

    Parameters
    ----------
    go : gm.graph
        Le graphe chargé.

    Returns
    -------
    dict
        {'BP': int, 'MF': int, 'CC': int}
    """
    racines = {
        'BP': 'GO:0008150', # Biological Process
        'MF': 'GO:0003674', # Molecular Function
        'CC': 'GO:0005575'  # Cellular Component
    }

    resultats = {}
    memo = {}

    # --- OPTIMISATION : Construction de l'index inversé ---
    # Structure : { Parent : [Liste des Enfants] }
    # Permet de descendre dans l'arbre (Parent -> Enfant) instantanément.
    print("   [Optimisation] Construction de l'index inversé...", end='', flush=True)

    reverse_graph = {}
    # On parcourt toutes les arêtes. Dans OBO : Source=Enfant, Cible=Parent.
    for source_node in go.edges:
        targets = go.edges[source_node] # targets = parents

        for target_node in targets:
            if target_node not in reverse_graph:
                reverse_graph[target_node] = []
            reverse_graph[target_node].append(source_node) # target a pour enfant source

    print(" Fait.")
    # ------------------------------------------------------

    def get_height(u):
        """Fonction récursive avec mémoïsation (Programmation dynamique)"""
        if u in memo: return memo[u]

        # Récupération optimisée des enfants
        enfants = reverse_graph.get(u, [])

        # FILTRAGE : On ne descend que vers d'autres termes GO
        go_enfants = []
        for c in enfants:
            if c in go.nodes:
                # On exclut les gènes (qui sont aussi des "enfants" dans le sens où ils pointent vers le terme)
                # On accepte 'GOTerm' ou None (tolérance) mais pas 'GeneProduct'
                if go.nodes[c].get('type') != 'GeneProduct':
                    go_enfants.append(c)

        # Cas de base : Feuille (pas d'enfants termes GO)
        if not go_enfants:
            memo[u] = 0
            return 0

        # Récursion : Hauteur = 1 + max(Hauteur des enfants)
        try:
            h = 1 + max(get_height(child) for child in go_enfants)
        except ValueError:
            h = 0

        memo[u] = h
        return h

    # Lancement du calcul pour les 3 racines
    for nom, root_id in racines.items():
        if root_id in go.nodes:
            try:
                memo = {} # Reset du cache par sécurité entre les ontologies
                resultats[nom] = get_height(root_id)
            except RecursionError:
                print(f"Erreur: Profondeur excessive pour {nom} (Cycle possible ?)")
                resultats[nom] = 0
        else:
            resultats[nom] = 0

    return resultats


##### Tests unitaires basiques #####
if __name__ == "__main__":
    print("# Gene Ontology module tests")

    # Définition du chemin de test (à adapter selon votre arborescence)
    # Recherche locale ou dans un sous-dossier 'data'
    test_paths = [
        "go-basic.obo",
        "Python/data/projet/go-basic.obo",
        "data/go-basic.obo"
    ]

    filename = None
    for p in test_paths:
        if os.path.exists(p):
            filename = p
            break

    if filename:
        print(f"Chargement de {filename}...")
        go = load_OBO(filename)
        print(f"Graph chargé : {len(go.nodes)} noeuds")

        if len(go.nodes) > 0:
            print("\nCalcul des profondeurs (max_depth)...")
            depths = max_depth(go)
            print("Profondeurs :", depths)

            # Vérification de cohérence (CC ne doit pas être 0)
            if depths.get('CC', 0) == 0:
                print("⚠️  Avertissement : Profondeur CC nulle. Vérifiez les relations part_of.")
    else:
        print("⚠️  Fichier 'go-basic.obo' introuvable pour le test.")
        print("   Veuillez placer le fichier dans le dossier courant ou modifier le chemin.")
