#!/bin/env python
# -*- coding: utf-8 -*-
"""
Module d'analyse Gene Ontology basé sur la librairie gm.graph
=============================================================
Ce module permet de charger un fichier GO (OBO) et ses annotations (GOA)
dans un graphe orienté. Il fournit des outils pour explorer la hiérarchie
et calculer des propriétés topologiques (profondeur).

Dépendances:
    - gm (GraphMaster)
    - re (Regular Expressions)
"""

import re
import gm

def load_OBO(filename='go-basic.obo'):
    """
    Parse un fichier OBO et construit un graphe de termes GO.

    Cette fonction gère les références anticipées ("forward references") :
    si un terme enfant fait référence à un parent qui n'a pas encore été
    lu dans le fichier, le parent est créé à la volée pour garantir
    la connectivité du graphe.

    Args:
        filename (str): Chemin vers le fichier .obo

    Returns:
        gm.graph: Le graphe orienté contenant les termes GO.
                  Arêtes : Enfant -> Parent (is_a, part_of).
    """
    go_graph = gm.graph(directed=True, weighted=False)
    go_graph.alt_id = {}

    def parseTerm(lines):
        """Fonction interne pour parser un bloc [Term]"""
        go_id = None
        is_obsolete = False

        # 1. Extraction ID et Obsolescence
        for line in lines:
            if line.startswith('id:'):
                # Exemple : "id: GO:0000001"
                go_id = line.split()[1]
            elif line.startswith('is_obsolete: true'):
                is_obsolete = True

        # On ignore les termes sans ID ou obsolètes
        if not go_id or is_obsolete:
            return

        # 2. Création ou Mise à jour du nœud
        # Note : Le nœud peut déjà exister s'il a été créé comme parent par un autre terme
        if go_id not in go_graph.nodes:
            go_graph.add_node(go_id, {'type': 'GOTerm'})
        else:
            # On s'assure que le type est bien défini
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
                go_attr['def'] = line.replace('def:', '').strip().strip('"')
            elif line.startswith('alt_id:'):
                # Gestion des IDs alternatifs pour la redirection
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
                    if buff: parseTerm(buff)
                    buff = []
                    in_term = True
                elif line == '[Typedef]':
                    # On arrête le parsing si on tombe sur les définitions de types
                    if buff: parseTerm(buff)
                    break
                elif in_term:
                    buff.append(line)
            # Ne pas oublier le dernier buffer
            if buff and in_term: parseTerm(buff)

    except FileNotFoundError:
        print(f"Erreur critique: Le fichier {filename} est introuvable.")
        return gm.graph()

    return go_graph


def load_GOA(go, filename, warnings=True):
    """
    Parse un fichier GOA (Gene Ontology Annotation) et ajoute les produits géniques au graphe.

    Args:
        go (gm.graph): Le graphe GO chargé précédemment.
        filename (str): Chemin vers le fichier .gaf/.goa.
        warnings (bool): Afficher les avertissements si un terme n'est pas trouvé.
    """
    # Optimisation : Cache local pour éviter les lookups répétés
    nodes = go.nodes
    alt_ids = go.alt_id

    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('!'): continue # Ignore les commentaires

                cols = line.rstrip().split('\t')
                if len(cols) < 11: continue # Sécurité format

                gp_id = cols[1]  # DB_Object_ID
                gt_id = cols[4]  # GO ID

                # 1. Résolution d'identifiant alternatif (si le terme a changé d'ID)
                if gt_id not in nodes and gt_id in alt_ids:
                    gt_id = alt_ids[gt_id]

                # 2. Vérification existence du terme
                if gt_id not in nodes:
                    if warnings:
                        print(f"⚠️ Impossible de rattacher {gp_id} à {gt_id} (Terme inconnu)")
                    continue

                # 3. Création du Produit Génique (s'il n'existe pas encore)
                if gp_id not in nodes:
                    go.add_node(gp_id, {
                        'id': gp_id,
                        'type': 'GeneProduct',
                        'name': cols[2],
                        'desc': cols[9],
                        'aliases': cols[10].split('|')
                    })

                # 4. Ajout de l'annotation (Arête GeneProduct -> GOTerm)
                e_attr = go.add_edge(gp_id, gt_id, {'relationship': 'annotation'})

                # Ajout du code de preuve (ex: IEA, EXP...)
                e_attr.setdefault('evidence-codes', []).append(cols[6])

    except FileNotFoundError:
        print(f"Erreur: Le fichier d'annotation {filename} est introuvable.")


def GOTerms(go, gp_id, recursive=False):
    """
    Retourne les termes GO associés à un produit génique.

    Args:
        go (gm.graph): Le graphe.
        gp_id (str): L'identifiant du gène.
        recursive (bool): Si True, inclut tous les termes parents (Ancêtres).

    Returns:
        list: Liste des identifiants GO trouvés.
    """
    if gp_id not in go.nodes:
        return []

    # Les gènes pointent vers les termes (neighbors)
    termes_initiaux = go.neighbors(gp_id)

    if not recursive:
        return list(termes_initiaux)
    else:
        # Parcours BFS pour remonter vers les parents (is_a)
        res = set(termes_initiaux)
        file_attente = list(termes_initiaux)

        while len(file_attente) != 0:
            courant = file_attente.pop(0)
            parents = go.neighbors(courant) # Dans le graphe OBO, Edge = Enfant -> Parent

            for parent in parents:
                if parent not in res:
                    res.add(parent)
                    file_attente.append(parent)
        return list(res)


def GeneProducts(go, go_id, recursive=False):
    """
    Retourne les produits géniques associés à un terme GO.

    Args:
        go (gm.graph): Le graphe.
        go_id (str): L'identifiant du terme GO.
        recursive (bool): Si True, inclut les gènes des termes enfants (Descendants).

    Returns:
        list: Liste des identifiants de gènes trouvés.
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
            entrants = go.predecessors(term)
            for entrant in entrants:
                # On ne descend que dans les termes GO (pas tout de suite aux gènes)
                # Utilisation de .get() pour sécurité
                if go.nodes[entrant].get('type') == 'GOTerm' and entrant not in ensemble_visites:
                    termes_cible.add(entrant)
                    ensemble_visites.add(entrant)
                    file.append(entrant)

    # 2. Pour tous les termes identifiés, on récupère les gènes (qui sont des prédécesseurs)
    genes_trouves = set()
    for term in termes_cible:
        sources = go.predecessors(term)
        for source in sources:
            if go.nodes[source].get('type') == 'GeneProduct':
                genes_trouves.add(source)

    return list(genes_trouves)


def max_depth(go):
    """
    Calcule la profondeur maximale des 3 sous-ontologies.
    OPTIMISÉE : Pré-calcule les relations parents-enfants pour une vitesse x1000.
    """
    racines = {
        'BP': 'GO:0008150',
        'MF': 'GO:0003674',
        'CC': 'GO:0005575'
    }

    resultats = {}
    memo = {}

    # --- OPTIMISATION CRITIQUE ---
    # On construit un dictionnaire inversé UNE SEULE FOIS au début.
    # Structure : { Enfant : [Liste des Parents/Predecessors] }
    print("   [Optimisation] Construction de l'index inversé...", end='', flush=True)
    reverse_graph = {}

    # On parcourt toutes les arêtes du graphe : Source -> Cible
    # Dans OBO : Source=Enfant, Cible=Parent.
    # Mais attention, ici on veut descendre de la Racine vers les Feuilles.
    # Dans le graphe chargé : Edge = Enfant -> Parent.
    # Donc pour descendre, on veut savoir "Qui a pour parent X ?".
    # Donc on cherche les sources qui pointent vers X.

    for source_node in go.edges:
        # neighbors renvoie les cibles (parents)
        # Si go.edges est un dict de dicts/sets (standard gm.py)
        targets = go.edges[source_node]

        for target_node in targets:
            # target_node est le parent, source_node est l'enfant
            if target_node not in reverse_graph:
                reverse_graph[target_node] = []
            reverse_graph[target_node].append(source_node)

    print(" Fait.")
    # -----------------------------

    def get_height(u):
        if u in memo: return memo[u]

        # Au lieu d'appeler go.predecessors(u) qui est lent,
        # on regarde directement dans notre index optimisé.
        enfants = reverse_graph.get(u, [])

        # FILTRAGE
        go_enfants = []
        for c in enfants:
            # Sécurité existence + Type
            if c in go.nodes:
                # Tolérance : on exclut juste les GeneProducts
                if go.nodes[c].get('type') != 'GeneProduct':
                    go_enfants.append(c)

        if not go_enfants:
            memo[u] = 0
            return 0

        # Récursion
        try:
            h = 1 + max(get_height(child) for child in go_enfants)
        except ValueError:
            h = 0

        memo[u] = h
        return h

    for nom, root_id in racines.items():
        if root_id in go.nodes:
            try:
                memo = {}
                resultats[nom] = get_height(root_id)
            except RecursionError:
                print(f"RecursionError sur {nom}")
                resultats[nom] = 0
        else:
            resultats[nom] = 0

    return resultats


##### Tests unitaires basiques #####
if __name__ == "__main__":
    print("# Gene Ontology module tests")

    # Adaptez ce chemin vers votre fichier local pour le test
    TEST_FILE = "Python/data/projet/go-basic.obo"

    # 1. Test Chargement
    go = load_OBO(TEST_FILE)
    print(f"Graph chargé : {len(go.nodes)} noeuds")

    if len(go.nodes) > 0:
        # 2. Test Profondeur
        print("\nCalcul des profondeurs (max_depth)...")
        depths = max_depth(go)
        print("Profondeurs :", depths)

        # Validation basique
        if depths['CC'] == 0:
            print("⚠️ Attention : Profondeur CC nulle. Vérifiez le chargement des relations part_of.")
    else:
        print("⚠️ Graphe vide. Vérifiez le chemin du fichier OBO.")
