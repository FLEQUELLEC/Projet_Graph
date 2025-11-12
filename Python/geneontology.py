#!/bin/env python
# -*- coding: utf-8 -*-
"""
Module d'analyse Gene Ontology basé sur la librairie gm.graph
=============================================================
Ce module permet de charger un fichier GO (OBO) et ses annotations (GOA)
dans un graphe orienté basé sur la classe graph de gm.py.
"""

import re
import gm

def load_OBO(filename='go-basic.obo'):
    """
    Parse un fichier OBO et construit un graphe de termes GO.

    Seules les relations 'is_a' et 'part_of' sont intégrées.
    Les termes obsolètes (is_obsolete: true) sont ignorés.
    """
    def parseTerm(lines):
        # ignore les termes obsolètes
        for l in lines:
            if l.startswith('is_obsolete: true'):
                return
        # création du nœud
        go_id = re_go_id.match(lines.pop(0)).group(1)
        go_graph.add_node(go_id, {'type': 'GOTerm'})
        go_attr = go_graph.nodes[go_id]
        # parsing des attributs
        for line in lines:
            if re_go_name.match(line):
                go_attr['name'] = re_go_name.match(line).group(1)
            elif re_go_namespace.match(line):
                go_attr['namespace'] = re_go_namespace.match(line).group(1)
            elif re_go_def.match(line):
                go_attr['def'] = re_go_def.match(line).group(1)
            elif re_go_alt_id.match(line):
                alt = re_go_alt_id.match(line).group(1)
                go_graph.alt_id[alt] = go_id
            elif re_go_is_a.match(line):
                parent_id = re_go_is_a.match(line).group(1)
                go_graph.add_edge(go_id, parent_id, {'relationship': 'is_a'})
            elif re_go_part_of.match(line):
                parent_id = re_go_part_of.match(line).group(1)
                go_graph.add_edge(go_id, parent_id, {'relationship': 'part_of'})

    # création du graphe orienté
    go_graph = gm.graph(directed=True, weighted=False)
    go_graph.alt_id = {}  # dictionnaire pour les identifiants alternatifs

    # regex
    re_go_id = re.compile(r'^id:\s+(GO:\d+)\s*$')
    re_go_name = re.compile(r'^name:\s+(.+)\s*$')
    re_go_namespace = re.compile(r'^namespace:\s+(.+)\s*$')
    re_go_def = re.compile(r'^def:\s+(.+)\s*$')
    re_go_alt_id = re.compile(r'^alt_id:\s+(GO:\d+)\s*$')
    re_go_is_a = re.compile(r'^is_a:\s+(GO:\d+)\s')
    re_go_part_of = re.compile(r'^relationship:\s+part_of\s+(GO:\d+)\s')

    # lecture du fichier
    with open(filename) as f:
        line = f.readline().rstrip()
        # saute le header
        while not line.startswith('[Term]'):
            line = f.readline().rstrip()
        buff = []
        while True:
            line = f.readline()
            if not line:
                parseTerm(buff)
                break
            line = line.rstrip()
            if line.startswith('[Term]'):
                parseTerm(buff)
                buff = []
            elif line.startswith('[Typedef]'):
                parseTerm(buff)
                break
            else:
                buff.append(line)
    return go_graph


def load_GOA(go, filename, warnings=True):
    """
    Parse un fichier GOA et ajoute les produits géniques annotés
    au graphe GO précédemment chargé.
    """
    with open(filename) as f:
        for line in f:
            if line.startswith('!'):  # ignorer les commentaires
                continue
            cols = line.rstrip().split('\t')
            gp_id = cols[1]
            gt_id = cols[4]

            # résolution d'identifiant alternatif
            if gt_id not in go.nodes:
                while gt_id not in go.nodes and gt_id in go.alt_id:
                    gt_id = go.alt_id[gt_id]

            # terme non trouvé
            if gt_id not in go.nodes:
                if warnings:
                    print(f"⚠️ Impossible de rattacher {gp_id} à {gt_id}")
                continue

            # création du produit génique
            if gp_id not in go.nodes:
                go.add_node(gp_id, {'id': gp_id, 'type': 'GeneProduct'})

            gp_attr = go.nodes[gp_id]
            gp_attr['name'] = cols[2]
            gp_attr['desc'] = cols[9]
            gp_attr['aliases'] = cols[10].split('|')

            # liaison produit → GO Term
            e_attr = go.add_edge(gp_id, gt_id, {'relationship': 'annotation'})
            if 'evidence-codes' not in e_attr:
                e_attr['evidence-codes'] = []
            e_attr['evidence-codes'].append(cols[6])


def GOTerms(go, gp_id, recursive=False):
    """
    Retourne les termes GO directement liés à un produit génique (successeurs).
    """
    if gp_id in go.nodes:
        return go.neighbors(gp_id)
    return None


def GeneProducts(go, go_id, recursive=False):
    """
    Retourne les produits géniques liés à un terme GO.
    (À compléter si besoin : prédecesseurs de type 'GeneProduct')
    """
    return None


def max_depth(go):
    """
    Calcule la profondeur maximale du graphe GO (TO DO).
    """
    return None


##### main → tests #####
if __name__ == "__main__":
    print("# Gene Ontology module tests")
    go = load_OBO("Python/data/go-basic.obo")
    print(f"{len(go.nodes)} termes chargés")
