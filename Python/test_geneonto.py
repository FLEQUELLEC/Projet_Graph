#!/bin/env python
# -*- coding: utf-8 -*-
"""
Script de test et de démonstration pour le projet Gene Ontology.
Organisme cible : Rattus norvegicus (Rat)
"""

import gm
import geneontology as go_lib
import os
import sys

def main():
    print("=================================================================")
    print("PROJET GRAPHE & BIOINFORMATIQUE - TEST COMPLET (R. norvegicus)")
    print("=================================================================")

    # --- 1. CONFIGURATION DES FICHIERS ---
    obo_file = "Python/data/projet/go-basic.obo"
    # Nom du fichier pour Florent (R. norvegicus)
    gaf_file = "Python/data/projet/122.R_norvegicus.goa"

    if not os.path.exists(obo_file):
        print(f"ERREUR CRITIQUE: Le fichier {obo_file} est introuvable.")
        print("Téléchargez-le ici: https://purl.obolibrary.org/obo/go/go-basic.obo")
        return

    if not os.path.exists(gaf_file):
        print(f"ERREUR : Le fichier {gaf_file} est introuvable.")
        print("Veuillez télécharger le fichier du Rat (R. norvegicus) dans le dossier courant.")
        return

    # --- 2. CHARGEMENT DE L'ONTOLOGIE ---
    print(f"\n[1/5] Chargement de l'ontologie {obo_file}...")
    go = go_lib.load_OBO(obo_file)
    print(f"   -> {go.nb_nodes()} termes GO chargés.")
    print(f"   -> {go.nb_edges()} relations (is_a, part_of).")

    # --- 3. CHARGEMENT DES ANNOTATIONS (GAF) ---
    print(f"\n[2/5] Chargement des annotations {gaf_file}...")
    go_lib.load_GOA(go, gaf_file, warnings=False)

    # Compter les gènes
    genes = [n for n in go.nodes if go.nodes[n].get('type') == 'GeneProduct']
    print(f"   -> {len(genes)} gènes (GeneProducts) ajoutés au graphe.")

    if len(genes) == 0:
        print("Erreur: Aucun gène chargé. Vérifiez le format du fichier GAF.")
        return

    # --- 4. CALCUL DE LA PROFONDEUR ---
    print(f"\n[3/5] Calcul de la profondeur maximale des ontologies...")
    depths = go_lib.max_depth(go)
    print("   Résultats :")
    for namespace, d in depths.items():
        print(f"   - {namespace:20} : {d}")

    # --- 5. TEST DES FONCTIONS DE RECHERCHE ---
    print(f"\n[4/5] Tests fonctionnels sur un gène aléatoire...")

    # On prend un gène au hasard (le 100ème pour éviter les premiers parfois bizarres)
    gene_test_id = genes[100] if len(genes) > 100 else genes[0]
    gene_name = go.nodes[gene_test_id].get('name', 'Inconnu')
    print(f"   Sujet de test : {gene_test_id} ({gene_name})")

    # A. GOTerms DIRECTS
    terms_direct = go_lib.GOTerms(go, gene_test_id, recursive=False)
    print(f"   A. Termes GO directs     : {len(terms_direct)}")

    # B. GOTerms RÉCURSIFS (Ancêtres)
    terms_all = go_lib.GOTerms(go, gene_test_id, recursive=True)
    print(f"   B. Termes GO (avec ancêtres) : {len(terms_all)}")

    if len(terms_direct) > 0:
        # On choisit un terme GO parmi ceux du gène pour faire le test inverse
        term_test_id = terms_direct[0]
        term_name = go.nodes[term_test_id].get('name', 'Inconnu')
        print(f"\n[5/5] Tests fonctionnels sur un terme GO : {term_test_id}")
        print(f"      Nom : {term_name}")

        # C. GeneProducts DIRECTS
        prods_direct = go_lib.GeneProducts(go, term_test_id, recursive=False)
        print(f"   C. Gènes annotés directement : {len(prods_direct)}")

        # D. GeneProducts RÉCURSIFS (Descendants)
        prods_all = go_lib.GeneProducts(go, term_test_id, recursive=True)
        print(f"   D. Gènes (avec descendants)  : {len(prods_all)}")

    print("\n=================================================================")
    print("TESTS TERMINÉS AVEC SUCCÈS")
    print("=================================================================")

if __name__ == "__main__":
    main()
