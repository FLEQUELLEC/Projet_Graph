# Projet : Traitement de Graphes et Réseaux Biologiques

Ce dépôt contient une bibliothèque Python dédiée à la manipulation de graphes et à l'analyse de réseaux biologiques, développée dans le cadre du Master Bioinformatique.

Le projet est divisé en deux modules principaux :
1.  **`gm.py`** : Une librairie généraliste de manipulation de graphes (orientés/non-orientés, pondérés).
2.  **`geneontology.py`** : Un module spécialisé pour parser et analyser les données de la Gene Ontology (GO).

---

## 🚀 Fonctionnalités

### 1. Bibliothèque de Graphes (`gm.py`)
Cette bibliothèque légère repose sur des structures de données natives (dictionnaires) pour la flexibilité et utilise **Polars** pour une lecture rapide des fichiers de données volumineux.

* **Structure** : Graphes orientés et non orientés, gestion des poids et attributs.
* **I/O** : Lecture optimisée de fichiers délimités (CSV, TSV) via Polars.
* **Parcours** :
    * **BFS (Parcours en Largeur)** : Calcul des plus courts chemins (non pondérés).
    * **DFS (Parcours en Profondeur)** : Classification des arêtes (arbres, retour, avant, transversales) et datation des sommets.
* **Topologie & Analyse** :
    * Composantes connexes.
    * Extraction de sous-graphes induits.
    * Détection de cycles (`is_acyclic`).
    * Tri topologique (`topological_sort`) pour les DAGs.

### 2. Module Gene Ontology (`geneontology.py`)
Permet d'intégrer des connaissances biologiques aux analyses de réseaux.

* **Parsing OBO** : Chargement de l'ontologie (relations `is_a`, `part_of`).
* **Parsing GOA/GAF** : Association des produits de gènes (Gene Products) aux termes GO.
* **Gestion des identifiants** : Résolution des `alt_id` et synonymes.
* **Analyse** :
    * Navigation dans le graphe GO (Parents/Enfants).
    * Calcul de profondeur et statistiques (en cours).

---

## 📦 Installation et Prérequis

Le projet nécessite **Python 3.10+** et les bibliothèques suivantes :

```bash
pip install pandas polars
```

Note : polars est utilisé pour garantir des performances élevées lors du chargement de fichiers d'interactions protéiques massifs (ex: STRINGdb).
📂 Structure du Projet
```Plaintext

.
├── Python/
│   ├── gm.py               # Bibliothèque principale (Graph Manipulation)
│   ├── geneontology.py     # Module d'analyse GO
│   ├── test.py             # Script de tests unitaires pour gm.py
│   ├── test_geneonto.py    # Script de tests pour le module GO
│   └── data/               # Fichiers sources (OBO, GAF, Interactions)
├── R/
│   └── script.R            # Scripts comparatifs (igraph)
└── README.md
```

💻 Utilisation
Manipulation de Graphes

```Python

import gm

# Création manuelle
g = gm.graph(directed=True)
g.add_edge('A', 'B', {'weight': 2})
g.add_edge('B', 'C')

# Lecture depuis un fichier
g_file = gm.graph.read_delim('Python/data/mon_reseau.tsv')

# Parcours BFS (Plus court chemin)
resultat = g.BFS('A', cible='C')
print(resultat['chemin'])

# Détection de cycle
if not g.is_acyclic():
    print("Attention, le graphe contient des cycles !")
```

Analyse Gene Ontology
```Python

import geneontology as gom

# 1. Charger l'ontologie (Structure)
go = gom.load_OBO('Python/data/go-basic.obo')

# 2. Charger les annotations (Associations Gène -> Fonction)
gom.load_GOA(go, 'Python/data/genome_annotations.gaf')

# 3. Explorer les voisins
termes_voisins = gom.GOTerms(go, 'GO:0006915')
```

🧪 Tests

Pour valider le fonctionnement des algorithmes, exécutez les scripts de test fournis :
Bash

# Tester les algos de graphes (BFS, DFS, Cycles...)
```bash
python Python/test.py
```

# Tester le chargement GO et les annotations
```bash
python Python/test_geneonto.py
```

📝 Auteur

Florent LE QUELLEC Projet réalisé dans le cadre du cours "Traitement de Graphes et Réseaux Biologiques".
