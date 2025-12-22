# Analyse de Graphes et Gene Ontology

Ce dépôt contient les bibliothèques Python développées dans le cadre de l'UE "Algorithmique des Graphes" (Master Bioinformatique, Université de Toulouse). Le projet vise à modéliser des réseaux biologiques et à explorer les données de la Gene Ontology (GO) pour l'organisme *Rattus norvegicus*.

## Description des Modules

Le projet s'articule autour de deux modules principaux :

### 1. `gm.py` (GraphMaster)

Une librairie généraliste pour la manipulation de graphes orientés et non-orientés. Elle privilégie l'utilisation de structures natives (dictionnaires) pour la flexibilité et intègre **Polars** pour l'optimisation des lectures de fichiers.

* **Structure de données** : Listes d'adjacence basées sur des dictionnaires imbriqués.
* **Parcours** : BFS (Largeur) et DFS (Profondeur).
* **Algorithmes** :
* Calcul des plus courts chemins.
* Détection de cycles (`is_acyclic`).
* Tri topologique (`topological_sort`) pour les graphes orientés acycliques (DAG).
* Calcul du coefficient de clustering et extraction de composantes connexes.

### 2. `geneontology.py`

Un module métier dédié à l'analyse biologique, héritant des capacités de `gm.py`.

* **Parsing** : Lecture robuste des fichiers OBO (structure de l'ontologie) et GAF (annotations).
* **Gestion des identifiants** : Prise en charge des références anticipées et des identifiants alternatifs.
* **Navigation** : Recherche récursive des ancêtres (parents) et descendants (enfants).
* **Topologie** : Calcul optimisé de la profondeur maximale des sous-ontologies via index inversé.

## Installation

Le projet nécessite **Python 3.10+**.

Installation des dépendances (utilisées pour le parsing rapide des fichiers CSV/TSV) :

```bash
pip install -r requirements.txt

```

*(Contenu de requirements.txt : `pandas`, `polars`)*

## Structure du Projet

```text
.
├── gm.py               # Librairie de manipulation de graphes
├── geneontology.py     # Module d'analyse GO
├── test_project.py     # Script de validation et statistiques
├── requirements.txt    # Dépendances
└── data/               # Dossier contenant les fichiers .obo et .goa

```

## Exemples d'utilisation

### Manipulation de Graphes (`gm.py`)

```python
import gm

# Création d'un graphe orienté
g = gm.graph(directed=True)
g.add_edge('A', 'B', {'weight': 2})
g.add_edge('B', 'C')

# Lecture depuis un fichier TSV
# g = gm.graph.read_delim('data/interactions.tsv')

# Parcours et analyse
if g.is_acyclic():
    print("Ordre topologique :", g.topological_sort())

chemin = g.BFS('A', cible='C')
print(f"Chemin le plus court : {chemin['chemin']}")

```

### Analyse Gene Ontology (`geneontology.py`)

```python
import geneontology as go_lib

# 1. Chargement des données
# Le fichier OBO définit la structure, le fichier GOA lie les gènes aux termes
go = go_lib.load_OBO('data/go-basic.obo')
go_lib.load_GOA(go, 'data/122.R_norvegicus.goa')

# 2. Récupération des annotations (avec récursion vers les ancêtres)
genes = go_lib.GeneProducts(go, 'GO:0006915', recursive=True)
print(f"Gènes associés au terme GO:0006915 : {len(genes)}")

# 3. Calcul des profondeurs maximales (Optimisé)
depths = go_lib.max_depth(go)
print("Profondeurs par ontologie :", depths)
# Ex: {'BP': 18, 'MF': 12, 'CC': 14}

```

## Auteur

**Florent LE QUELLEC**
Master Bioinformatique et Génomique Environnementale
Année universitaire 2024-2025
