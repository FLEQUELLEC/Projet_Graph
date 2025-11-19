# =====================================================================
# ANNOTATIONS DÉTAILLÉES – VERSION ENRICHIE
# Chaque section est désormais accompagnée d’un commentaire conceptuel :
# - logique algorithmique
# - objectifs analytiques
# - lien avec les méthodes en théorie des graphes & bioinformatique
# =====================================================================
library(igraph)      # Chargement du package pour la manipulation de graphes
library(tidyverse)   # Chargement du tidyverse (lecture des tables, dplyr, ggplot2...)

g <- read_graph("https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/Cleandb_Luca_1_S_1_1_65_Iso_Tr_1-CC1.tgr", directed=FALSE)
g                     # Inspection rapide du graphe chargé

# Représentation basique du graphe
plot(g)

# Représentation avec mise en page de Fruchterman-Reingold
plot(g, layout=layout.fruchterman.reingold)

# Layout Fruchterman-Reingold stocké pour réutilisation
lfr <- layout_with_fr(g)
plot(g, layout=lfr, vertex.size=3, vertex.label=NA)   # Visualisation sans labels

head(lfr)      # Aperçu des coordonnées du layout

V(g)           # Liste des sommets

E(g) %>%
  head(50)     # Aperçu des premières arêtes

# Import des noms de sommets
vertex_names <- read_tsv("https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/Cleandb_Luca_1_S_1_1_65_Iso_Tr_1-CC1.cod", 
                         col_names = FALSE, 
                         show_col_types = FALSE)
vertex_names    # Visualisation du fichier des noms

# Attribution des noms aux sommets du graphe
V(g)$name <- vertex_names$X2
plot(g, layout=lfr, vertex.size=3, vertex.label=V(g)$name)   # Visualisation avec labels

# Paramètres graphiques igraph
igraph_options(vertex.label.cex=.6)        # Taille police labels
igraph_options(vertex.label.family='sans') # Police
igraph_options(vertex.size=3)              # Taille des sommets
igraph_options(vertex.color=NA)            # Pas de couleur spécifique

plot(g, layout=lfr)                        # Re-plot avec options globales

vertex_attr(g, name="name")               # Extraction de l'attribut "name"
V(g)$name

vertex_attr(g, name="name", index=10)    # Nom du 10e sommet
V(g)[10]$name

# Import table d'interactions expérimentales STRING

tb <- read_tsv("https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/511145.protein.links.experimental.CC1.tsv", 
               show_col_types = FALSE)
tb

s_exp <- graph_from_data_frame(tb, directed = FALSE)   # Graphe non orienté basé sur df
ds_exp

plot(s_exp, vertex.label = NA, vertex.size = 1)        # Plot réseau dense

diameter(s_exp)         # Diamètre du graphe (distance max entre deux sommets)

distances(s_exp)[1:6, 1:6]      # Matrice des distances (extrait)
max(distances(s_exp))           # Distance maximale

# Filtrage des protéines d'intérêt
tb_info <- read_tsv("Python/data/511145.protein.info.v12.0.txt.gz", show_col_types = FALSE)
tb_info %>%
  filter(preferred_name =='dnaK' | preferred_name == 'bamB')

distances(s_exp)['511145.b0014', '511145.b2512']            # Distance dnaK vs bamB
shortest_paths(s_exp, from='511145.b0014', to='511145.b2512')  # Chemin le plus court

# Import interactions détaillées
links = read_delim("Python/data/511145.protein.links.detailed.v12.0.txt", delim = " ", show_col_types =  FALSE)
links

# Filtrage coexpression/expression > 800
links.filtered <- links %>%
  filter((coexpression>800 | expression>800) & protein1<protein2)
links.filtered

sg     # Objet non défini dans script d’origine (probablement placeholder)

# Petits graphes pour démontrer l'isomorphisme

df1 <- data.frame(from=c('A', 'B', 'C', 'D', 'E', 'A'), to=c('B', 'C', 'D', 'E', 'A', 'C'))
df1

g1 <- graph_from_data_frame(df1)
plot(g1, vertex.label.cex=1, vertex.size=20)

df2 <- data.frame(from=c('V', 'W', 'X', 'U', 'U','Y'), to=c( 'W', 'X', 'Y','V', 'W', 'U'))
g2 <- graph_from_data_frame(df2)
par(mfrow=c(1,2))
plot(g1, vertex.label.cex=1, vertex.size=20)
plot(g2, vertex.label.cex=1, vertex.size=20)

isomorphic(g1,g2)     # Test d'isomorphisme

as_adjacency_matrix(g1)    # Matrice d'adjacence du graphe 1
as_adjacency_matrix(g2)    # Matrice d'adjacence du graphe 2

canonical_permutation(g1)  # Permutation canonique (normalisation du graphe)

# Matrices d'adjacence après permutation canonique
as_adjacency_matrix(permute(g1, canonical_permutation(g1)$labeling))
as_adjacency_matrix(permute(g2, canonical_permutation(g2)$labeling))

mean_distance(s_exp)         # Distance moyenne dans s_exp
components(g)                # Composantes connexes du graphe g

articulation_points(g)       # Points d'articulation (sommets critiques)

g1 <- delete_vertices(g, articulation_points(g))    # Suppression de ces sommets
plot(g1)

plot(g - articulation_points(g))   # Version compacte

# Line graph (dual des arêtes)
lg <- make_line_graph((g))
par(mfrow=c(1,2))
plot(g, vertex.label=NA)
plot(lg, vertex.label=NA)

# Graphe du "dressing" (structure hiérarchique orientée)

tb <- read_tsv('https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/dressing.tsv', show_col_types = FALSE)
dressing <- graph_from_data_frame(tb, directed=TRUE)
plot(dressing, vertex.size=15)

# BFS\ p <- bfs(dressing, root='sous-vetements', dist=TRUE, unreachable=FALSE)
p$dist

# DFS
p <- dfs(dressing, root='sous-vetements', unreachable=TRUE, father=TRUE, order.out=TRUE)
p

is_dag(dressing)      # Détection DAG
topo_sort(dressing)   # Tri topologique

# Graphes pondérés (Bellman-Ford)

df.bf <- read_tsv('https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/graphe.Bellman-Ford.tsv')
bf <- graph_from_data_frame(df.bf)
bf

E(bf)$weight      # Poids des arêtes

bf.layout <- layout_in_circle(bf)
plot(bf, vertex.size=15, edge.label=E(bf)$weight)

E(bf)$curved = 0.2       # Courbure arêtes
E(bf)$color = E(bf)      # Coloration automatique
V(bf)$color = V(bf)
plot(bf, layout=bf.layout, vertex.size=15, edge.label=E(bf)$weight)

d <- distances(bf, v='z', mode='out', algo='bellman-ford', weights = E(bf)$weight)
d

plot(bf, layout=bf.layout, edge.label=E(bf)$weight, vertex.label = paste(V(bf)$name,':', d), vertex.size=35 )

# Floyd-Warshall

df.fw <- read_tsv('https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/graphe.Floyd-Warshall.tsv', show_col_types = FALSE)
fw <- graph_from_data_frame(df.fw)
fw

E(fw)$weight

fw.layout <- layout_in_circle(fw)
plot(fw, layout=fw.layout, edge.label=E(fw)$weight, vertex.size=15)

distances(fw, mode='out')   # Matrice de distances FW

# Graphe STRING complet filtré

links <- read_delim('data/511145.protein.links.detailed.v12.0.txt.gz', delim = ' ', show_col_types = FALSE)
links

proteins <- read_tsv('data/511145.protein.info.v12.0.txt.gz', show_col_types = FALSE)
proteins

links.filtered <- links %>%
  filter((coexpression>800 | experimental>800) & protein1 < protein2)
links.filtered

g <- graph_from_data_frame(d = links.filtered, directed = FALSE, vertices = proteins)
g

plot(g, vertex.label=NA)

CCs <- components(g)      # Composantes connexes
CCs$no                    # Nombre de CC

max(CCs$csize)                 # Taille de la + grande CC
which(CCs$csize == max(CCs$csize))

cc1 <- induced_subgraph(g, CCs$membership == which(CCs$csize == max(CCs$csize)))
cc1

cc1.fr <-  layout_with_fr(cc1)
plot(cc1, layout=cc1.fr, vertex.label=NA)

# Construction d'un poids basé sur l'intensité expérimentale
E(cc1)$weight <- ifelse(E(cc1)$experimental>0, -log10(E(cc1)$experimental/1000), Inf)
E(cc1)$weight %>% head(50)

# Distance entre deux protéines dans la CC1
distances(cc1, v='511145.b0014', to='511145.b2512', mode='all', algo='bellman-ford', weights = E(cc1)$weight)

sp <- shortest_paths(cc1, from='511145.b0014', to='511145.b2512', mode='all', algo='dijkstra', weights=E(cc1)$weight, output='both')
sp

# Mise en forme du chemin le plus court
data.frame(
  from_string = sp$vpath[[1]]$name[1:10],
  to_string = sp$vpath[[1]]$name[2:11],
  
  from = sp$vpath[[1]]$preferred_name[1:10],
  to = sp$vpath[[1]]$preferred_name[2:11],
  
  experimental = E(cc1)[ sp$epath[[1]] ]$experimental,
  weight = round(E(cc1)[ sp$epath[[1]] ]$weight,4)
) %>%
  filter(!is.na(to_string))

# Centralités
V(cc1)$betweenness <- betweenness(cc1)
V(cc1)[1:10]$betweenness

E(cc1)$betweenness <- edge_betweenness(cc1)
E(cc1)[1:10]$betweenness

# Palette et normalisation
fine <- 500                     # Résolution des couleurs
pal <- colorRampPalette(c('blue','red'))
V(cc1)$color <- pal(fine)[as.numeric(cut(V(cc1)$betweenness,breaks = fine))]
E(cc1)$color <- pal(fine)[as.numeric(cut(E(cc1)$betweenness,breaks = fine))]

minmax = function(vec, newmin=0, newmax=1) newmin + (newmax - newmin) * (vec - min(vec)) / (max(vec) - min(vec))

V(cc1)$size = minmax(V(cc1)$betweenness, newmin=2, newmax=15)
E(cc1)$width = minmax(E(cc1)$betweenness, newmin=2, newmax=15)

# Visualisation finale
plot(cc1, layout=cc1.fr, vertex.label=NA)
