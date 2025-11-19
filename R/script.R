library(igraph)
library(tidyverse)

g <- read_graph("https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/Cleandb_Luca_1_S_1_1_65_Iso_Tr_1-CC1.tgr", directed=FALSE)
g

#represnetation sous forme d'un graph

plot(g)

plot(g, layout=layout.fruchterman.reingold)



lfr <- layout_with_fr(g)
plot(g, layout=lfr, vertex.size=3, vertex.label=NA)

head(lfr)

V(g)

E(g) %>%
  head(50)

vertex_names <- read_tsv("https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/Cleandb_Luca_1_S_1_1_65_Iso_Tr_1-CC1.cod", 
                         col_names = FALSE, 
                         show_col_types = FALSE)
vertex_names

V(g)$name <- vertex_names$X2
plot(g, layout=lfr, vertex.size=3, vertex.label=V(g)$name)

igraph_options(vertex.label.cex=.6) # font size
igraph_options(vertex.label.family='sans')
igraph_options(vertex.size=3) 
igraph_options(vertex.color=NA) 

plot(g, layout=lfr)

vertex_attr(g, name="name")
V(g)$name

vertex_attr(g, name="name", index=10)
V(g)[10]$name


tb <- read_tsv("https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/511145.protein.links.experimental.CC1.tsv", 
               show_col_types = FALSE)
tb

s_exp <- graph_from_data_frame(tb, directed = FALSE)
s_exp

plot(s_exp, vertex.label = NA, vertex.size = 1)
diameter(s_exp)

distances(s_exp)[1:6, 1:6]
max(distances(s_exp))

p_info <- read_tsv("Python/data/511145.protein.info.v12.0.txt.gz", 
                   show_col_types = FALSE)
p_info %>%
  filter(preferred_name =='dnaK' | preferred_name == 'bamB')

distances(s_exp)['511145.b0014', '511145.b2512']
shortest_paths(s_exp, from='511145.b0014', to='511145.b2512')

links = read_delim("Python/data/511145.protein.links.detailed.v12.0.txt", delim = " ", show_col_types =  FALSE)
links
links.filtered <- links %>%
  filter((coexpression>800 | expression>800) & protein1<protein2)
links.filtered

sg 

df1 <- data.frame(from=c('A', 'B', 'C', 'D', 'E', 'A'), to=c('B', 'C', 'D', 'E', 'A', 'C'))
df1

g1 <- graph_from_data_frame(df1)
plot(g1, vertex.label.cex=1, vertex.size=20)

df2 <- data.frame(from=c( 'V', 'W', 'X', 'U', 'U','Y'), to=c( 'W', 'X', 'Y','V', 'W', 'U'))
g2 <- graph_from_data_frame(df2)
par(mfrow=c(1,2))
plot(g1, vertex.label.cex=1, vertex.size=20)
plot(g2, vertex.label.cex=1, vertex.size=20)

isomorphic(g1,g2)

as_adjacency_matrix(g1)
as_adjacency_matrix(g2)

canonical_permutation(g1)

as_adjacency_matrix(permute(g1, canonical_permutation(g1)$labeling))
as_adjacency_matrix(permute(g2, canonical_permutation(g2)$labeling))


mean_distance(s_exp)
components(g)

articulation_points(g)

g1 <- delete_vertices(g, articulation_points(g))
plot(g1)

plot(g - articulation_points(g))


lg <- make_line_graph((g))
par(mfrow=c(1,2))
plot(g, vertex.label=NA)
plot(lg, vertex.label=NA)


tb <- read_tsv('https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/dressing.tsv', 
               show_col_types = FALSE)
dressing <- graph_from_data_frame(tb, directed=TRUE)
plot(dressing, vertex.size=15)

p <- bfs(dressing, 
         root='sous-vetements', 
         dist=TRUE, 
         unreachable=FALSE)
p$dist

p <- dfs(dressing, 
         root='sous-vetements', 
         unreachable=TRUE, 
         father=TRUE, 
         order.out=TRUE)
p

is_dag(dressing)
topo_sort(dressing)

df.bf <- read_tsv('https://src.koda.cnrs.fr/roland.barriot/mbioinfo.graph.project/-/raw/main/data/graphe.Bellman-Ford.tsv')
bf <- graph_from_data_frame(df.bf)
bf