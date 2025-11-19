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