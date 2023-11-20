rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)
require(ggtree)
require(igraph)
require(Hmisc)

set.seed(66)
g <- sample_gnp(10, 5/10) %du% sample_gnp(10, 5/9)
g <- add_edges(g, c(1, 12))
g <- induced_subgraph(g, subcomponent(g, 1))

V(g)$color <- "steelblue3"
V(g)$color[1:10] <- "goldenrod"
E(g)$weight <- 2

g_plus <- g + 
  edge(4, 18, weight = 1) +
  edge(1, 16, weight = 1) +
  edge(1, 18, weight = 1) +
  edge(3, 16, weight = 1) +
  edge(3, 18, weight = 1) +
  edge(13, 3, weight = 1) +
  edge(18, 8, weight = 1) +
  edge(12, 4, weight = 1) +
  edge(6, 12, weight = 1) +
  edge(2, 16, weight = 1) +
  edge(2, 18, weight = 1) +
  edge(6, 16, weight = 1) +
  edge(6, 11, weight = 1) +
  edge(4, 13, weight = 1)

pdf("results/clique_classification_out/simulated_graphs.pdf",
    width = 5, height = 5)
plot(g, 
     vertex.label = NA,
     edge.arrow.size=0.02,
     vertex.size = 5,
     edge.width = E(g)$weight,
     vertex.color = V(g)$color,
     layout = layout_with_fr(g))

plot(g_plus, 
     vertex.label = NA,
     edge.arrow.size=0.02,
     vertex.size = 5,
     edge.width = E(g_plus)$weight,
     layout = layout_with_fr(g_plus))
dev.off()
