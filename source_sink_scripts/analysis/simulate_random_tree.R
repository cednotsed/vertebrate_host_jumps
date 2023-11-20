rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggtree)
require(ape)

set.seed(67)
tree <- rtree(10)

ggtree(tree,
       ladderize = T) +
  geom_tippoint(size = 5,
                pch = 21,
                color = "black",
                fill = "indianred") +
  geom_nodepoint(size = 5,
                 pch = 22,
                 color = "black",
                 fill = "blue")

ggsave("results/source_sink_analysis/simulated_tree.pdf", dpi = 600, width = 5, height = 3)
