rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggtreeExtra)

set.seed(66)
tree <- rtree(10)

# Match metadata to tips
n_human <- 5
tip_states <- c(rep("Human", n_human), rep("Animal", 10 - n_human))
tips <- tree$tip.label

# Actual ancestral reconstruction
x <- setNames(tip_states, tips)

fitER <- ape::ace(x, tree,
                  model="ER",
                  type="discrete",
                  method = "ML")

ancstats <- as.data.frame(fitER$lik.anc)
ancstats$node <- (1:tree$Nnode) + Ntip(tree)

# Plot ancestral reconstructions
dd <- data.frame(Accession = tree$tip.label, 
                 host = tip_states)

cols <- setNames(c("red", "blue"), sort(unique(x)))
pies <- nodepie(ancstats, 
                cols = 1:(ncol(ancstats) - 1),
                color = cols)

actual <- ggtree(tree, 
       branch.length = "none",
       ladderize = T,
       size = 0.001,
       color = "darkslategrey") %<+% dd + 
  geom_tippoint(aes(color = host)) + 
  scale_color_manual(values = cols) +
  theme(legend.position = "right") +
  geom_inset(pies, 
             width=0.08, 
             height=0.08) +
  labs(title = "Observed tree")

ggsave("results/source_sink_analysis/effects_of_permutation/actual_tree.pdf",
       width = 5, height = 5,
       plot = actual)

# Permuted ancestral reconstruction
tip_states <- sample(tip_states, length(tip_states), F)
x <- setNames(tip_states, tips)

fitER <- ape::ace(x, tree,
                  model="ER",
                  type="discrete",
                  method = "ML")

ancstats <- as.data.frame(fitER$lik.anc)
ancstats$node <- (1:tree$Nnode) + Ntip(tree)

# Plot ancestral reconstructions
dd <- data.frame(Accession = tree$tip.label, 
                 host = tip_states)

cols <- setNames(c("red", "blue"), sort(unique(x)))
pies <- nodepie(ancstats, 
                cols = 1:(ncol(ancstats) - 1),
                color = cols)

perm <- ggtree(tree, 
       branch.length = "none",
       ladderize = T,
       size = 0.001,
       color = "darkslategrey") %<+% dd + 
  geom_tippoint(aes(color = host)) + 
  scale_color_manual(values = cols) +
  theme(legend.position = "right") +
  geom_inset(pies, 
             width=0.08, 
             height=0.08) +
  labs(title = "Permuted tree")
ggsave("results/source_sink_analysis/effects_of_permutation/permuted_tree.pdf",
       width = 5, height = 5,
       plot = perm)

ggpubr::ggarrange(actual, perm)
ggsave("results/source_sink_analysis/effects_of_permutation/effects_of_permutation.pdf", 
       width = 8, height = 5)



