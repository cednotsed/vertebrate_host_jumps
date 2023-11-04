rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)
require(ggtree)
require(igraph)
require(Hmisc)

threshold <- 0.15
meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")

set.seed(66)

mash_path <- "results/mash_out/viral_family_subsets/Coronaviridae.220723.n5008.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/coronaviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Circoviridae.220723.n6481.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/circoviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Rhabdoviridae.220723.n2627.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/rhabdoviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Paramyxoviridae.220723.n2357.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/paramyxoviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Poxviridae.220723.n691.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/poxviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Adenoviridae.220723.n1258.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/adenoviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Picobirnaviridae.220723.n800.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/picobirnaviridae_data.csv"

mash_path <- "results/mash_out/viral_family_subsets/Genomoviridae.220723.n1109.tsv"
save_path <- "results/clique_classification_out/cytoscape_visualisation/genomoviridae_data.csv"

mat <- fread(mash_path) %>%
  as_tibble() %>%
  column_to_rownames("#query")

acc <- names(mat)
mat <- as.matrix(mat)

# Create network
g <- graph_from_adjacency_matrix(mat, 
                                 weighted = T, 
                                 mode = "undirected")

g_filt <- delete.edges(g, which(E(g)$weight > threshold))

comm <- cluster_infomap(g_filt, modularity = F)

# Match clusters with species metadata
meta.match <- tibble(accession = comm$names, 
                   cluster = comm$membership) %>%
  left_join(meta) %>%
  mutate(is_human = ifelse(host == "Homo sapiens", "Yes", "No"))

# Invert distances
E(g_filt)$weight <- 1 - E(g_filt)$weight

# Add community and host metadata
V(g_filt)$is_human <- meta.match$is_human
V(g_filt)$viral_clique <- meta.match$cluster

# Save file
csvs <- igraph::as_data_frame(g_filt, what = "both")

merged <- csvs$vertices %>%
  dplyr::rename(from = name) %>% 
  left_join(csvs$edges)

fwrite(merged, save_path)
# fwrite(csvs$edges, "results/clique_classification_out/coronaviridae/edgelist.csv")


meta %>%
  filter(accession == "MT814882.1")
