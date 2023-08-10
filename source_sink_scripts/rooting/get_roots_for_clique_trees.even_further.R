rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(treeio)
require(castor)
require(foreach)
require(phangorn)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)

meta <- fread(str_glue("results/clique_classification_out/final_cluster_metadata.220723.csv"))
tree_dir <- "data/trees/source_sink_mini_trees/with_buffer_outgroup.even_further/unrooted/"
tree_paths <- list.files(tree_dir, ".nwk")
# tree_path <- tree_paths[1]
tree_path

root_df <- fread("results/source_sink_analysis/source_sink_results.even_further.roots.curated.csv")

root_filt <- root_df %>%
  filter(root_resolved)

morsels <- foreach(clique_name = unique(root_filt$cluster)) %do% {
  tree_path <- tree_paths[grepl(str_glue("{clique_name}\\."), tree_paths)]
  
  # Get NJ tree
  tree <- read.tree(str_glue("{tree_dir}/{tree_path}"))
  tips <- tree$tip.label
  
  # Get outgroups
  outgroups <- deframe(root_filt %>%
    filter(cluster == clique_name) %>%
    select(outgroups))
  
  outgroups <- str_split(outgroups, ";", simplify = T)[1, ]
  
  outgrp_accs <- deframe(meta %>%
      filter(cluster %in% outgroups) %>%
      filter(accession %in% tips) %>%
      select(accession))
  
  not_clique <- deframe(meta %>%
                          filter(cluster != clique_name) %>%
                          filter(accession %in% tips) %>%
                          select(accession))
  
  # Root if monophyletic
  is_mono <- is.monophyletic(tree, outgrp_accs)
    
  if (is_mono) {
    tree <- root(tree, outgrp_accs, resolve.root = T)
    tree <- drop.tip(tree, not_clique)
    
    # Get root tip within clique
    root_index <- find_root(tree)
    first_descendants <- Descendants(tree, root_index, type = "children")
    tip_depths <- get_all_distances_to_root(tree, as_edge_count=T)[1:Ntip(tree)]
    shallowest_tips <- tree$tip.label[which(tip_depths == min(tip_depths))]
    root_tip <- shallowest_tips[1]
    
    return(tibble(cluster = clique_name, 
                  root_resolved = is_mono,
                  root_tip = root_tip))
  } else {
    print(str_glue("Eff! {clique_name} has a non-mono outgroup"))
  }
}

# Merge results
root_df %>%
  select(cluster, hosts, n_genomes, root_resolved, outgroups) %>%
  left_join(bind_rows(morsels)) %>%
  fwrite("results/source_sink_analysis/final_source_sink_roots.csv")

