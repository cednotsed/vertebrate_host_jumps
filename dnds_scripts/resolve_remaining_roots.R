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

root_df <- fread("results/source_sink_analysis/source_sink_results.roots.curated.csv")
meta <- fread(str_glue("results/clique_classification_out/final_cluster_metadata.220723.csv"))

previously_unresolved <- list.files("results/source_sink_analysis/clique_trees/unresolved_outgroups/")
previously_unresolved <- str_split(previously_unresolved, "\\.", simplify = T)[, 1]

to_do <- deframe(root_df %>% 
                   filter(root_resolved) %>%
                   select(cluster) %>%
                   filter(cluster %in% previously_unresolved))

to_do
tree_dir <- "data/trees/source_sink_mini_trees/with_buffer_outgroup/unrooted/"
tree_paths <- list.files(tree_dir, ".nwk", full.names = T)
tree_path <- tree_paths[1]
tree_path

morsels <- foreach(clique = to_do) %do% {
  # clique = to_do[37]
  print(clique)
  # Get NJ tree
  tree_path <- tree_paths[grepl(str_glue("{clique}\\."), tree_paths)]
  tree <- read.tree(str_glue("{tree_path}"))
  
  clique_accs <- deframe(meta %>%
    filter(cluster == clique) %>%
    select(accession))
  
  tips <- tree$tip.label
  
  # Get outgroups
  outgrp_string <- deframe(root_df %>%
    filter(cluster == clique) %>%
    select(outgroups))
  
  outgrp_cluster <- str_split(outgrp_string, ";")[[1]][1]
  
  outgrp_accs <- deframe(meta %>%
                           filter(accession %in% tips) %>%
                           filter(cluster %in% outgrp_cluster) %>%
                           select(accession))
  
  # Check if only single host
  tips_filt <- tips[tips %in% clique_accs]
  not_clique <- tips[!(tips %in% clique_accs)]
  
  n_hosts <- meta %>% 
    filter(accession %in% tips_filt) %>%
    filter(host_genus != "") %>%
    summarise(n_hosts = n_distinct(host_genus))
  
  is_multihost <- deframe(n_hosts) > 1
  
  if (is_multihost) {
    # Root if monophyletic
    is_mono <- is.monophyletic(tree, outgrp_accs)
    
    if (is_mono) {
      tree <- root(tree, outgrp_accs, resolve.root = T)
      tree <- drop.tip(tree, not_clique)
      
      # Get root tip within clique
      root_index <- find_root(tree)
      tip_depths <- get_all_distances_to_root(tree, as_edge_count=T)[1:Ntip(tree)]
      shallowest_tips <- tree$tip.label[which(tip_depths == min(tip_depths))]
      root_tip <- shallowest_tips[1]
      
      return(tibble(cluster = clique, 
                    root_resolved = is_mono,
                    outgroups = outgrp_string,
                    root_tip = root_tip))
    } else {
      print(str_glue("root problem with {clique}"))
      return(NULL)
    }
  }
}

# Merge results
root_df %>%
  left_join(bind_rows(morsels)) %>%
  fwrite("results/source_sink_analysis/source_sink_results.roots.curated.resolved.csv")

