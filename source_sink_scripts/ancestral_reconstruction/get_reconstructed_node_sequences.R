rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(treeio)
require(phytools)
require(castor)
require(Biostrings)
require(foreach)

out_dir <- "results/source_sink_analysis/ancestral_reconstructions/ancestral_sequences"
state_dir <- "data/trees/source_sink_mini_trees/without_outgroup.masked"
state_paths <- list.files(state_dir, ".state")

foreach(state_path = state_paths) %do% {
  # state_path = state_paths[1]
  clique_name <- str_split(state_path, "\\.")[[1]][1]
  print(clique_name)
  state_df <- fread(str_glue("{state_dir}/{state_path}"))
  
  node_list <- deframe(state_df %>%
    distinct(Node))
    
  node_seqs <- foreach(node_name = node_list,
                       .combine = "c") %do% {
    # node_name = node_list[1]
    state_filt <- deframe(state_df %>%
      filter(Node == node_name) %>%
      select(State))
    
    seq_string <- paste0(state_filt, collapse = "")
    temp_seq <- DNAStringSet(seq_string)
    names(temp_seq) <- node_name
    return(temp_seq)
  }
  
  writeXStringSet(node_seqs, str_glue("{out_dir}/{clique_name}.anc_seqs.fna"))
}


