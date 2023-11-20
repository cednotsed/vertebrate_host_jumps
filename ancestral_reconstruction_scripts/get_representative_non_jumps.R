rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

host_meta <- fread("data/metadata/parsed_host_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv") %>%
  left_join(host_meta)
good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  as_tibble() %>%
  arrange(clique_name, desc(n_traverses))

res_dir <- "results/source_sink_analysis/ancestral_reconstructions.same_host/temp_results"

file_list <- list.files(res_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  temp <- fread(file_name)
  if(nrow(temp) > 0) {
    return(temp)
  } 
}

res_df <- bind_rows(morsels)

res_df %>% 
  filter(clique_name %in% good_alns$clique_name) %>%
  filter(clique_name %in% jump_df$clique_name) %>%
  fwrite("results/mutational_load_out/non_jumps.csv")
