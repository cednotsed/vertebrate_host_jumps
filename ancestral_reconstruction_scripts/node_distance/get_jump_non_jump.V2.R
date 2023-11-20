rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)

genome_type <- fread("data/metadata/genome_type_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(genome_type)

good_alns <- fread("results/qc_out/good_alignments.csv")

# Get host counts
host_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  filter(host_genus != "") %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_genus))

genome_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(n_genomes = n_distinct(accession))

jump_dir <- "results/ancestral_reconstruction_out/node_distance/ancestral_reconstructions.node_distance/temp_results/"

jump_list <- list.files(jump_dir, full.names = T)

# Parse files
jump_morsels <- foreach(file_name = jump_list) %do% {
  temp <- fread(file_name)
  if(nrow(temp) > 0) {
    return(temp)
  } 
}

jump_df <- bind_rows(jump_morsels) %>%
  left_join(meta %>%
              distinct(anc_state = host,
                     anc_genus = host_genus,
                     anc_family = host_family,
                     anc_order = host_order,
                     anc_vertebrate = is_vertebrate)) %>%
  left_join(meta %>%
              distinct(tip_state = host,
                     tip_genus = host_genus,
                     tip_family = host_family,
                     tip_order = host_order,
                     tip_vertebrate = is_vertebrate)) %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  mutate(is_jump = T)

jump_df %>%
  fwrite("results/ancestral_reconstruction_out/host_jump_lists/node_distance.parsed_jumps.csv")


