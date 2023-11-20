rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(data.table)
require(tidyverse)
require(ape)
require(randomcoloR)
require(castor)
require(foreach)
require(adephylo)

genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  rename(clique_name = cluster)

# No. of genomes
genome_meta %>%
  nrow()

# No. of families
genome_meta %>%
  distinct(family) %>%
  nrow()

# No. of host orders
genome_meta %>%
  filter(is_vertebrate) %>%
  distinct(host_order) %>%
  filter(host_order != "") %>%
  nrow()

# Concordance with ICTV
ictv <- fread("data/metadata/ICTV_Master_Species_List_2022_MSL38.v2.060823.csv")
genome_meta %>%
  distinct(species) %>%
  summarise(prop = sum(species %in% ictv$Species) / n())

# No. of cliques
genome_meta %>%
  distinct(clique_name)

# No. of viral families
genome_meta %>%
  distinct(family) %>%
  filter(family != "")

# No. of host jumps
good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  filter(is_jump) %>%
  distinct(anc_name, tip_state, anc_state, clique_name)

jump_df %>%
  nrow()

jump_df %>%
  distinct(clique_name) %>%
  nrow()
