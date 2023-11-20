rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

meta <-fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
meta %>%
  group_by(cluster) %>%
  summarise(n = n_distinct(accession)) %>%
  filter(cluster %in% jump_df$clique_name) %>%
  arrange(n)
jump_df <- fread("results/dnds_out/all_jumps.dnds.diff_hosts.genus_counts.csv")
jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv")

  jump_df %>%
    filter(is_jump) %>%
  # distinct(clique_name, anc_name, tip_state, .keep_all = T) %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# MERS (alpha)
# https://www.nature.com/articles/nature12328
# JX869059.2
meta %>%
  filter(species == "Middle East respiratory syndrome-related coronavirus") %>%
  distinct(cluster) %>%
  filter(cluster %in% jump_df$clique_name)

meta %>%
  filter(cluster == "Coronaviridae_26") %>%
  distinct(host_genus)

# IBV (gamma)
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1007009
# DQ834384.1
meta %>%
  filter(grepl("Infectious bronchitis", genbank_title, ignore.case = T)) %>%
  distinct(cluster) %>%
  filter(cluster %in% jump_df$clique_name)

meta %>%
  filter(cluster == "Coronaviridae_6") %>%
  distinct(host_genus)

# SARS-CoV-2 (beta)
# Wuhan-hu-1
meta %>%
  filter(grepl("SARS-CoV-2", genbank_title, ignore.case = T)) %>%
  distinct(cluster) %>%
  filter(cluster %in% jump_df$clique_name)

meta %>%
  filter(cluster == "Coronaviridae_12") %>%
  distinct(host_genus)

# # PCoV (Delta)
# # https://www.nature.com/articles/s41467-022-29062-5#Sec8
# # OK546242.1
# meta %>%
#   filter(genus == "Deltacoronavirus") %>%
#   distinct(cluster, species) %>%
#   filter(cluster %in% jump_df$clique_name)
# 
# meta %>%
#   filter(cluster == "Coronaviridae_18") %>%
#   distinct(host_genus)
