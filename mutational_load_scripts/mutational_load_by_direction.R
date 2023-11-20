rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  rename(clique_name = cluster)

genome_counts <- meta %>%
  group_by(clique_name) %>%
  summarise(n = n_distinct(accession),
            prop_animal = sum(host_genus != "Homo") / n_distinct(accession)) %>%
  arrange(prop_animal)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  distinct(clique_name, anc_name, tip_state, .keep_all = T) %>%
  filter(anc_genus != "" & tip_genus != "")

host_list <- unique(c(jump_df$anc_genus, jump_df$tip_genus))
host_list
host_combn <- t(combn(host_list, 2, simplify = T)) %>%
  as_tibble()

foreach(i = seq(nrow(host_combn))) %do% {
  row <- host_combn[i, ]
  host1 <- row$V1
  host2 <- row$V2
  
  jump_df %>%
    filter(anc_genus == host1 & tip_genus == host2|
             anc_genus == host2 & tip_genus == host1) %>%
    mutate(direction = ifelse(anc_genus == host1, "Foward", "Reverse"))
  
}