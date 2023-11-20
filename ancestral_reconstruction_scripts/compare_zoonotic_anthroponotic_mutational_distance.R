rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(see)
require(ggpubr)

genome_type <- fread("data/metadata/genome_type_metadata.csv")
good_alns <- fread("results/qc_out/good_alignments.csv")

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(genome_type)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  filter(clique_name %in% good_alns$clique_name)

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

merged_df <- jump_df %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(host_counts)

# Get randomly chosen nodes
iter_df <- merged_df %>% distinct(clique_name, anc_name, anc_state, tip_state)

set.seed(69)

merged_filt <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  row <- iter_df[i,]
  temp <- merged_df %>%
    filter(clique_name == row$clique_name,
           anc_name == row$anc_name,
           anc_state == row$anc_state,
           tip_state == row$tip_state)
  
  temp %>%
    sample_n(1, replace = F)
}

zoo_filt <- merged_filt %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  left_join(genome_counts)

setNames((as.numeric(factor(zoo_filt$event_type)) - 1), zoo_filt$event_type)

logreg <- glm((as.numeric(factor(event_type)) - 1) ~ log(n_genomes) + family + log(patristic_dist),
    data = zoo_filt,
    family = "binomial")

sink("results/ancestral_reconstruction_out/compare_anthro_zoo.mutation.logreg.txt")
summary(logreg)
exp(coef(logreg)[["log(patristic_dist)"]])
sink()

# dN/dS
dnds_list <- list.files("results/dnds_out/all_jumps.temp_results/", full.names = T)
dnds_df <- foreach(file_name = dnds_list, .combine = c("bind_rows")) %do% {
  fread(file_name)
}

# Remove zeroes
dnds_filt <-  dnds_df %>%
  filter(is_jump) %>%
  mutate(ks = ifelse(ks == 0 | ks < 0, 0, ks),
         ka = ifelse(ka == 0, 0, ka)) %>%
  filter(!(ka == 0 | ks == 0)) %>%
  filter(ka != 0 & ks != 0) %>%
  mutate(kaks = ka / ks)

dnds_iter <- dnds_filt %>%
  distinct(anc_name, anc_state, tip_state, clique_name)

merged_dnds_filt <- foreach(i = seq(nrow(dnds_iter)), .combine = "bind_rows") %do% {
  row <- dnds_iter[i,]
  temp <- dnds_filt %>%
    filter(clique_name == row$clique_name,
           anc_name == row$anc_name,
           anc_state == row$anc_state,
           tip_state == row$tip_state)
  
  temp %>%
    sample_n(1, replace = F)
}

zoo_dnds_filt <- merged_dnds_filt %>%
  left_join(meta %>%
              distinct(anc_state = host, 
                       anc_genus = host_genus)) %>%
  left_join(meta %>%
              distinct(tip_state = host, 
                       tip_genus = host_genus)) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  left_join(genome_counts)

logreg <- glm((as.numeric(factor(event_type)) - 1) ~ family + log(kaks),
              data = zoo_dnds_filt,
              family = "binomial")

sink("results/ancestral_reconstruction_out/compare_anthro_zoo.dnds.logreg.txt")
summary(logreg)
exp(coef(logreg)[["log(kaks)"]])
sink()

