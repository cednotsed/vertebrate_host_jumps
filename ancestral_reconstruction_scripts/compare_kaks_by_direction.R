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

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(genome_type)

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

set.seed(69)

dnds_unique <- foreach(i = seq(nrow(dnds_iter)), .combine = "bind_rows") %do% {
  row <- dnds_iter[i,]
  temp <- dnds_filt %>%
    filter(clique_name == row$clique_name,
           anc_name == row$anc_name,
           anc_state == row$anc_state,
           tip_state == row$tip_state)
  
  temp %>%
    sample_n(1, replace = F)
}

# Get unique host pairs
host_list <- dnds_unique %>%
  distinct(anc_state, tip_state)

pair_list <- list()

host_morsels <- foreach(i = seq(nrow(host_list))) %do% {
  host_pair <- as.character(host_list[i, ])
  host_pair <- host_pair[order(host_pair)]
  return(host_pair)
}

pair_filt <- unique(host_morsels)

direction_df <- foreach(i = seq(length(pair_filt)), .combine = "bind_rows") %do% {
  row <- pair_filt[[i]]
  host1 <- row[1]
  host2 <- row[2]
  
  temp <- dnds_unique %>%
    filter(anc_state == host1 & tip_state == host2|
             anc_state == host2 & tip_state == host1) %>%
    mutate(direction = ifelse(anc_state == host1, "Foward", "Reverse"))
  
  if(n_distinct(temp$direction) == 2) {
    return(temp)
  }
}

lr1 <- lm(log(kaks) ~ clique_name + direction,
          data = direction_df)

sink("results/ancestral_reconstruction_out/compare_forward_reverse.kaks.ANOVA.txt")
anova(lr1)
sink()

## Zoonotic jumps only ##
zoo_filt <- direction_df %>%
  left_join(meta %>%
              distinct(anc_state = host, 
                       anc_genus = host_genus)) %>%
  left_join(meta %>%
              distinct(tip_state = host, 
                       tip_genus = host_genus)) %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state))

lr2 <- lm(log(kaks) ~ clique_name + event_type,
          data = zoo_filt)

sink("results/ancestral_reconstruction_out/compare_anthro_zoo.kaks.ANOVA.txt")
anova(lr2)
sink()

