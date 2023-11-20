rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(Hmisc)
require(ggpubr)
require(randomcoloR)
require(see)
require(ggrepel)

genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  dplyr::rename(clique_name = cluster)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  distinct(clique_name, anc_name, tip_state, .keep_all = T)

zoo_df <- jump_df %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state))

View(zoo_jumps)

zoo_jumps <- t(zoo_df %>%
  select(anc_genus, tip_genus))

set.seed(66)

perm_df <- foreach(i = seq(500), .combine = "bind_rows") %do% {
  perm_mat <- apply(zoo_jumps, 2, function(x) sample(x, replace = F))
  as.data.frame(t(perm_mat)) %>%
    rename(anc_genus = V1, tip_genus = V2) %>%
    mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
    summarise(n_anthro = sum(event_type == "Anthroponotic"),
              n_zoo = sum(event_type == "Zoonotic")) %>%
    mutate(ratio = n_anthro / n_zoo)
}

obs_df <- zoo_df %>%
  summarise(n_anthro = sum(event_type == "Anthroponotic"),
            n_zoo = sum(event_type == "Zoonotic")) %>%
  mutate(ratio = n_anthro / n_zoo)

p_val <- sum(perm_df$ratio > obs_df$ratio) / 1000


perm_df %>%
  ggplot(aes(x = ratio)) +
  geom_histogram() +
  geom_vline(xintercept = obs_df$ratio,
             color = "red",
             lty = "dashed") +
  geom_text(x = 0.9, y = 40, label = str_glue("p={p_val}"),
            color = "red")

ggsave("results/source_sink_analysis/structured_permutation_test.png")
