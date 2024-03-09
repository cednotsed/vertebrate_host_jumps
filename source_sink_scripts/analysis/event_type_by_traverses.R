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

host_meta <- fread("data/metadata/parsed_host_metadata.csv")
genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  dplyr::rename(clique_name = cluster) 

genome_meta

prop_count <- genome_meta %>%
  group_by(clique_name) %>%
  summarise(prop_human = sum(host_genus == "Homo") / n()) %>%
  filter(!(prop_human %in% c(1, 0)))

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  as_tibble() 

zoo_df <- jump_df %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  inner_join(host_meta)

# Allow only one hostA-hostB jump per node
iter_list <- zoo_df %>%
  distinct(anc_state, tip_state, clique_name, anc_name)

set.seed(66)

zoo_filt <- foreach(i = seq(nrow(iter_list)), .combine = "bind_rows") %do% {
  row <- iter_list[i, ]
  zoo_df %>%
    filter(anc_state == row$anc_state, 
           tip_state == row$tip_state,
           clique_name == row$clique_name,
           anc_name == row$anc_name) %>%
    sample_n(1, replace = F)
}

# Observed zoonotic proportion
obs_df <- zoo_filt %>%
  group_by(event_type) %>%
  summarise(n_jumps = n()) %>%
  mutate(prop = n_jumps / sum(n_jumps))

zoo_filt %>%
  group_by(event_type, n_traverses) %>%
  summarise(n_jumps = n()) %>%
  ungroup() %>%
  complete(n_traverses, event_type, fill = list(n_jumps = 0)) %>%
  arrange(n_traverses) %>% 
  ggplot(aes(x = n_traverses, y = n_jumps, fill = event_type)) +
  geom_bar(stat = "identity", position = "dodge",
           color = "black") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Nodes traversed", y = "No. host jumps")

ggsave("results/source_sink_analysis/anthroponotic_frequency.by_traverse.pdf",
       height = 2.5, width = 4)

zoo_filt %>%
  filter(event_type == "Anthroponotic" & n_traverses == 1|event_type == "Zoonotic") %>%
  group_by(event_type) %>%
  summarise(n_jumps = n())

zoo_filt %>%
  filter(n_traverses == 1 & event_type == "Zoonotic") %>% View()

zoo_filt %>%
  left_join(prop_count) %>%
  group_by(clique_name, prop_human) %>%
  summarise(prop_anthro = sum(event_type == "Anthroponotic") / n()) %>%
  ggplot(aes(x = prop_human, y = prop_anthro)) +
  geom_point()

zoo_filt %>%
  left_join(prop_count) %>%
  mutate(prop_group = case_when(prop_human < 0.25 ~ "<0.25",
                             prop_human >= 0.25 & prop_human <= 0.5 ~ "0.25-0.5",
                             prop_human >= 0.5 & prop_human <= 0.75 ~ "0.5-0.75",
                             prop_human > 0.75 ~ ">0.75")) %>%
  mutate(prop_group = factor(prop_group, c("<0.25", "0.25-0.5", "0.5-0.75", ">0.75"))) %>%
  group_by(prop_group) %>%
  summarise(prop_anthro = sum(event_type == "Anthroponotic") / n(),
            n_anthro = sum(event_type == "Anthroponotic"),
            n_cliques = n_distinct(clique_name)) %>%
  ggplot(aes(x = prop_group, y = prop_anthro)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.5,
             lty = "dashed",
             color = "red") +
  labs(x = "Prop. of human sequences") +
  ylim(0, 1) +
  geom_text(aes(label = str_glue("n={n_cliques}")))
  
