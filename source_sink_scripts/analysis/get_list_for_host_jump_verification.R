rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)

genome_type <- fread("data/metadata/genome_type_metadata.csv")

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv")) %>%
  left_join(genome_type)

jump_df <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.V2.csv")

# Counts
jump_df %>%
  filter(!is_jump) %>%
  distinct(clique_name) %>%
  nrow()

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

clique_counts <- jump_df %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- jump_df %>%
  # group_by(clique_name, is_jump) %>%
  # summarise(min_dist = min(patristic_dist)) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

merged_df %>% 
  mutate(anc_state = ifelse(anc_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", anc_state),
         tip_state = ifelse(tip_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", tip_state)) %>%
  filter(anc_state == "Homo sapiens" | tip_state == "Homo sapiens") %>%
  mutate(event_type = ifelse(anc_state != "Homo sapiens", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  ggplot(aes(y = log10(patristic_dist + 0.00000001), x = event_type,
             fill = event_type)) +
  geom_boxplot(position = position_nudge(x = 0.05, y = 0),
               width = 0.1,
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.1, y = 0), alpha = 1) +
  theme_classic() +
  coord_flip() +
  geom_pwc() +
  labs(x = "Is host jump?", y = "log10(min. distance)") +
  theme(legend.position = "none",
        text = element_text(color = "black")) +
  scale_fill_manual(values = c("steelblue3", "indianred3")) +
  facet_wrap(.~family)
