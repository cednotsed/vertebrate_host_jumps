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
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic"))

plot_df <- zoo_df %>% 
  separate(clique_name, c("family"), "_", remove = F) %>%
  group_by(family) %>%
  summarise(n_anthro = sum(event_type == "Anthroponotic"),
            n_zoo = sum(event_type == "Zoonotic"),
            prop_anthro = sum(event_type == "Anthroponotic") / n()) %>%
  pivot_longer(!c(prop_anthro, family), names_to = "event_type", values_to = "count") %>%
  arrange(desc(prop_anthro))


plot_df %>%
  mutate(event_type = ifelse(event_type == "n_anthro", "Anthroponotic", "Zoonotic")) %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family, y = count, fill = event_type)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           color = "black") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Viral family", y = "No. host jumps",
       fill = "Direction of jump")
  
ggsave("results/source_sink_analysis/anthroponotic_frequency.by_family.pdf",
       height = 3, width = 4)
ggsave("results/source_sink_analysis/anthroponotic_frequency.by_family.png",
       height = 5, width = 4)

# Without SC2 and flu
plot_df <- zoo_df %>% 
  filter(!(clique_name %in% c("Coronaviridae_12", "Orthomyxoviridae_1", "Orthomyxoviridae_2"))) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  group_by(family) %>%
  summarise(n_anthro = sum(event_type == "Anthroponotic"),
            n_zoo = sum(event_type == "Zoonotic"),
            prop_anthro = sum(event_type == "Anthroponotic") / n()) %>%
  pivot_longer(!c(prop_anthro, family), names_to = "event_type", values_to = "count") %>%
  arrange(desc(prop_anthro))


plot_df %>%
  mutate(event_type = ifelse(event_type == "n_anthro", "Anthroponotic", "Zoonotic")) %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family, y = count, fill = event_type)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Viral family", y = "No. host jumps",
       fill = "Direction of jump")

ggsave("results/source_sink_analysis/anthroponotic_frequency.by_family.no_SC2_flu.pdf",
       height = 3, width = 4)
ggsave("results/source_sink_analysis/anthroponotic_frequency.by_family.no_SC2_flu.png",
       height = 5, width = 4)
# plot_df %>% 
#   distinct(prop_anthro, family) %>%
#   summarise(prop_families = sum(prop_anthro >= 0.5),
#             total_families = n())
