rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)

genome_type <- fread("data/metadata/genome_type_metadata.csv")

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(genome_type)

# Cliques by viral family
plot_df <- meta %>%
  group_by(family, genome_type) %>%
  summarise(n_cliques = n_distinct(cluster)) %>%
  arrange(desc(n_cliques))

plot_df %>% distinct(family)

plot_df %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(y = family, x = n_cliques, fill = genome_type)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(genome_type), scales = "free") +
  theme_bw()

plot_df2 <- meta %>%
  group_by(cluster, family, genome_type) %>%
  summarise(human = sum(host_genus == "Homo") > 0,
            animal = sum(host_genus != "Homo" & host_genus != "") > 0,
            missing = sum(host_genus == "") > 0) %>%
  mutate(group = case_when(human & animal ~ "Both",
                           human & !animal ~ "Human only",
                           !human & animal ~ "Animal only",
                           !human & !animal & missing ~ "Missing only")) %>%
  group_by(family, genome_type, group) %>%
  summarise(n_cliques = n_distinct(cluster)) %>%
  mutate(family = factor(family, unique(plot_df$family)),
         group = factor(group, c("Animal only", "Human only", "Both", "Missing only")))

# Overall proportion by group
plot_df2 %>%
  group_by(group) %>%
  summarise(n = sum(n_cliques)) %>%
  mutate(prop = n / sum(n))

# Barplot by host group
plot_df2 %>%
  mutate(genome_type = case_when(grepl("ssRNA", genome_type) ~ "ssRNA",
                           grepl("dsRNA", genome_type) ~ "dsRNA",
                           grepl("ssDNA", genome_type) ~ "ssDNA",
                           grepl("dsDNA", genome_type) ~ "dsDNA")) %>%
  ggplot(aes(x = family, y = n_cliques, fill = group)) +
  geom_bar(stat = "identity", 
           position = "stack",
           color = "black") +
  facet_grid(. ~ genome_type, scales = "free", space = "free") +
  theme_bw() +
  labs(x = "Num. of viral cliques", y = "Viral family", fill = "Host") +
  scale_fill_manual(values = c("steelblue", "salmon", "goldenrod", "grey44")) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,
                                   face = "italic",
                                   angle = 90,
                                   vjust = 0.3,
                                   hjust = 1),
        axis.title = element_text(size = 14))
  
ggsave("results/clique_classification_out/clique_barplot_by_host_types.pdf", 
       dpi = 600,
       height = 5, width = 10)

ggsave("results/clique_classification_out/clique_barplot_by_host_types.png", 
       dpi = 600,
       height = 8, width = 10)


# Host count per viral clique
plot_df <- meta %>%
  filter(host_genus != "") %>%
  group_by(family, genome_type) %>%
  summarise(n_cliques = n_distinct(cluster)) %>%
  arrange(desc(n_cliques))

# Proportion of diversity involved in jumps
jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  distinct(clique_name) %>%
  mutate(contains_jump = T)

plot_df <- meta %>%
  distinct(clique_name = cluster, family) %>%
  left_join(jump_df) %>%
  mutate(contains_jump = replace_na(contains_jump, F)) %>%
  group_by(family) %>%
  summarise(prop = sum(contains_jump) / n() * 100) %>%
  arrange(desc(prop))

pal <- distinctColorPalette(n_distinct(plot_df$family))

plot_df %>%  
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family, y = prop, fill = family)) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14),
              axis.text.x = element_text(size = 14,
                                         face = "italic",
                                         angle = 90,
                                         vjust = 0.3,
                                         hjust = 1),
              axis.title = element_text(size = 14)) +
  labs(x = "Viral family", y = "% diversity involved in host jumps")

ggsave("results/meta_summary_out/perc_diversity_host_jumps.pdf",
       dpi = 600, 
       height = 4, 
       width = 8)
