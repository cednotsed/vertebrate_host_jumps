rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(see)

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
  group_by(clique_name, is_jump) %>%
  summarise(min_dist = min(patristic_dist)) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

merged_df %>%
  group_by(family) %>%
  summarise(n = n_distinct(clique_name)) %>%
  arrange(desc(n))

merged_df <- jump_df %>%
  group_by(clique_name, is_jump) %>%
  summarise(min_dist = min(patristic_dist)) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

type_list <- c("ssDNA", "dsDNA", "+ssRNA", "-ssRNA")

plots <- foreach(type = type_list) %do% {
  test <- merged_df %>% filter(genome_type == type)

  formula_string <- "n_hosts ~ log(n_genomes) + family + log(min_dist)"
  jump_mod <- glm(as.formula(formula_string),
                  data = test %>% filter(is_jump),
                  family = poisson())
  
  nonjump_mod <- glm(as.formula(formula_string),
                     data = test %>% filter(!is_jump),
                     family = poisson())
  
  slope_jump <- signif(coef(jump_mod)[["log(min_dist)"]], 3)
  z_jump <- signif(summary(jump_mod)$coefficient["log(min_dist)", "z value"], 3)
  df_jump <- nrow(test %>% filter(is_jump)) - nrow(coef(summary(jump_mod)))
  p_jump <- signif(summary(jump_mod)$coefficient["log(min_dist)", "Pr(>|z|)"], 3)
  
  slope_nonjump <- signif(coef(nonjump_mod)[["log(min_dist)"]], 3)
  z_nonjump <- signif(summary(nonjump_mod)$coefficient["log(min_dist)", "z value"], 3)
  df_nonjump <- nrow(test %>% filter(!is_jump)) - nrow(coef(summary(nonjump_mod)))
  p_nonjump <- signif(summary(nonjump_mod)$coefficient["log(min_dist)", "Pr(>|z|)"], 3)
  
  test %>%
    ggplot(aes(x = n_hosts, y = log10(min_dist), fill = is_jump, color = is_jump)) +
    geom_point(alpha = 0.9, color = "black", pch = 21) +
    geom_smooth(method = "lm", fill = "grey") +
    scale_fill_manual(values = c("steelblue3", "indianred3")) +
    scale_color_manual(values = c("steelblue3", "indianred3")) +
    labs(x = "Host range (No. host genera)", y = "log10(min dist.)",
         color = "Is host jump?",
         title = str_glue("{type}\nHost jumps: slope={slope_jump}, Z={z_jump}, d.f.={df_jump}, p={p_jump}\nNon-host jumps: slope={slope_nonjump}, Z={z_nonjump}, d.f.={df_nonjump}, p={p_nonjump}")) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none", 
          axis.title = element_text(face = "bold"))
}


plt <- egg::ggarrange(plots = plots)

ggsave("results/mutational_load_out/mutational_threshold.by_groups.pdf", 
       width = 14, height = 5, dpi = 600,
       plot = plt)
