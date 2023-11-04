rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(see)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv"))

genome_type <- fread("data/metadata/genome_type_metadata.csv")

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

genome_length <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(median_length = median(genome_length))

dnds_list <- list.files("results/dnds_out/all_jumps.by_gene.temp_results/", full.names = T)
dnds_df <- foreach(file_name = dnds_list, .combine = c("bind_rows")) %do% {
  fread(file_name)
}

# Remove zeroes
dnds_filt <-  dnds_df %>%
  mutate(ks = ifelse(ks == 0 | ks < 0, 0, ks),
         ka = ifelse(ka == 0, 0, ka)) %>%
  filter(!(ka == 0 | ks == 0)) %>%
  filter(ka != 0 & ks != 0) %>%
  mutate(kaks = ka / ks)

# Get min. kaks
iter_df <- dnds_filt %>%
  distinct(clique_name, is_jump)

min_df <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  row <- iter_df[i, ]
  dnds_filt %>%
    filter(clique_name == row$clique_name,
           is_jump == row$is_jump) %>%
    arrange(kaks) %>%
    head(1)
}

clique_counts <- dnds_filt %>%
  left_join(host_counts) %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- min_df %>%
  select(clique_name, is_jump, patristic_dist, min_kaks = kaks) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

type_list <- c("ssDNA", "dsDNA", "+ssRNA", "-ssRNA")

plots <- foreach(type = type_list) %do% {
  test <- merged_df %>% filter(genome_type == type)
  
  formula_string <- "n_hosts ~ log(n_genomes) + family + log(min_kaks)"
  jump_mod <- glm(as.formula(formula_string),
                  data = test %>% filter(is_jump),
                  family = poisson())
  
  nonjump_mod <- glm(as.formula(formula_string),
                     data = test %>% filter(!is_jump),
                     family = poisson())
  
  slope_jump <- signif(coef(jump_mod)[["log(min_kaks)"]], 3)
  z_jump <- signif(summary(jump_mod)$coefficient["log(min_kaks)", "z value"], 3)
  df_jump <- nrow(test %>% filter(is_jump)) - nrow(coef(summary(jump_mod)))
  p_jump <- signif(summary(jump_mod)$coefficient["log(min_kaks)", "Pr(>|z|)"], 3)
  
  slope_nonjump <- signif(coef(nonjump_mod)[["log(min_kaks)"]], 3)
  z_nonjump <- signif(summary(nonjump_mod)$coefficient["log(min_kaks)", "z value"], 3)
  df_nonjump <- nrow(test %>% filter(!is_jump)) - nrow(coef(summary(nonjump_mod)))
  p_nonjump <- signif(summary(nonjump_mod)$coefficient["log(min_kaks)", "Pr(>|z|)"], 3)
  
  test %>%
    ggplot(aes(x = n_hosts, y = log10(min_kaks), fill = is_jump, color = is_jump)) +
    geom_point(alpha = 0.9, color = "black", pch = 21) +
    geom_smooth(method = "lm", fill = "grey") +
    scale_fill_manual(values = c("steelblue3", "indianred3")) +
    scale_color_manual(values = c("steelblue3", "indianred3")) +
    labs(x = "Host range (No. host genera)", y = "log10(min dN/dS)",
         color = "Is host jump?",
         title = str_glue("{type}\nHost jumps: slope={slope_jump}, Z={z_jump}, d.f.={df_jump}, p={p_jump}\nNon-host jumps: slope={slope_nonjump}, Z={z_nonjump}, d.f.={df_nonjump}, p={p_nonjump}")) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none", 
          axis.title = element_text(face = "bold"))
}


plt <- egg::ggarrange(plots = plots,
                      )

ggsave("results/dnds_out/dnds_threshold.by_groups.pdf", 
       width = 14, height = 5, dpi = 600,
       plot = plt)

dn <- seq(0, 0.5, 0.01)

