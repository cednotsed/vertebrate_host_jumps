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

domestic <- c("Homo", "Sus", "Gallus", 
              "Equus", "Bos", "Anas",
              "Canis", "Felis")
jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  filter(anc_genus %in% domestic &
           tip_genus %in% domestic)

zoo_df <- jump_df %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state))

# Allow only one hostA-hostB jump per node
zoo_filt <- zoo_df %>%
  distinct(anc_state, tip_state, anc_name, clique_name, .keep_all = T)

zoo_filt %>%
  filter(clique_name == "Coronaviridae_12") %>%
  filter(event_type == "Zoonotic")

zoo_filt %>%
  group_by(family, event_type) %>%
  summarise(n = n())
# Observed zoonotic proportion
obs_df <- zoo_filt %>%
  group_by(event_type) %>%
  summarise(n_jumps = n()) %>%
  mutate(prop = n_jumps / sum(n_jumps))

# Bootstrap analysis
cl <- makeCluster(12)
registerDoParallel(cl)

set.seed(66)

boot_morsels <- foreach(i = seq(1000), .packages = c("foreach", "tidyverse")) %dopar% {
  # Resample all jumps
  temp_df <- zoo_filt %>%
    sample_n(nrow(zoo_filt), replace = T) %>%
    group_by(event_type) %>%
    summarise(n_jumps = n()) %>%
    pivot_wider(names_from = "event_type", values_from = "n_jumps")
  
  return(temp_df)
}

stopCluster(cl)

zoo_boot_df <- bind_rows(boot_morsels)

zoo_plot_df <- zoo_boot_df %>%
  pivot_longer(everything(), 
               names_to = "event_type", 
               values_to = "n_jumps")

CI_df <- zoo_plot_df %>%
  group_by(event_type) %>%
  summarise(low_bound = quantile(n_jumps, probs = 0.025), 
            high_bound = quantile(n_jumps, probs = 0.975))

# Add t test results
t_test <- t.test(zoo_boot_df$Anthroponotic , zoo_boot_df$Zoonotic,
                 paired = T, 
                 alternative = "two.sided")

t_val <- signif(t_test$statistic, 3)
deg_freedom <- t_test$parameter
p_val<- signif(t_test$p.value, 3)

zoo_plot_df %>%
  ggplot(aes(x = event_type, y = n_jumps, fill = event_type)) +
  geom_violin(alpha = 0.5) +
  geom_point(data = obs_df, 
             aes(x = event_type),
             size = 3) +
  geom_segment(data = CI_df,
               aes(x = event_type, xend = event_type, 
                   y = low_bound, yend = high_bound)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold")) +
  labs(x = "Transmission type", y = "No. of host jumps",
       title = str_glue("Paired t-test: t = {t_val}, d.f.={deg_freedom}, p={p_val}"))

ggsave("results/source_sink_analysis/anthroponotic_frequency.domestic_only.pdf", 
       dpi = 600, 
       width = 4,
       height = 3)


# Explore jumps by clique
# zoo_meta <- meta %>%
#   filter(accession %in% unique(zoo_df$tip_name))

# # Get viral species to clique mapping
# species_map_df <- foreach(clique = unique(zoo_df$clique_name), .combine = "bind_rows") %do% {
#   temp <- zoo_meta %>%
#     filter(cluster == clique) %>%
#     filter(species != "")
#   
#   tibble(clique_name = clique,
#          species_string = paste0(unique(temp$species), collapse = "; "))
# }
# 
# species_map_df %>% View()
#   # distinct(clique_name) %>%
#   View()
#   nrow()

# No. of distinct animals

## Proportion of SARS-CoV-2 jumps ##
# All jumps
zoo_filt %>%
  group_by(event_type) %>%
  summarise(n = n())

# CoV-2
zoo_filt %>%
  filter(clique_name == "Coronaviridae_12") %>%
  group_by(event_type) %>%
  summarise(n = n()) %>%
  mutate(prop = n / 383)

# MERS
zoo_filt %>%
  filter(clique_name == "Coronaviridae_26") %>%
  group_by(event_type) %>%
  summarise(n = n()) %>%
  mutate(prop = n / 383)

zoo_filt %>%
  filter(clique_name %in% c("Orthomyxoviridae_1", "Orthomyxoviridae_2")) %>%
  group_by(event_type) %>%
  summarise(n = n()) %>%
  mutate(prop = n / 383)

zoo_filt %>%
  filter(!(clique_name %in% c("Coronaviridae_26", "Coronaviridae_12",
                              "Orthomyxoviridae_1", "Orthomyxoviridae_2"))) %>%
  group_by(event_type) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
# zoo_filt %>%
#   filter(clique_name == "Orthomyxoviridae_2") %>%
#   group_by(event_type) %>%
#   summarise(n = n())

## Removing SARS-CoV-2 and Influenza A and MERS ##
# Bootstrap zoonotic proportion
zoo_filt2 <- zoo_filt %>%
  filter(!(clique_name %in% c("Coronaviridae_26", "Coronaviridae_12",
                              "Orthomyxoviridae_1", "Orthomyxoviridae_2")))

obs_df2 <- zoo_filt2 %>%
  group_by(event_type) %>%
  summarise(n_jumps = n()) %>%
  mutate(prop = n_jumps / sum(n_jumps))

cl <- makeCluster(12)
registerDoParallel(cl)

set.seed(66)

boot_morsels <- foreach(i = seq(1000), .packages = c("foreach", "tidyverse")) %dopar% {
  temp_df <- zoo_filt2 %>%
    sample_n(nrow(zoo_filt2), replace = T) %>%
    group_by(event_type) %>%
    summarise(n_jumps = n()) %>%
    pivot_wider(names_from = "event_type", values_from = "n_jumps")
  
  return(temp_df)
}

stopCluster(cl)

zoo_boot_df <- bind_rows(boot_morsels)

zoo_plot_df <- zoo_boot_df %>%
  pivot_longer(everything(), 
               names_to = "event_type", 
               values_to = "n_jumps")

CI_df <- zoo_plot_df %>%
  group_by(event_type) %>%
  summarise(low_bound = quantile(n_jumps, probs = 0.025), 
            high_bound = quantile(n_jumps, probs = 0.975))

# Add t test results
t_test <- t.test(zoo_boot_df$Anthroponotic , zoo_boot_df$Zoonotic,
                 paired = T, 
                 alternative = "two.sided")

t_val <- signif(t_test$statistic, 3)
deg_freedom <- t_test$parameter
p_val<- signif(t_test$p.value, 3)

zoo_plot_df %>%
  ggplot(aes(x = event_type, y = n_jumps, fill = event_type)) +
  geom_violin(alpha = 0.5) +
  geom_point(data = obs_df2, 
             aes(x = event_type),
             size = 3) +
  geom_segment(data = CI_df,
               aes(x = event_type, xend = event_type, 
                   y = low_bound, yend = high_bound)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold")) +
  labs(x = "Transmission type", y = "No. of host jumps",
       title = str_glue("Paired t-test: t = {t_val}, d.f.={deg_freedom}, p={p_val}"))

ggsave("results/source_sink_analysis/anthroponotic_frequency.remove_viruses.pdf",
       dpi = 600, 
       width = 4,
       height = 3)
# ggsave("results/source_sink_analysis/anthroponotic_frequency.pdf", 
#        dpi = 600, 
#        width = 4,
#        height = 3)
# zoo_filt %>%
#   group_by(event_type) %>%
#   group_by(event_type) %>%
#   summarise(n = n()) %>%
#   left_join(species_map_df) %>%
#   mutate(species_string = str_glue("{clique_name}\n({species_string})")) %>%
#   group_by(clique_name, species_string) %>%
#   summarise(n_anthro = sum(event_type == "Anthroponotic"),
#             n_zoo = sum(event_type == "Zoonotic")) %>%
#   arrange(desc(n_zoo)) %>%
#   ggplot(aes(x = n_zoo, y = n_anthro)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_classic() +
#   geom_text_repel(aes(label = clique_name))
# 
# zoo_filt %>%
#   filter(clique_name == "Coronaviridae_26") %>%
#   View()
# 
# zoo_filt %>%
#   group_by(event_type) %>%
#   summarise(n_jumps = n_distinct(anc_name))
# 
# zoo_filt %>%
#   filter(clique_name == "Orthomyxoviridae_1") %>% 
#   group_by(event_type) %>%
#   summarise(n_jumps = n_distinct(anc_name))
# 
# zoo_filt %>%
#   filter(clique_name == "Coronaviridae_26") %>%
#   left_join(species_map_df) %>% View()
#   filter(event_type == "Zoonotic") %>%
#   summarise(n = n())
# 
# zoo_filt %>%
#   group_by(clique_name) %>%
#   summarise(n_anthro = sum(event_type == "Anthroponotic"),
#             n_zoo = sum(event_type == "Zoonotic")) %>%
#   arrange(desc(n_anthro)) %>%
#   ggplot(aes(x = n_zoo, y = n_anthro)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_classic() +
#   geom_text_repel(aes(label = clique_name)) +
#   labs(x = "Zoonotic jumps", y = "Anthroponotic jumps")
# 
# ggsave("results/source_sink_analysis/relative_frequency.pdf", dpi = 600, width = 6, height = 5)
# 
# 
# zoo_filt %>%
#   group_by(clique_name, anc_state, tip_state, event_type) %>%
#   summarise(n_events = n_distinct(anc_name)) %>%
#   filter(n_events >= 1) %>%
#   ungroup() %>%
#   summarise(n_anthro = sum(event_type == "Anthroponotic"),
#             n_zoo = sum(event_type == "Zoonotic"))
# 
# zoo_df
# zoo_filt %>%
#   group_by(distinct)
#   filter(clique_name == "Coronaviridae_4") %>%
#   filter(event_type == "Anthroponotic") %>%
#   select(accession = tip_name, anc_state, tip_state) %>%
#   left_join(meta) %>%
#   # select(genbank_title) %>%
#   View()
# 
#  
# meta %>%
#   filter(accession == "HG994852.1")
# # Check bias
# genome_counts <- genome_meta %>%
#   filter(host != "Homo sapiens") %>%
#   group_by(clique_name, host) %>%
#   summarise(n_genomes = n()) %>%
#   ungroup()
# 
# human_counts <- genome_meta %>%
#   group_by(clique_name) %>%
#   summarise(human_genomes = sum(host == "Homo sapiens")) %>%
#   ungroup()
# 
# count_df <- genome_counts %>%
#   left_join(human_counts) %>%
#   mutate(genome_ratio = n_genomes / human_genomes) %>%
#   filter(human_genomes != 0)
# 
# bias_df <- zoo_filt %>% 
#   mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
#   select(host, event_type, clique_name) %>%
#   left_join(count_df)
# # filter(n_genomes > 20)
# 
# logreg <- glm(event_type ~ genome_ratio,
#               data = bias_df %>% 
#                 mutate(event_type = ifelse(event_type == "Anthroponotic", T, F)) %>%
#                 filter(), 
#               family = "binomial")
# summary(logreg)
# 
# coeffs <- summary(logreg)$coefficients
# slope <- signif(coeffs["genome_ratio", "Estimate"], 3)
# or <- signif(exp(slope), 3)
# z_val <- signif(coeffs["genome_ratio", "z value"], 3)
# p_val <- signif(coeffs["genome_ratio", "Pr(>|z|)"], 3)
# deg_freedom <- nrow(bias_df) - 2
# 
# bias_df %>% 
#   ggplot(aes(x = event_type, y = log10(genome_ratio), fill = event_type)) +
#   labs(x = "Event type", y = "Log10(animal:human genome count)",
#        title = str_glue("Logistic regression: OR(Anthroponotic) = {or}; two-sided z-test: z = {z_val}, d.f. = {deg_freedom}, p = {p_val}")) +
#   theme_classic() +
#   geom_boxplot(position = position_nudge(x = 0.05, y = 0),
#                width = 0.1,
#                outlier.shape = NA,
#                alpha = 0.3) +
#   geom_violinhalf(position = position_nudge(x = 0.1, y = 0), alpha = 1) +
#   theme_classic() +
#   coord_flip() +
#   theme(legend.position = "none")
# 
# # ggsave("results/source_sink_analysis/bias_in_ancestral_reconstruction.pdf", width = 5, height = 3)
