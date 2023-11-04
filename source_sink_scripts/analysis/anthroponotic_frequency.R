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

host_meta <- fread("data/metadata/parsed_host_metadata.csv")
genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  mutate(host = ifelse(grepl("sp\\.", host), 
                       gsub(" sp\\.", "", host), 
                       host)) %>%
  mutate(host = ifelse(host %in% c("Homo", "Homo sapiens"), 
                       "Homo sapiens", 
                       host)) %>%
  filter(host != "") %>%
  dplyr::rename(clique_name = cluster)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name)

# Count distinct nodes
iter_df <- jump_df %>% distinct(clique_name, anc_state, tip_state)

jump_filt <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  # i = 1
  row <- iter_df[i, ]
  anc <- row$anc_state
  tip <- row$tip_state
  clique <- row$clique_name
  
  temp_filt <- jump_df %>%
    filter(clique_name == clique,
           anc_state == anc,
           tip_state == tip) %>%
    distinct(anc_name, .keep_all = T)
  
  return(temp_filt)
}

# No. of distinct host jumps
jump_filt %>%
  nrow()

# Prop. of human jumps
jump_filt %>%
  left_join(host_meta %>% select(anc_state = host, host1_genus = host_genus)) %>%
  left_join(host_meta %>% select(tip_state = host, host2_genus = host_genus)) %>%
  filter(host1_genus == "Homo"| host2_genus == "Homo") %>%
  nrow()

zoo_df <- jump_df %>%
  mutate(anc_state = ifelse(anc_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", anc_state),
         tip_state = ifelse(tip_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", tip_state)) %>%
  filter(anc_state == "Homo sapiens" | tip_state == "Homo sapiens") %>%
  mutate(event_type = ifelse(anc_state != "Homo sapiens", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  inner_join(host_meta)

# Check bias
genome_counts <- genome_meta %>%
  filter(host != "Homo sapiens") %>%
  group_by(clique_name, host) %>%
  summarise(n_genomes = n()) %>%
  ungroup()

human_counts <- genome_meta %>%
  group_by(clique_name) %>%
  summarise(human_genomes = sum(host == "Homo sapiens")) %>%
  ungroup()
  
count_df <- genome_counts %>%
  left_join(human_counts) %>%
  mutate(genome_ratio = n_genomes / human_genomes) %>%
  filter(human_genomes != 0)

bias_df <- zoo_df %>% 
  select(host, event_type, clique_name) %>%
  left_join(count_df)
  # filter(n_genomes > 20)

logreg <- glm(event_type ~ genome_ratio,
              data = bias_df %>% 
                mutate(event_type = ifelse(event_type == "Anthroponotic", T, F)) %>%
                filter(), 
              family = "binomial")
summary(logreg)

coeffs <- summary(logreg)$coefficients
slope <- signif(coeffs["genome_ratio", "Estimate"], 3)
or <- signif(exp(slope), 3)
z_val <- signif(coeffs["genome_ratio", "z value"], 3)
p_val <- signif(coeffs["genome_ratio", "Pr(>|z|)"], 3)
deg_freedom <- nrow(bias_df) - 2

bias_df %>% 
  ggplot(aes(x = event_type, y = log10(genome_ratio), fill = event_type)) +
  labs(x = "Event type", y = "Log10(animal:human genome count)",
       title = str_glue("Logistic regression: OR(Anthroponotic) = {or}; two-sided z-test: z = {z_val}, d.f. = {deg_freedom}, p = {p_val}")) +
  theme_classic() +
  geom_boxplot(position = position_nudge(x = 0.05, y = 0),
               width = 0.1,
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.1, y = 0), alpha = 1) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none")

ggsave("results/source_sink_analysis/bias_in_ancestral_reconstruction.pdf", width = 5, height = 3)

homo_plot_df <- zoo_df %>%
  left_join(count_df, by = c("clique_name", "host")) %>%
  # filter(n_genomes > 1) %>%
  group_by(host_order) %>%
  summarise(prop_anthro = sum(event_type == "Anthroponotic") / n(),
            n_cliques = n_distinct(clique_name)) %>%
  filter(host_order != "")

pal <- randomcoloR::distinctColorPalette(n_distinct(homo_plot_df$host_order))

# Asymmetry plot
homo_plot_df %>%
  ggplot(aes(x = host_order, y = prop_anthro)) +
  geom_point(aes(fill = host_order,
                 size = n_cliques),
             color = "black",
             pch = 21) +
  scale_fill_manual(values = pal, guide = "none") +
  geom_hline(yintercept = 0.5, lty = "dashed", color = "black") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.5,
                                   face = "italic"),
        axis.title = element_text(face = "bold")) +
  labs(x = "Host order", 
       y = "Prop. of anthroponotic jumps",
       size = "No. of viral cliques") 

ggsave("results/source_sink_analysis/anthroponotic_assymetry.by_host_order.png", 
       dpi = 300, 
       height = 6,
       width = 8)

# No. of total anthroponotic jumps
zoo_df %>%
  summarise(prop_anthro = sum(event_type == "Anthroponotic") / n())
# Bootstrap zoonotic proportion
zoo_filt <- zoo_df %>%
  distinct(clique_name, anc_state, tip_state, event_type, .keep_all = T)

zoo_filt
# zoo_df %>%
#   group_by(clique_name) %>%
#   summarise(is_anthro = any(event_type == "Anthroponotic"),
#             is_zoo = any(event_type == "Zoonotic")) %>%
#   ungroup() %>%
#   summarise(prop_anthro = sum(is_anthro) / n(),
#             prop_zoo = sum(is_zoo) / n())

obs_df <- zoo_filt %>%
  group_by(event_type) %>%
  summarise(n_jumps = n()) %>%
  mutate(prop = n_jumps / sum(n_jumps))

obs_df
# 
# clique_list <- unique(zoo_filt$clique_name)
# 
# set.seed(66)
# cl <- makeCluster(12)
# registerDoParallel(cl)
# 
# boot_morsels <- foreach(i = seq(1000),
#                         .packages = c("foreach", "tidyverse")) %dopar% {
#                           clique_temp <- sample(clique_list, length(clique_list), replace = T)
#                           
#                           temp_morsels <- foreach(clique_name = clique_list) %do% {
#                             zoo_filt %>%
#                               filter(clique_name == clique_name)
#                           }
#                           
#                           temp_df <- bind_rows(temp_morsels) %>%
#                             sample_n(nrow(zoo_filt), replace = T) %>%
#                             group_by(event_type) %>%
#                             summarise(n_jumps = n()) %>%
#                             pivot_wider(names_from = "event_type", values_from = "n_jumps")
#                           
#                           return(temp_df)
#                         }
# 
# stopCluster(cl)
# 
# zoo_boot_df <- bind_rows(boot_morsels)
# 
# zoo_plot_df <- zoo_boot_df %>%
#   pivot_longer(everything(), 
#                names_to = "event_type", 
#                values_to = "n_jumps")
# 
# CI_df <- zoo_plot_df %>%
#   group_by(event_type) %>%
#   summarise(low_bound = quantile(n_jumps, probs = 0.025), 
#             high_bound = quantile(n_jumps, probs = 0.975))
# 
# # Add t test results
# t_test <- t.test(zoo_boot_df$Anthroponotic , zoo_boot_df$Zoonotic,
#                  paired = T, 
#                  alternative = "two.sided")
# 
# t_val <- signif(t_test$statistic, 3)
# deg_freedom <- t_test$parameter
# p_val<- signif(t_test$p.value, 3)
# 
# zoo_plot_df %>%
#   ggplot(aes(x = event_type, y = n_jumps, fill = event_type)) +
#   geom_violin(alpha = 0.5) +
#   geom_point(data = obs_df, 
#              aes(x = event_type),
#              size = 3) +
#   geom_segment(data = CI_df,
#                aes(x = event_type, xend = event_type, 
#                    y = low_bound, yend = high_bound)) +
#   theme_classic() +
#   theme(legend.position = "none",
#         axis.title = element_text(face = "bold")) +
#   labs(x = "Transmission type", y = "No. of host jumps",
#        title = str_glue("Paired t-test: t = {t_val}, d.f.={deg_freedom}, p={p_val}"))
# 
# ggsave("results/source_sink_analysis/anthroponotic_frequency.pdf", 
#        dpi = 600, 
#        width = 4,
#        height = 3)

zoo_df %>%
  filter(event_type == "Anthroponotic") %>%
  filter(host_order == "Passeriformes")
genome_meta %>% filter(clique_name == "Flaviviridae_18")
