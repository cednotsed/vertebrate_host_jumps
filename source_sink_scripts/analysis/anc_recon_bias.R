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

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
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

zoo_df <- jump_df %>%
  mutate(anc_state = ifelse(anc_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", anc_state),
         tip_state = ifelse(tip_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", tip_state)) %>%
  filter(anc_state == "Homo sapiens" | tip_state == "Homo sapiens") %>%
  mutate(event_type = ifelse(anc_state != "Homo sapiens", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  inner_join(host_meta)

iter_df <- zoo_df %>% distinct(clique_name, anc_state, tip_state, event_type)

zoo_filt <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  # i = 1
  row <- iter_df[i, ]
  anc <- row$anc_state
  tip <- row$tip_state
  clique <- row$clique_name
  
  temp_filt <- zoo_df %>%
    filter(clique_name == clique,
           anc_state == anc,
           tip_state == tip) %>%
    distinct(anc_name, .keep_all = T)
  
  return(temp_filt)
}

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

anthro_df <- zoo_filt %>%
  group_by(clique_name) %>%
  summarise(is_anthro = sum(event_type == "Anthroponotic") > 0)

bias_df <- zoo_filt %>% 
  select(host, event_type, clique_name) %>%
  left_join(count_df) %>%
  left_join(anthro_df) %>%
  separate(clique_name, c("family"), "_", remove = F)

logreg <- glm(event_type ~ is_anthro + log10(genome_ratio),
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
  ggplot(aes(x = event_type, y = genome_ratio, fill = event_type)) +
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

ggsave("results/source_sink_analysis/bias_in_ancestral_reconstruction.distinct_nodes.pdf", width = 5, height = 3)
