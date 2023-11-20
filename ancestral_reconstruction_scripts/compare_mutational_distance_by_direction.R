rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(see)
require(ggpubr)

genome_type <- fread("data/metadata/genome_type_metadata.csv")

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  dplyr::rename(clique_name = cluster) %>%
  left_join(genome_type)

genome_counts <- meta %>%
  group_by(clique_name) %>%
  summarise(n_genomes = n())

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump)

# Get randomly chosen nodes
iter_df <- jump_df %>% 
  distinct(clique_name, anc_name, anc_state, tip_state)

set.seed(69)

merged_filt <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  row <- iter_df[i,]
  temp <- jump_df %>%
    filter(clique_name == row$clique_name,
           anc_name == row$anc_name,
           anc_state == row$anc_state,
           tip_state == row$tip_state)
  
  temp %>%
    sample_n(1, replace = F)
}

# Get unique host pairs
host_list <- merged_filt %>%
  distinct(anc_state, tip_state)

pair_list <- list()

host_morsels <- foreach(i = seq(nrow(host_list))) %do% {
  host_pair <- as.character(host_list[i, ])
  host_pair <- host_pair[order(host_pair)]
  return(host_pair)
}

pair_filt <- unique(host_morsels)

mutation_df <- foreach(i = seq(length(pair_filt)), .combine = "bind_rows") %do% {
  row <- pair_filt[[i]]
  host1 <- row[1]
  host2 <- row[2]
  
  temp <- merged_filt %>%
    filter(anc_state == host1 & tip_state == host2|
             anc_state == host2 & tip_state == host1) %>%
    mutate(direction = ifelse(anc_state == host1, "Foward", "Reverse"))
  
  if(n_distinct(temp$direction) == 2) {
    return(temp)
  }
}

mutation_parsed <- mutation_df %>%
  left_join(genome_counts)

lr1 <- lm(log(patristic_dist) ~ clique_name + direction,
   data = mutation_parsed)

sink("results/ancestral_reconstruction_out/compare_forward_reverse.mutation.ANOVA.txt")
anova(lr1)
sink()

## Zoonotic jumps only ##
zoo_filt <- merged_filt %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  left_join(genome_counts)

lr2 <- lm(log(patristic_dist) ~ clique_name + event_type,
          data = zoo_filt)

sink("results/ancestral_reconstruction_out/compare_anthro_zoo.mutation.ANOVA.txt")
anova(lr2)
sink()

# mutation_df %>%
#   ggplot(aes(x = direction, y = log10(patristic_dist), fill = direction)) +
#   geom_boxplot() +
#   geom_pwc() +
#   theme_bw() +
#   scale_fill_manual(values = c("goldenrod", "olivedrab4")) +
#   labs(x = "Jump direction", y = "Log10(mutational dist.)") +
#   theme(legend.position = "none")

# logreg <- glm((as.numeric(factor(event_type)) - 1) ~ family + log(kaks),
#               data = zoo_dnds_filt,
#               family = "binomial")
# 
# sink("results/ancestral_reconstruction_out/compare_anthro_zoo.dnds.logreg.txt")
# summary(logreg)
# exp(coef(logreg)[["log(kaks)"]])
# sink()

meta %>%
  filter(clique_name == "Coronaviridae_12") %>%
  filter(host_genus != "") %>%
  summarise(n_distinct(host_genus))
  distinct(species)
