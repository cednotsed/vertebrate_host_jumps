rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(data.table)
require(tidyverse)
require(ape)
require(randomcoloR)
require(castor)
require(foreach)
require(adephylo)

genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  dplyr::rename(clique_name = cluster)

host_meta <- genome_meta %>% distinct(host, host_species, host_genus)

genome_counts <- genome_meta %>%
  group_by(host_species) %>%
  summarise(n_genomes = n_distinct(accession)) %>%
  arrange(desc(n_genomes))

total_genome_counts <- genome_meta %>%
  group_by(clique_name) %>%
  summarise(n_genomes = n_distinct(accession)) %>%
  arrange(desc(n_genomes))

host_counts <- genome_meta %>%
  filter(host_species != "") %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_species))

contree <- read.tree("data/supertree/olival_cytb_supertree.tree")
contree$tip.label <- gsub("_", " ", contree$tip.label)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  filter(is_jump) %>%
  left_join(host_meta %>%
              distinct(anc_state = host, anc_species = host_species)) %>%
  left_join(host_meta %>%
              distinct(tip_state = host, tip_species = host_species)) %>%
  filter(anc_species != "" & tip_species != "") %>%
  filter(anc_species %in% contree$tip.label & tip_species %in% contree$tip.label) %>%
  group_by(clique_name, family, anc_species, tip_species) %>%
  summarise(min_dist = min(patristic_dist))

dist_df <- foreach(i = seq(nrow(jump_df)), .combine = "bind_rows") %do% {
  row <- jump_df[i, ]
  row %>% mutate(host_dist = get_pairwise_distances(contree, row$anc_species, row$tip_species))
}

merged_df <- dist_df %>%
  left_join(genome_counts %>% select(anc_species = host_species, anc_genomes = n_genomes)) %>%
  left_join(genome_counts %>% select(tip_species = host_species, tip_genomes = n_genomes)) %>%
  left_join(host_counts) %>%
  left_join(total_genome_counts) %>%
  mutate(mean_effort = 0.5 * (anc_genomes + tip_genomes)) %>%
  mutate(involves_humans = (anc_species == "Homo sapiens" | tip_species == "Homo sapiens")) %>%
  ungroup()

count_df <- merged_df %>%
  group_by(clique_name) %>%
  summarise(n_jumps = n())
  # filter(clique_name != "Coronaviridae_12")
  # filter(n_jumps > 5)

merged_filt <- merged_df %>%
  mutate(involves_humans = (anc_species == "Homo sapiens" | tip_species == "Homo sapiens")) %>%
  filter(clique_name %in% count_df$clique_name)
  
merged_filt %>%
  ggplot(aes(y = log10(min_dist), x = host_dist, color = involves_humans)) +
  geom_point() +
  theme_classic() +
  labs(x = "Min mutational dist.", y = "Host distance",
       fill = "Involves humans") +
  geom_smooth(method = "loess")

lr1 <- lm(log(min_dist) ~ log(anc_genomes) + log(tip_genomes) + family + log(host_dist),
          data = merged_filt %>%
            filter(!involves_humans))
summary(lr1)

lr2 <- lm(log(min_dist) ~ log(anc_genomes) + log(tip_genomes) + family + log(host_dist),
          data = merged_filt %>%
            filter(involves_humans))
summary(lr2)
summary(lr2)

anova(lr2)
hist(lr2$residuals)
merged_df %>%
  filter(log10(min_dist) > -5) %>%
  ggplot(aes(log10(min_dist))) +
  geom_histogram()

hist(lr$residuals)
shapiro.test(lr$residuals)
hist(log(merged_df$tip_genomes))

test <- merged_filt %>%
  left_join(genome_meta %>% distinct(anc_species = host_species, anc_family = host_family)) %>%
  left_join(genome_meta %>% distinct(tip_species = host_species, tip_family = host_family)) %>%
  mutate(cross_family = anc_family == tip_family)


test_lr <- lm(min_dist ~ clique_name + cross_family,
   data = test)
summary(test_lr)
