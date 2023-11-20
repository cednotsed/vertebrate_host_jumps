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

dnds_list <- list.files("results/dnds_out/all_jumps.temp_results/", full.names = T)
dnds_df <- foreach(file_name = dnds_list, .combine = c("bind_rows")) %do% {
  fread(file_name)
}

# Remove zeroes
dnds_filt <-  dnds_df %>%
  mutate(ks = ifelse(ks == 0 | ks < 0, 0, ks),
         ka = ifelse(ka == 0, 0, ka)) %>%
  # filter(!(ka == 0 | ks == 0)) %>%
  filter(ka != 0 & ks != 0) %>%
  mutate(kaks = ka / ks)

jump_df <- dnds_filt %>%
  filter(is_jump) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(host_meta %>%
              distinct(anc_state = host, anc_species = host_species)) %>%
  left_join(host_meta %>%
              distinct(tip_state = host, tip_species = host_species))

parsed <- jump_df %>%
  filter(anc_species != "" & tip_species != "") %>%
  filter(anc_species %in% contree$tip.label & tip_species %in% contree$tip.label) %>%
  group_by(clique_name, family, anc_species, tip_species) %>%
  summarise(min_kaks = min(kaks))

dist_df <- foreach(i = seq(nrow(parsed)), .combine = "bind_rows") %do% {
  row <- parsed[i, ]
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

merged_df %>%
  mutate(involves_humans = (anc_species == "Homo sapiens" | tip_species == "Homo sapiens")) %>%
  ggplot(aes(x = host_dist, y = log(min_kaks), color = involves_humans)) +
  geom_point() +
  theme_classic() +
  labs(x = "log(min. dN/dS)", y = "host distance",
       fill = "Involves humans") +
  geom_smooth(method = "lm")


lr1 <- lm(log(min_kaks) ~ log(anc_genomes) + log(tip_genomes) + family + log(host_dist),
          data = merged_df %>%
            filter(!involves_humans))

summary(lr1)

lr2 <- lm(log(min_kaks) ~ log(anc_genomes) + log(tip_genomes) + family + log(host_dist),
          data = merged_df %>%
            filter(involves_humans))

summary(lr2)

merged_df %>%
  filter(log10(min_dist) > -5) %>%
  ggplot(aes(log10(min_dist))) +
  geom_histogram()

hist(lr$residuals)
shapiro.test(lr$residuals)
hist(log(merged_df$tip_genomes))
