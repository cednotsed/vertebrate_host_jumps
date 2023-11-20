rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(data.table)
require(tidyverse)
require(ape)
require(randomcoloR)
require(castor)
require(foreach)
require(adephylo)

node_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/node_distance.parsed_jumps.csv") %>%
  distinct(clique_name, anc1_name, anc2_name, .keep_all = T)
tip_df <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump)

plot_df <- node_df %>% 
  rename(anc_name = anc2_name,
         node_dist = patristic_dist) %>%
  right_join(tip_df %>% rename(tip_dist = patristic_dist)) %>% 
  select(clique_name, anc_name, anc1_name, tip_name, tip_dist, node_dist) %>%
  separate(clique_name, c("family"), "_", remove = F)

plot_df %>%
  ggplot(aes(x = tip_dist, y = node_dist)) +
  geom_point() +
  geom_smooth()

cor.test(plot_df$tip_dist, plot_df$node_dist)

# lr <- lm(tip_dist ~ n_genomes + family + node_dist,
#    data = plot_df)
# 
# summary(lr)
# 
# plot_df %>%
#   filter(family == "Retroviridae") %>%
#   ggplot(aes(x = tip_dist, y = node_dist)) +
#   geom_point() +
#   geom_smooth()
