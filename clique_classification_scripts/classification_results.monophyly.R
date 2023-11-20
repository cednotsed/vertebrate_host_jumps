rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  mutate(species = tolower(species)) 

tax_meta <- fread("data/VIRION.v0.2.1_240922/TaxonomyVirus.csv") %>%
  # filter(ICTVRatified) %>%
  select(species = Virus)

file_dir <- "results/clique_classification_out/monophyly/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)
}

bind_rows(morsels) %>%
  group_by(t) %>%
  summarise(prop_monophyly = sum(mono) / sum(total)) %>%
  ggplot(aes(x = t, y = prop_monophyly)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  labs(x = "Mash distance threshold", y = "Prop. monophyletic cliques") +
  geom_vline(xintercept = 0.15,
             lty = "dashed",
             color = "red")

ggsave("results/clique_classification_out/monophyly_results.pdf", width = 5, height = 3)


# bind_rows(morsels) %>%
#   filter(t == 0.15) %>%
#   arrange(prop_mono) %>%
#   group_by(t) %>%
#   summarise(median_monophyly = median(prop_mono, na.rm = T)) %>%
#   ggplot(aes(x = t, y = median_monophyly, color = family)) +
#   geom_point() +
#   geom_line() +
#   theme_classic() +
#   labs(x = "Mash distance threshold", y = "Median prop. monophyletic clades") +
#   geom_vline(xintercept = 0.15,
#              lty = "dashed",
#              color = "red")

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)
}

bind_rows(morsels) %>%
  summarise(median_mono = median(prop_mono, na.rm = T))
