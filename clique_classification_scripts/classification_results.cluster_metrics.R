rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  mutate(species = tolower(species)) 

file_dir <- "results/clique_classification_out/clustering_metrics/"
file_list <- list.files(file_dir, full.names = T)
file_list <- file_list[grepl("ICTV", file_list)]

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)
}

bind_rows(morsels) %>%
  # filter(type == "species") %>%
  group_by(t, type) %>%
  summarise(median_ari = median(ari, na.rm = T), 
            median_ami = median(ami, na.rm = T)) %>%
  pivot_longer(!c(t, type), names_to = "metric", values_to = "value") %>%
  mutate(metric = ifelse(metric == "median_ami", "Median AMI", "Median ARI")) %>%
  ggplot(aes(x = t, y = value, color = type, shape = metric)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0.15,
             color = "red",
             lty = "dashed") +
  theme_classic() +
  labs(x = "Mash distance threshold",
       y = "Value",
       shape = "Metric",
       color = "Taxon level")

ggsave("results/clique_classification_out/clustering_metrics.ICTV.pdf", width = 5, height = 3)

