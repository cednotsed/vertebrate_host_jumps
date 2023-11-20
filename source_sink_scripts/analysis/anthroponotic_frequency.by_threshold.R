rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  rename(clique_name = cluster)

genome_counts <- meta %>%
  group_by(clique_name) %>%
  summarise(n = n_distinct(accession),
            prop_animal = sum(host_genus != "Homo") / n_distinct(accession)) %>%
  arrange(prop_animal)

zoo_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  filter(anc_genus == "Homo"| tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  distinct(clique_name, anc_name, tip_state, .keep_all = T) %>%
  left_join(genome_counts)

# genome_counts %>%
#   filter(prop_human > 0.9)
#   ggplot(aes(x = cluster, y = prop_human)) +
#   geom_bar(stat = "identity")

threshold_df <- foreach(t = seq(0, 1, 0.05), .combine = "bind_rows") %do% {
  zoo_df %>%
    filter(prop_animal > t) %>%
    summarise(n_anthro = sum(event_type == "Anthroponotic"),
              n_zoo = sum(event_type == "Zoonotic"),
              n_cliques = n_distinct(clique_name)) %>%
    mutate(threshold = t)
}

threshold_df %>%
  mutate(ratio = n_anthro / n_zoo) %>%
  filter(threshold < 0.20) %>%
  ggplot(aes(x = threshold, y = ratio, size = n_cliques)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 1, 
             lty = "dashed",
             color = "red") +
  geom_text(aes(label = str_glue("n={n_cliques}")),
            size = 3,
            vjust = 2) +
  labs(x = "Prop. animal genomes > t", y = "Anthroponotic:zoonotic ratio")

ggsave("results/source_sink_analysis/ratio_by_threshold.pdf",
       dpi = 600,
       width = 5, height = 5)
  

meta %>%
  filter(clique_name == "Coronaviridae_12") %>%
  summarise(n = sum(host_genus != "Homo")) 
species_mapping <- meta %>%
  distinct(species, clique_name)
genome_counts %>%
  filter(clique_name %in% zoo_df$clique_name) %>%
  filter(prop_animal < 0.20) %>%
  mutate(prop_animal = round(prop_animal, 2)) %>%
  left_join(species_mapping) %>%
  fwrite("results/source_sink_analysis/low_prop_cliques.csv")

zoo_df %>%
  filter(clique_name == "Herpesviridae_17")
  # distinct(anc_state)
