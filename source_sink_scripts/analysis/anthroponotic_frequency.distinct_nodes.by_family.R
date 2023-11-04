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

zoo_df
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

plot_df <- zoo_filt %>% 
  separate(clique_name, c("family"), "_", remove = F) %>%
  group_by(family) %>%
  summarise(n_anthro = sum(event_type == "Anthroponotic"),
            n_zoo = sum(event_type == "Zoonotic"),
            prop_anthro = sum(event_type == "Anthroponotic") / n()) %>%
  pivot_longer(!c(prop_anthro, family), names_to = "event_type", values_to = "count") %>%
  arrange(desc(prop_anthro))


plot_df %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family, y = count, fill = event_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
