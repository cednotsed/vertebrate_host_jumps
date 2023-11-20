rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(see)
require(ggpubr)
require(igraph)

genome_type <- fread("data/metadata/genome_type_metadata.csv")

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  dplyr::rename(clique_name = cluster) %>%
  left_join(genome_type)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  filter(anc_genus != "" & tip_genus != "") %>%
  filter(anc_genus != tip_genus)

# Get randomly chosen nodes
merged_filt <- jump_df %>% 
  distinct(clique_name, anc_name, anc_state, tip_state, .keep_all = T)

el <- merged_filt %>%
  select(anc_genus, tip_genus, clique_name) %>%
  group_by(anc_genus, tip_genus) %>%
  summarise(n = n_distinct(clique_name))

g <- graph_from_data_frame(el)

# Add community and host metadata
meta.match <- tibble(host_genus = names(V(g))) %>%
  left_join(meta %>% distinct(host_genus, host_order, host_class))

V(g)$host_order <- meta.match$host_order
V(g)$host_class <- meta.match$host_class

csvs <- igraph::as_data_frame(g, what = "both")

merged <- csvs$vertices %>%
  dplyr::rename(from = name) %>% 
  left_join(csvs$edges)
merged %>% View()
fwrite(merged, "results/source_sink_analysis/host_jump_edgelist.csv")