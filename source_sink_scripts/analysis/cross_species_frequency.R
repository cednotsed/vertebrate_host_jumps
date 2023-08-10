setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(Hmisc)
require(ggpubr)
require(randomcoloR)
require(igraph)
require(ggrepel)

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

jump_df <- fread("results/source_sink_analysis/putative_host_jumps.csv") %>%
  mutate(prop_traverse = n_traverses / total_depth)

# # Get unique jumps for each host pair and direction
# jump_list <- jump_df %>%
#   distinct(clique_name, anc_state, tip_state)
# 
# jump_morsels <- foreach(i = seq(nrow(jump_list))) %do% {
#   # i = 1
#   row <- jump_list[i, ]
#   clique_name <- row$clique_name
#   anc_host <- row$anc_state
#   tip_host <- row$tip_state
#   jump_df %>% 
#     filter(clique_name == clique_name,
#            anc_state == anc_host, 
#            tip_state == tip_host) %>%
#     arrange(prop_traverse) %>%
#     head(1)
# }

# Check bias
genome_counts <- genome_meta %>%
  left_join(host_meta) %>%
  group_by(host_genus) %>%
  summarise(n_genomes = n()) %>%
  ungroup()

edge_list <- jump_df %>%
  # distinct(clique_name, anc_genus, tip_genus) %>%
  # group_by(anc_order, tip_order) %>%
  # summarise(n_cliques = n_distinct(clique_name)) %>%
  filter(anc_genus != "" & tip_genus != "") %>%
  select(anc_genus, tip_genus)

g <- graph_from_data_frame(edge_list, directed = TRUE)

# plot(g, 
#      layout = layout.fruchterman.reingold(g),
#      edge.width = E(g)$n_cliques)
plot_df <- tibble(host_genus = names(degree(g, mode = "out")),
                  indegree = degree(g, mode = "in"),
                  outdegree = degree(g, mode = "out")) %>%
  mutate(error = abs(indegree - outdegree) > 50,
         node_type = case_when(outdegree - indegree > 0 ~ "Source",
                               outdegree - indegree < 0 ~ "Sink",
                               outdegree == indegree ~ "Neutral")) %>%
  left_join(genome_counts)
  
plot_filt <- plot_df %>% filter(error)

plot_df %>%
  ggplot(aes(x = indegree, y = outdegree, 
             size = n_genomes, 
             color = node_type)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1,
                lty = "dashed",
                color = "grey40") +
    geom_text_repel(mapping = aes(x = indegree, y = outdegree, label = host_genus),
                    data = plot_filt,
                    size = 3,
                    color = "black") +
  scale_color_manual(values = c("black", "steelblue4", "salmon")) +
  theme_classic() +
  # xlim(0, 60) +
  labs(x = "Viral cliques received", y = "Viral cliques given",
       color = "Host type", 
       size = "No. genomes")

ggsave("results/source_sink_analysis/node_degrees.all_jump.png", width = 6, height = 5)
