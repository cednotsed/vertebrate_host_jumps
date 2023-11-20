rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)

host_meta <- fread("data/metadata/parsed_host_metadata.csv")
genome_type <- fread("data/metadata/genome_type_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv") %>%
  left_join(host_meta) %>%
  left_join(genome_type)

good_alns <- fread("results/qc_out/good_alignments.csv")

# Get host counts
host_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  filter(host_genus != "") %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_genus))

genome_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(n_genomes = n_distinct(accession))

jump_dir <- "results/source_sink_analysis/ancestral_reconstructions/temp_results"
nonjump_dir <- "results/source_sink_analysis/ancestral_reconstructions.same_host/temp_results"

jump_list <- list.files(jump_dir, full.names = T)

# Parse files
jump_morsels <- foreach(file_name = jump_list) %do% {
  temp <- fread(file_name)
  if(nrow(temp) > 0) {
    return(temp)
  } 
}

jump_df <- bind_rows(jump_morsels) %>%
  left_join(host_meta %>%
              select(anc_state = host,
                     anc_genus = host_genus,
                     anc_family = host_family,
                     anc_order = host_order,
                     anc_vertebrate = is_vertebrate)) %>%
  left_join(host_meta %>%
              select(tip_state = host,
                     tip_genus = host_genus,
                     tip_family = host_family,
                     tip_order = host_order,
                     tip_vertebrate = is_vertebrate)) %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  mutate(is_jump = T)

# Get matching non jumps
nonjump_list <- list.files(nonjump_dir, full.names = T)

nonjump_morsels <- foreach(file_name = nonjump_list) %do% {
  temp <- fread(file_name)
  if(nrow(temp) > 0) {
    return(temp)
  } 
}

length(jump_list) == length(nonjump_list)

nonjump_df <- bind_rows(nonjump_morsels) %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  # Remove tips that are involved in host jumps
  filter(!(tip_name %in% jump_df$tip_name)) %>%
  left_join(host_meta %>%
              select(anc_state = host,
                     anc_genus = host_genus,
                     anc_family = host_family,
                     anc_order = host_order,
                     anc_vertebrate = is_vertebrate)) %>%
  left_join(host_meta %>%
              select(tip_state = host,
                     tip_genus = host_genus,
                     tip_family = host_family,
                     tip_order = host_order,
                     tip_vertebrate = is_vertebrate)) %>%
  filter(clique_name %in% good_alns$clique_name)
  # filter(tip_vertebrate & anc_vertebrate)

# Get matching non-jumps
cl <- makeCluster(12)
registerDoParallel(cl)

tip_list <- deframe(nonjump_df %>%
                      distinct(tip_name))

set.seed(69)

subsample_morsels <- foreach(tip = tip_list, .packages = c("tidyverse")) %dopar% {
  nonjump_temp <- nonjump_df %>%
    filter(tip_name == tip) %>%
    sample_n(1, replace = F)
}

stopCluster(cl)

nonjump_filt <- bind_rows(subsample_morsels) %>%
  mutate(is_jump = F)

# Merge non-jumps and jumps
merged_df <- bind_rows(nonjump_filt, jump_df) %>%
  left_join(host_counts) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  # filter(clique_name %in% nonjump_filt$clique_name & 
  #          clique_name %in% jump_df$clique_name) %>%
  filter(tip_vertebrate & anc_vertebrate)

# # Check if all cliques have jumps and nonjumps
# type_counts <- merged_df %>% 
#   group_by(clique_name) %>% 
#   summarise(n = n_distinct(is_jump)) %>%
#   arrange(n) %>%
#   filter(n > 1)
# 
# merged_filt <- merged_df %>%
#   filter(clique_name %in% type_counts$clique_name)
# 
# table(merged_filt$is_jump)

fwrite(merged_df, "results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.V2.csv")

# res_df %>%
#   filter(is.na(host_order1)|is.na(host_order2)) %>% 
#   distinct(anc_state, tip_state, host_order1, host_order2) %>% 
#   View()


