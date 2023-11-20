rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(doParallel)

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

jump_dir <- "results/source_sink_analysis/permutation/permutation_ancestral_reconstructions/temp_results"
nonjump_dir <- "results/source_sink_analysis/permutation/permutation_ancestral_reconstructions.same_host/temp_results"

all_jump_list <- list.files(jump_dir, full.names = T)
all_nonjump_list <- list.files(nonjump_dir, full.names = T)

all_jump_idx <- str_split(all_jump_list, "\\.", simplify = T)[, 2]
all_nonjump_idx <- str_split(all_nonjump_list, "\\.", simplify = T)[, 3]
length(all_jump_idx)

foreach(idx = seq(100)) %do% {
  jump_list <- all_jump_list[all_jump_idx == idx]
  nonjump_list <- gsub(jump_dir, nonjump_dir, jump_list)
  
  # Check that all files are present
  if (length(jump_list) == length(nonjump_list) &
      length(jump_list) == 248) {
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
    nonjump_morsels <- foreach(file_name = nonjump_list) %do% {
      temp <- fread(file_name)
      if(nrow(temp) > 0) {
        return(temp)
      } 
    }
    
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
    
    fwrite(merged_df, 
           str_glue("results/mutational_load_out/host_jump_lists/permutations/jump_non_jump.diff_genus.{idx}.csv"))
  }
}

merged_df %>%
  filter(is_jump) %>%
  mutate(direction = case_when(anc_genus == "Homo" & tip_genus != "" ~ "Anthroponotic",
                               anc_genus != "" & tip_genus == "Homo"~ "Zoonotic")) %>% 
  filter(!is.na(direction)) %>%
  group_by(direction) %>%
  summarise(n = n())

meta %>%
  group_by(cluster) %>%
  summarise(prop_human = sum(host_genus == "Homo") / n()) %>%
  arrange(desc(prop_human)) %>%
  filter(prop_human != 1,
         prop_human != 0) %>% 
  ggplot(aes(prop_human)) +
  geom_histogram()
