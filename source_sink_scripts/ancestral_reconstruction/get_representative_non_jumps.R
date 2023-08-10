rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

host_meta <- fread("data/metadata/parsed_host_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(host_meta)

jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  as_tibble() %>%
  arrange(clique_name, desc(n_traverses))

res_dir <- "results/source_sink_analysis/ancestral_reconstructions.same_host/temp_results"

file_list <- list.files(res_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  temp <- fread(file_name)
  if(nrow(temp) > 0) {
    return(temp)
  } 
}

res_df <- bind_rows(morsels)
res_df %>%
  fwrite("results/mutational_load_out/non_jumps.csv")
# Before full analysis is completed
available <- unique(res_df$clique_name)
length(available)

jump_df <- jump_df %>% filter(clique_name %in% available)

# Check if all are non-host jumps
res_df %>% 
  filter(anc_state != tip_state) %>%
  nrow()

# Get traverses with at least 50% of same host
max_df <- res_df %>%
  group_by(clique_name, tip_name, tip_state) %>%
  summarise(max_traverses = max(n_traverses))

res_filt <- res_df %>%
  left_join(max_df) %>%
  # mutate(prop_max_traversed = max_traverses) %>%
  filter(max_traverses > 4)

# Match host jumps to non host jumps
match_morsels <- foreach(clique = unique(jump_df$clique_name)) %do% {
  jump_filt <- jump_df %>%
    filter(clique_name == clique) %>%
    distinct(tip_name)
  jump_count <- n_distinct(jump_filt$tip_name)
  
  res_temp <- res_filt %>%
    filter(clique_name == clique) 
  
  nonjump_accs <- deframe(res_temp %>%
    distinct(tip_name))
  
  nonjump_count <- length(nonjump_accs)
  
  if (nonjump_count > jump_count) {
    accs <- sample(nonjump_accs, jump_count, replace = F)
  } else {
    accs <- nonjump_accs
  }
  
  crumbs <- foreach(acc = accs) %do% {
    res_temp %>%
      filter(tip_name == acc) %>% 
      sample_n(1, replace = F)
  }
  
  bind_rows(crumbs)
} 

matched <- bind_rows(match_morsels)
nrow(matched)

matched %>%
  fwrite("results/mutational_load_out/non_jumps.csv")
  # res_filt %>%
  #   filter(clique_name == clique) %>%
  #   distinct(tip_name) %>%
  #   sample_n(ifelse(jump_filt))
  
  # non_jump <- res_filt %>%
  #   filter(!(tip_name %in% tips_used)) %>%
  #   filter(clique_name == clique) %>%
  #   filter(n_traverses %in% seq(ifelse(n_trav <= 1, 1, n_trav - 1), n_trav + 1))
  # 
  # if (nrow(non_jump > 0)) {
  #   rep_temp <- non_jump %>%
  #     sample_n(size = 1, replace = F)
  #   tips_used <- c(tips_used, rep_temp$tip_name)
  #   return(rep_temp)
  # } else {
  #   return(NULL)
  # }


# res_df %>%
#   filter(tip_vertebrate & anc_vertebrate) %>%
#   fwrite("results/source_sink_analysis/putative_host_jumps.csv")