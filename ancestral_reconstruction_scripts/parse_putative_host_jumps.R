rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)

host_meta <- fread("data/metadata/parsed_host_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv") %>%
  left_join(host_meta)

res_dir <- "results/source_sink_analysis/ancestral_reconstructions/temp_results"

file_list <- list.files(res_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  temp <- fread(file_name)
  if(nrow(temp) > 0) {
    return(temp)
  } 
}

res_df <- bind_rows(morsels) %>%
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
                     tip_vertebrate = is_vertebrate))

res_df %>%
  filter(tip_vertebrate & anc_vertebrate) %>%
  fwrite("results/mutational_load_out/putative_host_jumps.csv")
# res_df %>%
#   filter(is.na(host_order1)|is.na(host_order2)) %>% 
#   distinct(anc_state, tip_state, host_order1, host_order2) %>% 
#   View()


