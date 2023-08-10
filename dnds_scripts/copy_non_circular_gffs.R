rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(Hmisc)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

circular_cliques <- deframe(meta %>% 
                              filter(is_circular) %>%
                              distinct(cluster))

clique_list <- list.files("results/dnds_out/prokka_out")
to_do <- clique_list[!(clique_list %in% circular_cliques)]
print(length(to_do))

for(clique_name in to_do) {
  # clique_name = to_do[1]
  print(clique_name)
  input_dir <- str_glue("results/dnds_out/agat_fix/{clique_name}")
  gffs <- list.files(input_dir, ".gff", recursive = T, full.names = T)
  save_dir <- str_glue("results/dnds_out/agat_filt/{clique_name}")
  
  dir.create(save_dir)
  file.copy(gffs, save_dir, recursive = T)
  
}
