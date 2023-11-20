rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)
require(Biostrings)

genome_type <- fread("data/metadata/genome_type_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")

dnds_list <- list.files("results/dnds_out/all_jumps.by_gene.temp_results/", full.names = T)

dnds_df <- foreach(file_name = dnds_list, .combine = "bind_rows") %do% {
  if(file.size(file_name) > 246) {
    fread(file_name)
  }
}
