rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)

fna_dir <- "data/genomes/all_viruses.220723.filt"
fna_list <- list.files(fna_dir, "\\.fna")
acc_list <- gsub("\\.fna", "", fna_list)

meta <- fread("data/metadata/all_viruses.220723.filt.csv")

meta %>%
  filter(!(accession %in% acc_list)) %>% 
  select(accession) %>%
  fwrite("data/metadata/missing_genomes.accessions_only.txt",
         col.names = F,
         eol = "\n")
