rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(Hmisc)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
circular_families <- c("Anelloviridae", "Circoviridae", "Genomoviridae",
                       "Hepadnaviridae", "Smacoviridae")

gen_list <- list.files("data/genomes/source_sink_mini_trees/without_outgroup/", 
                       "\\.fna",
                       full.names = T)

for(gen_path in gen_list) {
  clique_name <- str_split(gen_path, "/")[[1]]
  clique_name <- str_split(clique_name[length(clique_name)], "\\.")[[1]][1]
  family_name <- str_split(clique_name, "_")[[1]][1]
  
  print(clique_name)

  fna <- readDNAStringSet(gen_path)
  
  # Concatenate genomes if circular
  if (family_name %in% circular_families) {
    fna_parsed <- DNAStringSet(paste0(as.character(fna), as.character(fna)))
    names(fna_parsed) <- names(fna)
  } else {
    fna_parsed <- fna
  }
  
  # Make clique genome directory
  # dir.create(str_glue("data/genomes/dnds_analysis/{clique_name}"))
  
  # Split fasta
  for(i in seq(length(fna_parsed))) {
    genome <- fna_parsed[i]
    genome_name <- names(genome)
    writeXStringSet(genome, str_glue("data/genomes/dnds_analysis/{clique_name}/{genome_name}.fna"))
  }
}
