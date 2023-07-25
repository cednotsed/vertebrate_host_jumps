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

clique_list <- list.files("results/dnds_out/agat_fix")
to_do <- clique_list[!(clique_list %in% circular_cliques)]

for(clique_name in to_do) {
  # clique_name = to_do[1]
  print(clique_name)
  gff_list <- list.files(str_glue("results/dnds_out/prokka_out/{clique_name}"),
                         ".gff")
  
  save_dir <- str_glue("results/dnds_out/kill_lists/{clique_name}")
  dir.create(save_dir)
  
  for(gff_name in gff_list) {
    # gff_name = gff_list[1]
    # Get features
    gff <- read.csv(str_glue("results/dnds_out/prokka_out/{clique_name}/{gff_name}"),
                    sep = "\t",
                    header = F,
                    comment.char = "#") %>%
      filter(!is.na(V5)) %>%
      separate(V9, c("id"), ";") %>%
      mutate(id = gsub("ID=", "", id))
    
    # Get actual genome length before concat
    fna_name <- unique(gff$V1)
    fna_filt <- fna[fna_name]
    gen_length <- width(fna_filt)
    
    # Identify features with start pos > actual length
    to_kill <- gff %>%
      filter(V4 > gen_length) %>%
      select(id)
    
    fwrite(to_kill, str_glue("{save_dir}/{fna_name}.txt"),
           col.names = F,
           eol = "\n")
  }
}
