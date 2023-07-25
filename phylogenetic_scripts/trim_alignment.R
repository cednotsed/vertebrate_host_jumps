rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

# aln_paths <- list.files("data/alignments/source_sink_mini_trees/with_buffer_outgroup", 
#                         "\\.aln", 
#                         full.names = T)

aln_paths <- list.files("data/alignments/source_sink_mini_trees/without_outgroup", 
                        "\\.aln", 
                        full.names = T)
aln_paths <- aln_paths[!grepl("log|insertions|to_align|trim_to", aln_paths)]
# aln_paths <- aln_paths[grepl("Herpesviridae_29", aln_paths)]
aln_paths

for (aln_path in aln_paths) {
  print(aln_path)
  aln <- read.dna(aln_path, 
                  format = "fasta",
                  as.matrix = T)
  n_aln <- nrow(aln)
  
  pos_to_keep <- foreach(i = seq(ncol(aln)), .combine = "c") %do% {
    small_gap <- sum(as.character(aln[, i]) %in% c("-", "n")) / n_aln <= 0.1
    
    if(small_gap) {
      return(i)
    } else {
      return(NA)
    }
  }
  
  pos_to_keep <- pos_to_keep[!is.na(pos_to_keep)]
  trimmed <- aln[, pos_to_keep]
  
  # Parse reference name
  ref_name <- rownames(trimmed)[1]
  parsed_ref_name <- str_split(ref_name, " ")[[1]][1]
  rownames(trimmed)[1] <- parsed_ref_name
  
  # Save trimmed alignment
  save_name <- gsub(".aln", str_glue(".trim_to_{ncol(trimmed)}pos.aln"), aln_path)
  write.FASTA(trimmed, save_name)
  rm(trimmed)
}

