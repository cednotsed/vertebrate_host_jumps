rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(ape)
require(foreach)

# aln_paths <- list.files("data/alignments/source_sink_mini_trees/with_buffer_outgroup", 
#                         "\\.aln", 
#                         full.names = T)

aln_paths <- list.files("data/alignments/source_sink_mini_trees/without_outgroup", 
                        "\\.aln", 
                        full.names = T)
aln_paths <- aln_paths[!grepl("log|insertions|to_align|trim_to|masked", aln_paths)]
# aln_paths <- aln_paths[grepl("Herpesviridae_29", aln_paths)]
aln_paths

for (aln_path in aln_paths) {
  # aln_path = aln_paths[1]
  print(aln_path)
  aln <- readDNAStringSet(aln_path)
  aln_mat <- as.matrix(aln)
  
  # Get gapped positions
  gappy_pos <- apply(aln_mat, 2,
                     function(x) {sum(x %in% c("-", "N")) / length(x) > 0.1})
  n_unmasked <- sum(!gappy_pos)
  
  # Mask gapped positions
  aln_mat[, gappy_pos] <- ifelse(aln_mat[, gappy_pos] == "-", "-", "N")
  ncol(aln_mat) == unique(width(aln))
  masked_aln <- DNAStringSet(apply(aln_mat, 1, paste0, collapse = ""))

  # Save masked alignment
  save_name <- gsub(".aln", str_glue(".masked_to_{n_unmasked}pos.aln"), aln_path)
  save_name <- gsub("/without_outgroup/", "/without_outgroup.masked/", save_name)
  writeXStringSet(masked_aln, save_name)
  rm(trimmed)
}

