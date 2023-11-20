rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(Biostrings)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
ref_master <- fread("data/metadata/dnds_analysis/ref_master.csv")
anc_seq_list <- list.files("results/ancestral_reconstruction_out/ancestral_sequences/", 
                           full.names = T)
tip_seq_list <- list.files("data/alignments/source_sink_mini_trees/without_outgroup.masked/",
                           full.names = T)
jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(anc_name != "Root")

for(i in seq(nrow(ref_master))) {
  clique <- ref_master[i, ]$clique_name
  ref_name <- ref_master[i, ]$ref_name
  virus <- ref_master[i, ]$virus
  start <- ref_master[i, ]$start_pos
  end <- ref_master[i, ]$end_pos
  
  jump_filt <- jump_df %>%
    filter(clique_name == clique)
  
  # Get alignment files
  anc_fna <- readDNAStringSet(anc_seq_list[grepl(clique, anc_seq_list)])
  tip_fna <- readDNAStringSet(tip_seq_list[grepl(clique, tip_seq_list)])

  aln <- foreach(j = seq(nrow(jump_filt)), .combine = "c") %do% {
    row <- jump_filt[j, ]
    
    # Get interleaved tip and anc seqs
    anc_seq <- anc_fna[row$anc_name]
    tip_seq <- tip_fna[row$tip_name]
    
    # rename anc_seq to avoid dupl. names
    names(anc_seq) <- str_glue("{row$anc_name}|{row$tip_name}")
    return(c(tip_seq, anc_seq))
  }
  
  # Trim to reference
  aln_mat <- as.matrix(aln)
  to_keep <- aln_mat[ref_name, ] != "-"
  aln_filt <- aln_mat[, to_keep]
  
  # Extract spike
  spike <- aln_filt[ , start:end]
  
  # Convert back to DNAss
  spike_fna <- DNAStringSet(apply(spike, 1, paste0, collapse = ""))
  
  # Replace gaps with Ns
  spike_parsed <- Biostrings::chartr("-", "N", spike_fna)
  
  writeXStringSet(spike_parsed,
                  str_glue("results/dnds_out/family_gene_alignments_V2/{virus}_spike.aln"))
}

jump_df %>%
  filter(clique_name == "Coronaviridae_6") %>%
  filter(is_jump)
