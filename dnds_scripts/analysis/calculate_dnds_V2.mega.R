rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(MSA2dist)
require(doParallel)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name)

clique_list <- deframe(jump_df %>%
  distinct(clique_name))

jump_list <- jump_df %>%
  distinct(clique_name, is_jump) %>%
  filter(grepl("^A|^B|^C|^F|^G|^H|^O", clique_name)) %>%
  arrange(clique_name)

done <- list.files("results/dnds_out/all_jumps.mega.temp_results/")

cl <- makeCluster(16)
registerDoParallel(cl)

morsels <- foreach(i = seq(nrow(jump_list))) %do% {
  # i = 42
  clique <- jump_list[i, ]$clique_name
  is_jump <- jump_list[i, ]$is_jump
  print(clique)
  gene_list <- list.files(str_glue("results/dnds_out/mini_gene_alignments_sorted/{clique}"),
                          full.names = T)
  
  tip_names <- jump_df %>%
    filter(clique_name == clique,
           is_jump == is_jump)
  
  tip_regex <- paste0(tip_names$tip_name, collapse = "|")
  file_subset <- gene_list[grepl(tip_regex, gene_list)]
  
  start <- Sys.time()
  concat_cds <- foreach(gene_file = file_subset, 
                        .combine = "cbind",
                        .packages = c("tidyverse", "Biostrings", "MSA2dist")) %dopar% {
    as.matrix(readDNAStringSet(gene_file, use.names = F))
  }
  
  end <- Sys.time()
  print(end - start)
  
  concat_seq <- DNAStringSet(apply(concat_cds, 1, paste0, collapse = ""))
  names(concat_seq) <- c("tip", "anc")
  dnds <- dnastring2kaks(concat_seq, model = "Li")
  ka <- as.numeric(dnds[, "ka"])
  ks <- as.numeric(dnds[, "ks"])
  
  temp <- tibble(clique_name = clique, 
                 is_jump = is_jump,
                 n_genes = length(file_subset),
                 cds_length = width(concat_seq)[1],
                 ka = ka,
                 ks = ks)
  fwrite(temp, str_glue("results/dnds_out/all_jumps.mega.temp_results/{clique}_{i}.csv"))
  
  return(temp)  
}

stopCluster(cl)

bind_rows(morsels)
