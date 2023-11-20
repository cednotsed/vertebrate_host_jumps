rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(MSA2dist)
require(doParallel)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  filter(anc_name != "Root")

clique_list <- deframe(jump_df %>%
                         distinct(clique_name))

clique_list <- deframe(jump_df %>%
  distinct(clique_name))

# Do missing runs
# done <- list.files("results/dnds_out/all_jumps.temp_results/")
# done <- gsub(".csv", "", done)
# 
# clique_list <- clique_list[!(clique_list %in% done)]
# 
# print(length(clique_list))

cl <- makeCluster(12)
registerDoParallel(cl)

morsels <- foreach(clique = clique_list) %do% {
  # clique = clique_list[51]
  # clique = "Sedoreoviridae_19"
  jump_temp <- jump_df %>%
    filter(clique_name == clique)
  
  gene_list <- list.files(str_glue("results/dnds_out/mini_gene_alignments_sorted/{clique}"),
                          full.names = T)
  print(clique)
  
  crumbs <- foreach(i = seq(nrow(jump_temp)),
                    .packages = c("tidyverse", "Biostrings", "foreach", "data.table", "MSA2dist")) %dopar% {
    # i = 41
    row <- jump_temp[i, ]
    is_jump <- row$is_jump
    anc_name <- row$anc_name
    tip_name <- row$tip_name
    anc_state <- row$anc_state
    tip_state <- row$tip_state
    n_traverses <- row$n_traverses
    patristic_dist <- row$patristic_dist
    
    file_subset <- gene_list[grepl(tip_name, gene_list)]
    
    # start <- Sys.time()
    if(!identical(file_subset, character(0))) {
      concat_cds <- foreach(gene_file = file_subset, .combine = "cbind") %do% {
        # gene_file = file_subset[8]
        as.matrix(readDNAStringSet(gene_file, use.names = F))
      }
      
      # end <- Sys.time()
      # print(end - start)
      
      concat_seq <- DNAStringSet(apply(concat_cds, 1, paste0, collapse = ""))
      names(concat_seq) <- c("tip", "anc")
      dnds <- dnastring2kaks(concat_seq, model = "Li")
      ka <- as.numeric(dnds[, "ka"])
      ks <- as.numeric(dnds[, "ks"])
      
      temp <- tibble(clique_name = clique, 
                     is_jump = is_jump,
                     anc_name = anc_name,
                     tip_name = tip_name,
                     anc_state = anc_state,
                     tip_state = tip_state,
                     n_traverses = n_traverses,
                     patristic_dist = patristic_dist,
                     n_genes = length(file_subset),
                     cds_length = width(concat_seq)[1],
                     ka = ka,
                     ks = ks)
      
      return(temp)
    } else {
      print(str_glue("no gene annotations for {clique}:{tip_name}"))
      return(NULL)
    }
  }
  
  merged <- bind_rows(crumbs)
  fwrite(merged, str_glue("results/dnds_out/all_jumps.temp_results/{clique}.csv"))
  
  return(merged)
  
}

stopCluster(cl)

# Post-hoc merging
file_list <- list.files("results/dnds_out/all_jumps.temp_results/", full.names = T)
file_list

res <- foreach(file = file_list, .combine = "bind_rows") %do% {
  fread(file)
}

fwrite(res, "results/dnds_out/all_jumps.dnds.diff_hosts.genus_counts.csv")
