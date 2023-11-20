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

cl <- makeCluster(16)
registerDoParallel(cl)

morsels <- foreach(clique = clique_list) %do% {
  # clique = clique_list[8]
  gene_list <- list.files(str_glue("results/dnds_out/mini_gene_alignments_sorted/{clique}"),
                          full.names = F)
  print(clique)
  
  crumbs <- foreach(gene_file = gene_list,
                    .packages = c("tidyverse", "Biostrings", "foreach", "data.table", "MSA2dist")) %dopar% {
    # gene_file = gene_list[73]
    # print(gene_file)
    .GlobalEnv$clique <- clique                  
    file_string <- str_split(gene_file, "-")[[1]]
    anc_name <- file_string[2]
    tip_name <- file_string[3]
    gene_annot <- gsub("\\.fna", "", file_string[4])
    
    gene_aln <- readDNAStringSet(str_glue("results/dnds_out/mini_gene_alignments_sorted/{clique}/{gene_file}"))
    
    dnds <- dnastring2kaks(gene_aln, model = "Li")
    ka <- as.numeric(dnds[, "ka"])
    ks <- as.numeric(dnds[, "ks"])
    
    temp <- tibble(anc_name = anc_name,
                   tip_name = tip_name,
                   gene_length = width(gene_aln)[1],
                   gene_annot = gene_annot,
                   ka = ka,
                   ks = ks)
      
    return(temp)
  }
  
  merged <- bind_rows(crumbs) %>%
    inner_join(jump_df)
  
  fwrite(merged, str_glue("results/dnds_out/all_jumps.by_gene.temp_results/{clique}.csv"))
  
  return(merged)
  
}

stopCluster(cl)

res <- bind_rows(morsels)
fwrite(res, "results/dnds_out/all_jumps.by_gene.diff_hosts.genus_counts.csv")
