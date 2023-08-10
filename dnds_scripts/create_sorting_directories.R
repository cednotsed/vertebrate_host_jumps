rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name)

clique_list <- jump_df %>%
  distinct(clique_name)

for(clique in clique_list$clique_name) {
  dir.create(str_glue("results/dnds_out/mini_gene_alignments_sorted/{clique}"))
}

  # fwrite("results/dnds_out/clique_list.txt",
  #        col.names = F,
  #        eol = "\n")


