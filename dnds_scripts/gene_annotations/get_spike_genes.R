rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)
require(Biostrings)

genome_type <- fread("data/metadata/genome_type_metadata.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")

# Get host counts
host_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  filter(host_genus != "") %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_genus))

genome_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(n_genomes = n_distinct(accession))

genome_length <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(median_length = median(genome_length))

dnds_list <- list.files("results/dnds_out/all_jumps.by_gene.temp_results/", full.names = T)
dnds_df <- foreach(file_name = dnds_list, .combine = "bind_rows") %do% {
  if(file.size(file_name) > 246) {
    fread(file_name)
  }
}

dnds_filt <-  dnds_df

clique_counts <- dnds_filt %>%
  left_join(host_counts) %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- dnds_filt %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

cov_df <- merged_df %>%
  filter(grepl("Coronaviridae_12", clique_name))

cov_split <- cov_df %>% 
  # distinct(clique_name)
  mutate(gene_type = case_when(grepl("ORF1|1a|1b|polyprotein|non_struct|ns", gene_annot, ignore.case = T) ~ "ORF1ab",
                               grepl("internal|I_protein|i_protein|protein_i", gene_annot, ignore.case = T) ~ "I protein",
                               grepl("matrix|membrane|protein_M|M_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("M") ~ "Membrane",
                               grepl("N2|N_protein|nucelo|nucleo|nucloe|envelope|E_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("E") ~ "Envelope",
                               grepl("nucelo|nucleo|nucloe", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("N") ~ "Nucleocapsid",
                               # grepl("internal|I_protein|matrix|N_protein|nucelo|nucleo|nucloe|membrane|protein_M|M_protein|envelope|E_protein", gene_annot, ignore.case = T) |
                               #   gene_annot %in% c("M", "E", "N")~ "Structural",
                               grepl("S_protein|spike|glyco|surface|hema|hemma|esterase", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("S", "HE") ~ "Entry",
                               grepl("ORF3a|3a|3b|3c|protein_3", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("orf3a") ~ "ORF3",
                               grepl("4a|4b|4c", gene_annot, ignore.case = T) ~ "ORF4",
                               grepl("5a|5b|5c", gene_annot, ignore.case = T) ~ "ORF5",
                               grepl("sars6|6a|6b|protein_6|6_protein", gene_annot, ignore.case = T) ~ "ORF6",
                               grepl("7a|7b|protein_7", gene_annot, ignore.case = T) ~ "ORF7",
                               grepl("8a|8b|protein_8", gene_annot, ignore.case = T) ~ "ORF8",
                               grepl("9a|9b|protein_9", gene_annot, ignore.case = T) ~ "ORF9",
                               grepl("32_kda|32kda", gene_annot, ignore.case = T) ~ "32kDa protein")) 
  # filter(!is.na(gene_type)) %>%
  # filter(gene_type != "ORF9")

spike <- cov_df %>%
  # filter(gene_type == "Entry") %>%
  filter((grepl("S_protein|spike|surface", gene_annot, ignore.case = T) |
           gene_annot %in% c("S")) &
           !grepl("hema", gene_annot))


# gene_list <- list.files("results/dnds_out/mini_gene_alignments_sorted/Coronaviridae_12", full.names = T)
# ref_path <- "results/dnds_out/family_gene_alignments/SC2-spike.fna"

# gene_list <- list.files("results/dnds_out/mini_gene_alignments_sorted/Coronaviridae_18", full.names = T)
# ref_path <- "results/dnds_out/family_gene_alignments/PDCoV-spike.fna"

gene_list <- list.files("results/dnds_out/mini_gene_alignments_sorted/Coronaviridae_26", full.names = T)
ref_path <- "results/dnds_out/family_gene_alignments/MERS-spike.fna"
# 
# gene_list <- list.files("results/dnds_out/mini_gene_alignments_sorted/Coronaviridae_6", full.names = T)
# ref_path <- "results/dnds_out/family_gene_alignments/IBV-spike.fna"

spike_paths <- gene_list[grepl("spike|surface", gene_list)]

spike_fna <- foreach(gene_path = spike_paths, .combine = "c") %do% {
  # print(gene_path)
  temp <- readDNAStringSet(gene_path)
  return(Biostrings::chartr("-", "N", temp))
}

ref <- readDNAStringSet(ref_path)
spike_fna <- c(ref, spike_fna)

save_path <- gsub("spike.fna", "spike.all_genomes.fna", ref_path)
writeXStringSet(spike_fna, save_path)
# clique_list <- deframe(spike %>%
#                          distinct(clique_name))
names(spike_fna)
# cl <- makeCluster(16)
# registerDoParallel(cl)
# 
# spike_fasta <- foreach(clique = clique_list, .combine = "c") %do% {
#   # clique = clique_list[1]
#   gene_list <- list.files(str_glue("results/dnds_out/mini_gene_alignments_sorted/{clique}"),
#                           full.names = T)
#   
#   gene_paths <- gene_list[grepl("S_protein|S1_glyco|spike|surface_glycoprotein|S\\.fna", gene_list, ignore.case = T)]
#   
#   crumbs <- foreach(gene_path = gene_paths, .combine = "c", .packages = c("Biostrings")) %dopar% {
#     print(gene_path)
#     temp <- readDNAStringSet(gene_path)
#     return(Biostrings::chartr("-", "N", temp))
#   }
#   
#   return(crumbs)
# }
# 
# writeXStringSet(spike_fasta, "results/dnds_out/family_gene_alignments/CoV_spike.fna")
# stopCluster(cl)
# 
# n_distinct(spike$tip_name)
#   
