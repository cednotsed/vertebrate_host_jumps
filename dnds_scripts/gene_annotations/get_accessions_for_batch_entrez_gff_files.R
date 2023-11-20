rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(Biostrings)

# Get all accessions
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv")

meta %>%
  select(accession) %>%
  fwrite(str_glue("results/dnds_out/ncbi_annotation/all_gff.accessions.txt"),
         col.names = F,
         eol = "\n")

# Download new files only
done_df <- fread("results/dnds_out/ncbi_annotation/all_gene_annotations.061123.parsed.gff3")

meta %>% 
  distinct(accession) %>%
  filter(!(accession %in% unique(done_df$V1))) %>%
  fwrite("results/dnds_out/ncbi_annotation/to_download.accessions.2.txt",
         col.names = F,
         eol = "\n")

to_do <- accs %>%
  filter(!(tip_name %in% done))
  # fwrite("results/dnds_out/ncbi_annotation/to_download.accessions.2.txt",
  #        col.names = F,
  #        eol = "\n")

# fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")
# good_alns <- fread("results/qc_out/good_alignments.csv")
# jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
#   filter(clique_name %in% good_alns$clique_name)
# 
# accs <- jump_df %>% 
#   distinct(tip_name)

# 
# test <- jump_df %>%
#   separate(clique_name, c("family"), "_") %>%
#   filter(family %in% c("Spinareoviridae", "Sedoreoviridae"))
# jump_df

# 
# chunk_list <- split(meta$accession,           
#                     ceiling(seq_along(meta$accession) / 9999))
# 
# for (i in seq(length(chunk_list))) {
#   print(i)
#   chunk <- chunk_list[[i]]
#   print(length(chunk))
#   fwrite(tibble(chunk), 
#          str_glue("results/dnds_out/ncbi_annotation/gff_chunks.accessions.{i}.txt"),
#          col.names = F,
#          eol = "\n")
# }
# 
# 
# gff_list <- list.files("results/dnds_out/ncbi_annotation/", "sequence.chunk", full.names = T)
# gff_list <- gff_list[grepl("parsed", gff_list)]
# morsels <- foreach(gff = gff_list) %do% {
#   fread(gff)
# }
# 
# bind_rows(morsels) %>%
#   filter(V3 == "CDS") %>%
#   fwrite("results/dnds_out/ncbi_annotation/all_gene_annotations.290723.n53xxx.parsed.gff3",
#          eol = "\n",
#          col.names = F)
# meta$accession[!(meta$accession %in% done)]
# meta %>%
#   filter(!(accession %in% done)) %>%
#   select(accession) %>%
#   fwrite("results/dnds_out/ncbi_annotation/missing.txt",
#          col.names = F,
#          eol = "\n")
# 
# bind_rows(morsels) %>%
#   # filter(V3 == "CDS") %>%
#   group_by(V1) %>%
#   summarise(n = n()) %>%
#   filter(n == 2)
# fna_filt <- fna[accs$tip_name]

# writeXStringSet(fna_filt, "results/dnds_out/vapid_annotation/for_vapid/references.fna")


# dummy_meta <- accs %>%
#   dplyr::rename(strain = tip_name) %>%
#   mutate("collection-date" = "Unknown",
#          country = "Unknown",
#          coverage = "Unknown")
# 
# fwrite(dummy_meta, "results/dnds_out/vapid_annotation/for_vapid/dummy_meta.csv")


# host_meta %>%
#   filter(host == "Mustela sp.")
