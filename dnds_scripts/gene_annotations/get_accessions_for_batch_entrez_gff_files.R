rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(Biostrings)

# fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

jump_df <- fread("results/source_sink_analysis/putative_host_jumps.csv")

accs <- jump_df %>% 
  distinct(tip_name)

# Download new files only
done <- fread("results/dnds_out/ncbi_annotation/gene_annotations.290723.n8404.parsed.gff3")$V1

accs %>%
  filter(!(tip_name %in% done)) %>%
  fwrite("results/dnds_out/ncbi_annotation/to_download.accessions.2.txt",
         col.names = F,
         eol = "\n")


# # Get all accessions
# meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
# meta %>%
#   select(accession) %>%
#   fwrite(str_glue("results/dnds_out/ncbi_annotation/all_gff.accessions.txt"),
#          col.names = F,
#          eol = "\n")
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
