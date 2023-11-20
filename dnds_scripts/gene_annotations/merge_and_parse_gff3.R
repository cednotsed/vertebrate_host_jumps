rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(Biostrings)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv")

gff_list <- list.files("results/dnds_out/ncbi_annotation/gff_chunks", "sequence.chunk", full.names = T)
gff_list <- gff_list[grepl("parsed", gff_list)]

morsels <- foreach(gff = gff_list) %do% {
  temp <- fread(gff)
  
  print(n_distinct(temp$V1))
  return(temp)
}

done <- bind_rows(morsels)

# Get only CDS
done_filt <- done %>%
  filter(V3 == "CDS")

n_annots <- n_distinct(done_filt$V1)

done_filt %>%  
  fwrite(str_glue("results/dnds_out/ncbi_annotation/all_gene_annotations.290723.n{n_annots}.parsed.gff3"),
         eol = "\n",
         col.names = F)

# Check if all done
meta$accession[!(meta$accession %in% unique(done$V1))]

# If not get missing accessions
meta %>%
  filter(!(accession %in% unique(done$V1))) %>%
  select(accession) %>%
  fwrite("results/dnds_out/ncbi_annotation/missing.txt",
         col.names = F,
         eol = "\n")
