rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(data.table)
require(tidyverse)
require(Biostrings)

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")
fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

length(fna)

for(fam in unique(meta$family)) {
  # fam <- "Anelloviridae"
  meta_filt <- meta %>%
    filter(family == fam)
  
  temp_fna <- fna[names(fna) %in% meta_filt$accession]  
  writeXStringSet(temp_fna, str_glue("data/genomes/viral_family_subsets/{fam}.220723.n{length(temp_fna)}.fna"))
}

print("PLEASE RMBR TO DOS2UNIX!!!!!")

# Check difference
family_fnas <- list.files("data/genomes/viral_family_subsets/", full.names = T)
for (path in family_fnas) {
  old_path <- gsub("viral_family_subsets", "viral_family_subsets_OLD", path)
  new <- readDNAStringSet(path)
  old <- readDNAStringSet(old_path)
  is_match <- all(names(new) %in% names(old))
  is_count <- length(new) == length(old)
  print(str_glue("{path}: {is_match} {is_count}"))
}
