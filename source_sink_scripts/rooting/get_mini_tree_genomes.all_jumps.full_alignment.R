rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(Hmisc)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")

parsed <- meta %>%
  filter(host_genus != "")

count_df <- parsed %>%
  group_by(cluster) %>%
  summarise(n_genomes = n(),
            genus_range = n_distinct(host_genus),
            is_human = sum(host_genus == "Homo") > 0)

parsed_filt <- parsed %>% 
  left_join(count_df) %>%
  filter(genus_range > 1) %>%
  filter(n_genomes >= 10)

# Write clique genomes
fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

for(cluster_name in unique(parsed_filt$cluster)) {
  print(cluster_name)
  family_name <- str_split(cluster_name, "_")[[1]][1]
  
  cluster_accs <- deframe(meta %>%
                            filter(cluster == cluster_name) %>%
                            select(accession))
  
  fna_filt <- fna[names(fna) %in% cluster_accs]
  
  if (length(fna_filt) == length(cluster_accs)) {
    writeXStringSet(fna_filt, str_glue("data/genomes/source_sink_mini_trees/full_alignments/{cluster_name}.n{length(fna_filt)}.fna"))
  } else {
    print(str_glue("Eff!problem with {cluster_name} mini-tree"))
  }
}

print(str_glue("You should have {length(unique(parsed_filt$cluster))} fna files"))
