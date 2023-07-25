rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(ggtree)
require(treeio)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)

meta <- fread(str_glue("results/clique_classification_out/final_cluster_metadata.220723.csv"))

all_genomes <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

aln_list <- list.files("data/alignments/source_sink_mini_trees/with_buffer_outgroup/", 
                       "trim_to",
                       full.names = T)

clique_list <- str_split(aln_list, "/", simplify = T)
clique_list <- clique_list[, ncol(clique_list)]
clique_list <- str_split(clique_list, "\\.", simplify = T)[, 1]
clique_list

morsels <- foreach(clique_name = clique_list) %do% {
  aln_path <- aln_list[grepl(str_glue("{clique_name}\\."), aln_list)]
  aln <- read.dna(aln_path, format = "fasta")
  aln_length <- ncol(aln)
  
  # Get average genome length for clique
  clique_acc <- deframe(meta %>%
    filter(cluster == clique_name) %>%
    select(accession))
  
  fna_filt <- all_genomes[clique_acc]
  avg_genome_length <- mean(width(fna_filt))
  
  tibble(viral_clique = clique_name, prop_length_retained = aln_length / avg_genome_length)
  
}

res <- bind_rows(morsels)

res %>%
  ggplot(aes(x = prop_length_retained)) +
  geom_histogram() +
  labs(x = "Prop. of genome length retained",
       y = "No. of viral cliques") +
  theme_bw()

ggsave("results/source_sink_analysis/prop_aligned.png", dpi = 300, width = 5, height = 5)
