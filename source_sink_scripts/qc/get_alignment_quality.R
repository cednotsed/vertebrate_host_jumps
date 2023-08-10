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
require(ape)

meta <- fread(str_glue("results/clique_classification_out/final_cluster_metadata.220723.csv"))

all_genomes <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

aln_list <- list.files("data/alignments/source_sink_mini_trees/without_outgroup.masked/", 
                       "masked_to",
                       full.names = T)

clique_list <- str_split(aln_list, "/", simplify = T)
clique_list <- clique_list[, ncol(clique_list)]
clique_list <- str_split(clique_list, "\\.", simplify = T)[, 1]
clique_list

length_df <- meta %>%
  group_by(cluster) %>%
  summarise(median_length = median(genome_length))

morsels <- foreach(clique_name = clique_list) %do% {
  # clique_name = clique_list[1]
  
  # Get length of alignment
  aln_path <- aln_list[grepl(str_glue("{clique_name}\\."), aln_list)]
  aln <- read.dna(aln_path, format = "fasta")
  # aln_length <- ncol(aln)
  
  median_length <- deframe(length_df %>%
    filter(cluster == clique_name))
  
  # Get prop alignment kept
  pos_string <- str_split(aln_path, "/")[[1]][5]
  pos_string <- str_split(pos_string, "\\.")[[1]][3]
  n_unmasked <- as.numeric(gsub("masked_to_|pos", "", pos_string))
  
  prop_unmasked <- n_unmasked / median_length
  
  tibble(viral_clique = clique_name, prop_length_retained = prop_unmasked)
}

res <- bind_rows(morsels)

meta %>%
  filter(cluster == "Hepadnaviridae_1") %>%
  distinct(genus)
res %>% arrange(prop_length_retained)
res %>%
  ggplot() +
  geom_histogram(aes(x = prop_length_retained),
                 bins = 50) +
  labs(x = "Prop. of genome length retained",
       y = "No. of viral cliques") +
  theme_bw()

ggsave("results/source_sink_analysis/prop_alignment_unmasked.png", dpi = 300, width = 5, height = 5)

res %>%
  filter(prop_length_retained < 0.8)
# Save
res %>%
  dplyr::rename(clique_name = viral_clique) %>%
  filter(prop_length_retained > 0.8) %>%
  fwrite("results/qc_out/good_alignments.csv")
