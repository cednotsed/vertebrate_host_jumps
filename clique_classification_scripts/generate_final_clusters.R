rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)
require(ggtree)
require(igraph)
require(Hmisc)

threshold <- 0.15
meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  select(-host_genus, -host_species)

set.seed(66)

to_do <- unique(meta$family)
mash_list <- list.files("results/mash_out/viral_family_subsets/", 
                        "\\.tsv",
                        full.names = T)

fam_morsels <- foreach(fam = to_do) %do% {
  print(fam)
  
  mash_path <- mash_list[grepl(fam, mash_list)]
  
  mat <- fread(mash_path) %>%
    as_tibble() %>%
    column_to_rownames("#query")
  
  acc <- names(mat)
  mat <- as.matrix(mat)
  
  # Create network
  g <- graph_from_adjacency_matrix(mat, 
                                   weighted = T, 
                                   mode = "undirected")
  
  g_filt <- delete.edges(g, which(E(g)$weight > threshold))
  
  comm <- cluster_infomap(
    g_filt,
    modularity = F
  )
  
  # Match clusters with species metadata
  meta.match <- tibble(accession = comm$names, 
                       cluster = comm$membership) %>%
    left_join(meta) %>%
    mutate(cluster = str_glue("{fam}_{cluster}"))

  return(meta.match)
}

# # Append columns
to_add <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  select(accession, molecule_type, is_segmented, is_circular, host)

final_meta <- bind_rows(fam_morsels) %>%
  separate(host, c("host_genus", "host_species"), " ") %>%
  mutate(host_genus = capitalize(host_genus)) %>%
  left_join(to_add)
  # filter(!(grepl("dae|nae", host_genus) & !grepl("naeus|naem|naes|naeo", host_genus)))

final_meta %>%
  group_by(family) %>%
  summarise(n = n_distinct(cluster)) %>%
  arrange(desc(n))

fwrite(final_meta, "results/clique_classification_out/final_cluster_metadata.220723.new.csv")

final_meta %>%
  distinct(accession) %>%
  nrow()

final_meta %>%
  nrow()

# # Parse host columns and add host metadata
# old_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.OLD.csv")
# host_meta <- fread("data/metadata/parsed_host_metadata.csv")
# 
# # Parse host with more than two words
# word_count <- old_meta %>%
#   mutate(original_host = host) %>%
#   mutate(host = gsub("  ", " ", host)) %>%
#   mutate(n_words = str_count(host, "\\w+"))
# 
# short_names <- word_count %>%
#   filter(n_words <= 2) %>%
#   mutate(host = gsub(" sp\\.| subsp\\.| subgen\\.", "", host))
# 
# long_names <- word_count %>%
#   filter(n_words > 2) %>%
#   separate(host, c("host1", "host2"), " ", remove = F) %>%
#   mutate(host = str_glue("{host1} {host2}")) %>%
#   mutate(host = gsub(" sp\\.| subsp\\.| subgen\\.", "", host))
# 
# combined <- bind_rows(long_names, short_names) %>%
#   mutate(host = ifelse(!(host %in% host_meta$host),
#                        word(host, 1),
#                        host)) %>%
#   left_join(host_meta) %>%
#   select(-host1, -host2, -n_words)
# 
# fwrite(combined, "results/clique_classification_out/final_cluster_metadata.220723.csv",
#        eol = "\n")
