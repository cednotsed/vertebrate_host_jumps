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

# Get host list column
host_morsels <- foreach(clique_name = unique(parsed_filt$cluster)) %do% {
  print(clique_name)
  
  host_list <- deframe(parsed_filt %>%
                         filter(cluster == clique_name) %>%
                         distinct(host_genus) %>%
                         filter(!(host_genus %in% c("", "Homo"))))
  
  tibble(cluster = clique_name, 
         hosts = paste0(host_list, collapse = ";"))
}

# Previous analysis
# previous <- fread("results/source_sink_analysis/source_sink_results.curated.csv") %>%
#   mutate(to_keep = ifelse(to_keep == "", NA, to_keep))

to_write <- bind_rows(host_morsels) %>%
  left_join(parsed_filt %>% distinct(cluster, n_genomes))

fwrite(to_write, "results/source_sink_analysis/source_sink_results.csv")

# Write clique genomes
fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

mash_list <- list.files("results/mash_out/viral_family_subsets/", 
                        "\\.tsv",
                        full.names = T)

for(cluster_name in unique(parsed_filt$cluster)) {
  print(cluster_name)
  family_name <- str_split(cluster_name, "_")[[1]][1]
  
  cluster_accs <- deframe(meta %>%
                            filter(cluster == cluster_name) %>%
                            select(accession))
  
  # Get +-20 closest isolates
  buffer <- 20
  mash_path <- mash_list[grepl(family_name, mash_list)]
  
  family_dist <- fread(mash_path) %>%
    rename(`#query` = "query")
  
  to_compute <- family_dist$query[!(family_dist$query %in% cluster_accs)]
  
  dist_filt <- family_dist %>%
    filter(query %in% cluster_accs) %>%
    select(all_of(c("query", to_compute))) %>%
    column_to_rownames("query")
  
  col_min <- apply(dist_filt, 2, min)
  col_min_df <- tibble(acc = names(col_min),
                       min_dist = col_min) %>%
    arrange(min_dist) %>%
    head(buffer)

  to_keep <- c(col_min_df$acc, cluster_accs)
  
  fna_filt <- fna[names(fna) %in% to_keep]
  
  if (length(fna_filt) == length(to_keep)) {
    writeXStringSet(fna_filt, str_glue("data/genomes/source_sink_mini_trees/with_buffer_outgroup/{cluster_name}.n{length(fna_filt)}.fna"))
  } else {
    print(str_glue("Eff!problem with {cluster_name} mini-tree"))
  }
}
