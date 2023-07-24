rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)
require(ggtree)
require(igraph)

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  mutate(species = tolower(species))

mash_list <- list.files("results/mash_out/viral_family_subsets/", 
                        pattern = "\\.tsv",
                        full.names = T)

done <- list.files("results/clique_classification_out/monophyly/", ".csv")
done <- gsub(".monophyly.csv", "", done)
done

mat_list <- list.files("results/mash_out/viral_family_subsets/", "\\.tsv", full.names = T)

to_do <- meta %>% 
  distinct(family) %>%
  filter(!(family %in% done))

fam_morsels <- foreach(fam = unique(to_do$family)) %do% {
  print(fam)
  mat_path <- mash_list[grepl(fam, mash_list)]
  mat <- fread(mat_path) %>%
    as_tibble() %>%
    column_to_rownames("#query")
  
  acc <- names(mat)
  mat <- as.matrix(mat)
  
  # Plot NJ tree
  print("Computing NJ...")
  cm <- data.matrix(mat)
  tree <- bionj(cm)
  
  write.tree(tree, str_glue("data/trees/viral_family_subsets/{fam}.220723.mash_dist.nwk"))
  
  # Iterate through thresholds
  print("Estimating monophyly...")
  mono_morsels <- foreach(threshold = seq(0, 0.5, 0.05)) %do% {
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
      left_join(meta)
    
    clus_morsels <- foreach(cluster_no = unique(meta.match$cluster),
            .combine = "c") %do% {
      is.monophyletic(tree, 
                      meta.match[meta.match$cluster == cluster_no, ]$accession)          
    }
    
    prop_mono <- sum(clus_morsels) / length(clus_morsels)
    
    tibble(family = fam, t = threshold, prop_mono = prop_mono)
  }
  
  temp_df <- bind_rows(mono_morsels)
  
  fwrite(temp_df, str_glue("results/clique_classification_out/monophyly/{fam}.monophyly.csv"))
  
  return(temp_df)
}
