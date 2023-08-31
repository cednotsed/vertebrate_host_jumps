rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(igraph)
require(fossil)
require(aricode)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv") %>%
  mutate(species = tolower(species))

tax_meta <- fread("data/metadata/ICTV_Master_Species_List_2022_MSL38.v2.060823.csv") %>% 
  select(species = Species) %>%
  mutate(species = tolower(species))

done <- list.files("results/clique_classification_out/clustering_metrics/", "ICTV.csv")
done <- gsub(".AMI_ARI.ICTV.csv", "", done)
done

mat_list <- list.files("results/mash_out/viral_family_subsets/", "\\.tsv", full.names = T)
mat_list

to_do <- meta %>% 
  distinct(family) %>%
  filter(!(family %in% done))

# cl <- makeCluster(4)
# registerDoParallel(cl)

fam_morsels <- foreach(fam = to_do$family,
                       .packages = c("igraph", "foreach",
                                     "tidyverse", "aricode",
                                     "data.table")) %do% {
  mat_path <- mat_list[grepl(fam, mat_list)]
                                       
  mat <- fread(mat_path) %>%
    as_tibble() %>%
    column_to_rownames("#query")
  
  acc <- names(mat)
  
  mat <- as.matrix(mat)

  # # For subsampling
  # subsampled_index <- sample(seq(nrow(mat)), 500, replace = F)
  # dist_mat <- mat[subsampled_index, subsampled_index]
  # dist_mat[1:5, 1:5]
  # mat <- mat[1:5, 1:5]

  # Create network
  g <- graph_from_adjacency_matrix(mat, 
                                   weighted = T, 
                                   mode = "undirected")
  
  species_filt <- deframe(tibble(accession = acc) %>%
    left_join(meta) %>%
    group_by(species) %>%
    summarise(n = n()) %>% 
    inner_join(tax_meta) %>%
    select(species))
  
  morsels <- foreach(threshold = seq(0, 0.5, 0.05)) %do% {
    # threshold <- 0.15
    g_filt <- delete.edges(g, which(E(g)$weight > threshold))
    
    comm <- cluster_infomap(
      g_filt,
      modularity = F
    )
    
    # Match clusters with species metadata
    meta.match <- tibble(accession = comm$names, 
                         cluster = comm$membership) %>%
      left_join(meta) %>%
      filter(species %in% species_filt)
    
    # Calculate NMI
    ami_species <- AMI(meta.match$species, meta.match$cluster)
    ami_genus <- AMI(meta.match$genus, meta.match$cluster)
    ari_species <- ARI(meta.match$species, meta.match$cluster)
    ari_genus <- ARI(meta.match$genus, meta.match$cluster)
    
    # Results
    bind_rows(
      tibble(t = threshold, 
           type = "species",
           ami = ami_species,
           ari = ari_species,
           n_taxa = n_distinct(meta.match$species),
           n_clusters = n_distinct(meta.match$cluster)),
      tibble(t = threshold, 
             type = "genus",
             ami = ami_genus,
             ari = ari_genus,
             n_taxa = n_distinct(meta.match$genus),
             n_clusters = n_distinct(meta.match$cluster)))
  }
  
  plot_df <- bind_rows(morsels) %>%
    mutate(family = fam)
    
  fwrite(plot_df, 
         str_glue("results/clique_classification_out/clustering_metrics/{fam}.AMI_ARI.ICTV.csv"))
  
  return(plot_df)
}


