rm(list = ls())
setwd("c:/git_repos/viral_sharing/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)

tax_meta <- fread("data/VIRION.v0.2.1_240922/TaxonomyVirus.csv") %>%
  filter(ICTVRatified) %>%
  select(species = Virus)

meta <- fread("data/metadata/all_viruses.140423.filt.QCed.csv") %>%
  mutate(species = tolower(species)) %>%
  inner_join(tax_meta)

to_keep <- deframe(fread("data/metadata/all_viruses.140423.filt.QCed.csv") %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  filter(n > 500) %>%
  distinct(family))
  
family_list <- unique(meta$family)[unique(meta$family) %in% to_keep]

morsels <- foreach(family_name = family_list) %do% {
  # family_name <- family_list[2]
  dist_df <- fread(str_glue("results/phylogenetic_out/mash_out/viral_family_subsets/{family_name}.140423.filt.QCed.formatted.tsv")) %>%
    column_to_rownames("#query")
    
  sp_list <- deframe(meta %>%
    filter(family == family_name) %>%
    distinct(species))
  
  crumbs <- foreach(sp = sp_list) %do% {
    # sp <- sp_list[1]
    accs <- deframe(meta %>%
      filter(family == family_name,
             species == sp) %>%
      distinct(accession))
    
    median_dist <- median(as.matrix(dist_df[accs , accs]))
    
    tibble(species = sp, median_dist = median_dist)
  }
  
  bind_rows(crumbs) %>%
    mutate(family = family_name)
}

bind_rows(morsels) %>%
  # filter(median_dist > 0.15) %>% View()
  # distinct(family)
  ggplot(aes(x = median_dist)) +
  geom_histogram() +
  geom_vline(xintercept = 0.15) +
  facet_grid(rows = vars(family))
