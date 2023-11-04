rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(taxizedb)
require(doParallel)

dat <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))

all_tax_match <- gsub(" sp\\.", "", unique(dat$host)) %>% 
  name2taxid(db = "ncbi", out_type = "summary") %>%
  mutate(rank = taxid2rank(id)) %>%
  dplyr::rename(host = name, taxid = id, rank = rank) %>%
  filter(rank %in% c("no rank", "species group", "subspecies", 
                     "species", "subgenus", "genus", 
                     "subfamily", "family", "superfamily", 
                     "tribe", "parvorder", "infraorder",
                     "order", "infraclass", "class",
                     "suborder", "subphylum", "phylum")) %>%
  distinct()

# unique(dat$host) %>% 
#   name2taxid(db = "ncbi", out_type = "summary") %>%
#   mutate(rank = taxid2rank(id)) %>%
#   distinct(rank) %>% View()
# cl <- makeCluster(12)
# registerDoParallel(cl)

# Get order and vertebrate status
morsels <- foreach(i = seq(nrow(all_tax_match))) %do% {
  row <- all_tax_match[i, ]
  taxid <- row$taxid
  taxonomy <- classification(taxid)[[1]]
  is_vert <- ifelse("Vertebrata" %in% taxonomy$name, T, F)
  # subp <- ifelse(identical(subp, character(0)), "Unknown", subp)
  
  species <- deframe(taxonomy %>%
                       filter(rank == "species") %>%
                       select(name))
  genus <- deframe(taxonomy %>%
                     filter(rank == "genus") %>%
                     select(name))
  family <- deframe(taxonomy %>%
                      filter(rank == "family") %>%
                      select(name))
  ord <- deframe(taxonomy %>%
                   filter(rank == "order") %>%
                   select(name))
  
  cls <- deframe(taxonomy %>%
                   filter(rank == "class") %>%
                   select(name))
  
  tibble(host = row$host, 
         taxid = taxid, 
         is_vertebrate = is_vert,
         species = ifelse(identical(species, character(0)), NA, species),
         genus = ifelse(identical(genus, character(0)), NA, genus),
         family = ifelse(identical(family, character(0)), NA, family),
         order = ifelse(identical(ord, character(0)), NA, ord),
         class = ifelse(identical(cls, character(0)), NA, cls))
}

host_meta <- bind_rows(morsels) %>%
  dplyr::rename(host_class = class,
         host_order = order,
         host_family = family,
         host_genus = genus,
         host_species = species) %>%
  mutate(is_dup = duplicated(host)) %>%
  filter(!is_dup) %>%
  select(-is_dup)

fwrite(host_meta, "data/metadata/parsed_host_metadata.csv")  


# # Fix metadata
# host_meta <- fread("data/metadata/parsed_host_metadata.csv")
# host_meta %>%
#   mutate(is_dup = duplicated(host)) %>%
#   filter(!is_dup) %>%
#   select(-is_dup) %>%
#   fwrite("data/metadata/parsed_host_metadata.csv")

