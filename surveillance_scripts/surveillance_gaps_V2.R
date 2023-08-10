rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.gt1000nt.220723.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(host != "") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv"))
genome_meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")

to_keep <- c("Poxviridae", "Herpesviridae", "Picornaviridae",
             "Peribunyaviridae", "Birnaviridae", "Circoviridae",
             "Rhabdoviridae", "Hantaviridae", "Polyomaviridae",
             "Paramyxoviridae", "Bornaviridae", "Filoviridae",
             "Nairoviridae", "Spinareoviridae", "Sedoreoviridae",
             "Phenuiviridae", "Arenaviridae", "Parvoviridae",
             "Redondoviridae", "Adenoviridae", "Togaviridae",
             "Orthomyxoviridae", "Arteriviridae", "Coronaviridae",
             "Tobaniviridae", "Nanhypoviridae", "Olifoviridae",
             "Gresnaviridae", "Nanghoshaviridae", "Cremegaviridae",
             "Caliciviridae", "Asfarviridae", "Alloherpesviridae",
             "Astroviridae", "Flaviviridae", "Retroviridae",
             "Nyamiviridae", "Metaviridae", "Hepadnaviridae",
             "Genomoviridae", "Smacoviridae", "Picobirnaviridae",
             "Reoviridae", "Anelloviridae", "Papillomaviridae",
             "Adomaviridae", "Hepeviridae", "Kolmioviridae",
             "Pneumoviridae", "Amnoonviridae", "Nodaviridae",
             "Matonaviridae", "Pseudoviridae")

plusRNA <- c("Picornaviridae", "Togaviridae", "Arteriviridae",
             "Coronaviridae", "Tobaniviridae", "Nanhypoviridae",
             "Olifoviridae", "Gresnaviridae", "Nanghoshaviridae",
             "Cremegaviridae", "Caliciviridae", "Astroviridae",
             "Flaviviridae", "Hepeviridae", 
             "Nodaviridae", "Matonaviridae")
minusRNA <- c("Peribunyaviridae", "Rhabdoviridae", "Hantaviridae",
              "Paramyxoviridae", "Bornaviridae", "Filoviridae",
              "Nairoviridae", "Phenuiviridae", "Orthomyxoviridae",
              "Nyamiviridae", "Kolmioviridae", "Pneumoviridae",
              "Amnoonviridae")
ambiRNA <- c("Arenaviridae")
dsRNA <- c("Birnaviridae", "Spinareoviridae", "Sedoreoviridae",
           "Picobirnaviridae", "Reoviridae")
ssDNA <- c("Circoviridae", "Parvoviridae", "Redondoviridae",
           "Genomoviridae", "Smacoviridae", "Anelloviridae")
dsDNA <- c("Poxviridae", "Herpesviridae", "Polyomaviridae",
           "Adenoviridae", "Asfarviridae", "Alloherpesviridae",
           "Papillomaviridae", "Adomaviridae")
ssRNA_RT <- c("Retroviridae", "Pseudoviridae", "Metaviridae")
dsDNA_RT <- c("Hepadnaviridae")

meta_filt <- meta %>%
  filter(!is.na(host_order),
         family %in% to_keep) %>%
  filter(host_order != "") %>%
  filter(host_family != "Hominidae") %>%
  mutate(genome_type = case_when(family %in% plusRNA ~ "+ssRNA",
                                 family %in% minusRNA ~ "-ssRNA",
                                 family %in% ambiRNA ~ "+/-ssRNA",
                                 family %in% dsRNA ~ "dsRNA",
                                 family %in% ssDNA ~ "ssDNA",
                                 family %in% dsDNA ~ "dsDNA",
                                 family %in% minusRNA ~ "-ssRNA",
                                 family %in% ssRNA_RT ~ "ssRNA (RT)",
                                 family %in% dsDNA_RT ~ "dsDNA (RT)"))

meta_filt %>%
  # filter(is_vertebrate) %>%
  group_by(host_genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

genome_meta %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv")) %>%
  distinct(species)
filter(is_vertebrate) %>%
  # group_by(family) %>%
  summarise(n = n_distinct(family)) %>% View()
# group_by(host_order) %>%
# filter(host_order != "") %>%
# filter(is_vertebrate) %>%
summarise(n = n()) %>%
  arrange(desc(n)) %>% View()

