rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(randomcoloR)

dat <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.gt1000nt.220723.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))

host_meta <- fread("data/metadata/parsed_host_metadata.csv") %>%
  select(-taxid, -host_species) %>%
  distinct() %>%
  mutate(is_dup = duplicated(host)) %>%
  filter(!is_dup) %>%
  select(-is_dup)

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

dat_filt <- dat %>%
  filter(family %in% to_keep) %>%
  filter(grepl("complete genome", genbank_title, ignore.case = T)) %>%
  filter(host %in% host_meta$host) %>%
  left_join(host_meta) %>%
  filter(is_vertebrate)

fam_counts <- dat_filt %>%
  filter(is_vertebrate) %>%
  group_by(family) %>%
  summarise(n_genomes = n(),
            n_hosts = n_distinct(host_order)) %>%
  arrange(desc(n_genomes)) %>%
  head(20)

pal <- distinctColorPalette(n_distinct(fam_counts$family))

fam_counts %>%
  dplyr::rename(viral_family = "family") %>%
  mutate(viral_family = factor(viral_family, unique(fam_counts$family))) %>%
  ggplot(aes(x = viral_family, 
             y = n_hosts, 
             size = log10(n_genomes),
             fill = viral_family)) +
  geom_point(pch = 21, color = "black") +
  scale_fill_manual(values = pal) + 
  guides(fill = F) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   face = "italic")) +
  labs(x = "Viral family", 
       y = "No. of host orders",
       size = "Log10(complete genomes)")


country_counts <- dat_filt %>%
  filter(country != "") %>%
  filter(!grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
                genbank_title, ignore.case = T)) %>%
  group_by(country) %>%
  summarise(n_genomes = n()) %>%
  arrange(desc(n_genomes)) %>%
  head(20)

country_by_fam <- dat_filt %>%
  filter(!grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
                genbank_title, ignore.case = T)) %>%
  filter(family %in% fam_counts$family) %>%
  filter(country %in% country_counts$country) %>%
  group_by(country, family) %>%
  summarise(n_genomes = n())

country_by_fam %>%
  mutate(country = factor(country, unique(country_counts$country))) %>%
  ggplot(aes(x = country, y = n_genomes, fill = family)) +
  geom_bar(stat = "identity", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = pal)
  