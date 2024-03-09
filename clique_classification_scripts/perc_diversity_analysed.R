rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)
require(Hmisc)
require(maps)
require(randomcoloR)

parsed_meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.parsed_host.csv")
# genome_meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")
ictv_meta <- fread("data/metadata/ICTV_Master_Species_List_2022_MSL38.v2.060823.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))
clique_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")

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

plot_df <- parsed_meta %>%
  filter(family %in% to_keep) %>%
  mutate(family = ifelse(family %in% clique_meta$family, family, "Others")) %>%
  filter(species != "") %>%
  group_by(family) %>%
  distinct(species) %>%
  mutate(is_analysed = species %in% clique_meta$species) %>%
  group_by(family) %>%
  summarise(n_analysed = sum(is_analysed), 
            n_total = n(),
            prop_analysed = sum(is_analysed) / n()) %>%
  arrange(desc(n_analysed))
  
plot_df %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family)) +
  geom_bar(aes(y = n_total),
               stat = "identity") +
  geom_bar(aes(y = n_analysed),
           stat = "identity",
           fill = "olivedrab") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.3,
                                   face = "italic")) +
  labs(x = "Viral family", y = "Viral species (NCBI Virus)") 

ggsave("results/clique_classification_out/perc_diversity_analysed.pdf", 
       dpi = 600, 
       width = 6, 
       height = 3.5)
plot_df %>%
  summarise(all_analysed = sum(n_analysed),
            all_total = sum(n_total)) %>%
  mutate(all_analysed / all_total)

length(to_keep)
# Proportion of viruses vertebrates
parsed_meta %>%
  group_by(is_vertebrate) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(prop = n / sum(n))

# Proportion of viruses by type
type_counts <- parsed_meta %>%
  filter(family != "Pleolipoviridae") %>% 
  mutate(molecule_type = ifelse(molecule_type %in% c("DNA", "RNA", "", 
                                                     "unknown"),
                                "unknown",
                                molecule_type)) %>%
  # filter(!grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2",
  #               genbank_title,
  #               ignore.case = T)) %>%
  # filter(molecule_type != "dsDNA; ssDNA") %>%
  mutate(molecule_type = case_when(grepl("dsDNA", molecule_type) ~ "dsDNA",
                                   grepl("ssDNA", molecule_type) ~ "ssDNA",
                                   grepl("dsRNA", molecule_type) ~ "dsRNA",
                                   grepl("ssRNA", molecule_type) ~ "ssRNA",
                                   molecule_type == "unknown" ~ "unknown")) %>%
  group_by(molecule_type) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n, na.rm = T) * 100)

type_counts %>%
  arrange(desc(prop))

# ictv_meta %>%
#   dplyr::rename(molecule_type = genome_composition) %>%
#   filter(molecule_type != "dsDNA; ssDNA") %>%
#   mutate(molecule_type = case_when(grepl("dsDNA", molecule_type) ~ "dsDNA",
#                                    grepl("ssDNA", molecule_type) ~ "ssDNA",
#                                    grepl("dsRNA", molecule_type) ~ "dsRNA",
#                                    grepl("ssRNA", molecule_type) ~ "ssRNA")) %>%
#   group_by(molecule_type) %>%
#   summarise(n_species = n_distinct(species)) %>%
#   left_join(type_counts) %>%
#   mutate(ratio = n / n_species) %>%
#   arrange(desc(ratio))


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