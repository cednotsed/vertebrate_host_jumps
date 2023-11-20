rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)
require(Hmisc)
require(maps)
require(randomcoloR)

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

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv")) %>%
  select(accession, host, host_species,
         host_genus, host_family, genbank_title, 
         country, family, species, 
         collection_date) %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, "")) %>%
  separate(collection_date, c("Y", "M", "D"), "-", remove = F) %>%
  mutate(Y = as.numeric(Y),
         M = as.numeric(M),
         D = as.numeric(D))


# Overall prop missing
meta %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, "")) %>%
  filter(host_species != "Homo sapiens",
         host_genus != "Homo") %>%
  mutate(host_genus = ifelse(host_genus == "",
                               "Missing", host_genus)) %>%
  summarise(prop_missing = sum(host_genus == "Missing") / n())

# Overall date missing
meta %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, "")) %>%
  filter(host_species != "Homo sapiens",
         host_genus != "Homo") %>%
  summarise(perc_Y_missing = sum(is.na(Y)) / n() * 100,
            perc_M_missing = sum(is.na(M)) / n() * 100,
            perc_D_missing = sum(is.na(D)) / n() * 100,
            n_genomes = n())

plot_df <- meta %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, "")) %>%
  filter(host_species != "Homo sapiens",
         host_genus != "Homo") %>%
  mutate(host_genus = ifelse(host_genus == "",
                             "Missing", host_genus)) %>%
  group_by(family) %>%
  summarise(perc_host_missing = sum(host_genus == "Missing") / n() * 100,
            perc_Y_missing = sum(is.na(Y)) / n() * 100,
            perc_M_missing = sum(is.na(M)) / n() * 100,
            perc_D_missing = sum(is.na(D)) / n() * 100,
            n_genomes = n()) %>%
  filter(n_genomes >= 100) %>%
  arrange(desc(perc_host_missing)) %>%
  filter(family %in% to_keep)

pal <- distinctColorPalette(n_distinct(plot_df$family))
pal_list <- setNames(pal, unique(plot_df$family))

plot_df %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family, y = perc_host_missing, fill = family)) +
  geom_bar(stat = "identity",
           color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.3, 
                                   size = 14,
                                   hjust = 1, 
                                   face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  scale_fill_manual(values = pal_list) +
  labs(x = "Viral family", y = "Unresolved genus (%)")

ggsave(str_glue("results/meta_summary_out/missing_genus_barplot.pdf"), 
       dpi = 600,
       width = 8, height = 4) 

plot_df2 <- plot_df %>%
  arrange(desc(perc_Y_missing))

plot_df2 %>%
  mutate(family = factor(family, unique(plot_df2$family))) %>%
  ggplot(aes(x = family, y = perc_Y_missing, fill = family)) +
  geom_bar(stat = "identity", 
           color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.3, 
                                   size = 14,
                                   hjust = 1, 
                                   face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  scale_fill_manual(values = pal_list) +
  labs(x = "Viral family", y = "Collection year missing (%)")

ggsave(str_glue("results/meta_summary_out/missing_dates_barplot.pdf"), 
       dpi = 600,
       width = 8, height = 4) 

meta %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, "")) %>%
  filter(host_species != "Homo sapiens",
         host_genus != "Homo") %>%
  filter(family == "Kolmioviridae") %>%
  filter(is.na(Y))
n_distinct(plot_df$family)
n_distinct(plot_df2$family)

plot_df3 <- meta %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, "")) %>%
  filter(host_species != "Homo sapiens",
         host_genus != "Homo") %>%
  mutate(host_genus = ifelse(is.na(host_genus) | host_genus == "",
                             "Missing", host_genus)) %>%
  filter(country != "") %>%
  group_by(country) %>%
  summarise(perc_host_missing = sum(host_genus == "Missing") / n() * 100,
            perc_Y_missing = sum(is.na(Y)) / n() * 100,
            perc_M_missing = sum(is.na(M)) / n() * 100,
            perc_D_missing = sum(is.na(D)) / n() * 100,
            n_genomes = n()) %>%
  mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
  mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
  dplyr::rename(region = country)

world_map <- map_data("world")

map_df <- subset(world_map, region != "Antarctica") %>%
  left_join(plot_df3 %>%
              mutate(missing_host = case_when(perc_host_missing < 25 ~ "<25%",
                                               perc_host_missing >= 25 & perc_host_missing <= 50 ~ "25-50%",
                                               perc_host_missing > 50 & perc_host_missing <= 75 ~ "50-75%",
                                               perc_host_missing > 75 ~ ">75%")) %>%
              mutate(missing_host = factor(missing_host, c("<25%", "25-50%", "50-75%", ">75%"))) %>%
              mutate(missing_year = case_when(perc_Y_missing < 25 ~ "<25%",
                                              perc_Y_missing >= 25 & perc_Y_missing <= 50 ~ "25-50%",
                                              perc_Y_missing > 50 & perc_Y_missing <= 75 ~ "50-75%",
                                              perc_Y_missing > 75 ~ ">75%")) %>%
              mutate(missing_year = factor(missing_year, c("<25%", "25-50%", "50-75%", ">75%"))))

ggplot() +
  geom_map(data = map_df, 
           map = map_df,
           aes(x = long, y = lat, map_id = region,
               fill = missing_year),
           size = 0.25,
           color = "black") +
  scale_fill_viridis_d(na.value = "grey") +
  # scale_fill_manual(values = c("steelblue", ""))
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Collection year missing")

ggsave(str_glue("results/meta_summary_out/missing_year_by_country.pdf"), 
       dpi = 600,
       width = 9, height = 5) 

ggplot() +
  geom_map(data = map_df, 
           map = map_df,
           aes(x = long, y = lat, map_id = region,
               fill = missing_host),
           size = 0.25,
           color = "black") +
  scale_fill_viridis_d(na.value = "grey") +
  # scale_fill_manual(values = c("steelblue", ""))
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Unresolved genus")

# plot_df3 %>% filter(region == "USA")
ggsave(str_glue("results/meta_summary_out/missing_genus_by_country.pdf"), 
       dpi = 600,
       width = 9, height = 5) 

require(UpSetR)
