rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)
require(Hmisc)
require(maps)

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv"))

genome_meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")
ictv_meta <- fread("data/metadata/ICTV_Master_Species_List_2022_MSL38.v2.060823.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))

# Proportion of viruses by type

type_counts <- meta %>%
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


# Proportion of SARS-CoV-2
meta %>%
  summarise(prop = sum(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2",
                       genbank_title,
                       ignore.case = T) / n()))

meta %>%
  summarise(prop = sum(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2",
                             genbank_title,
                             ignore.case = T) / n()))
meta %>%
  nrow()
# Proportion of human sequences
meta %>%
  filter(!is.na(taxid)) %>%
  summarise(n = sum(taxid == 9606| host_species == "Homo sapiens") / n())

# Proportion of animal orders
meta %>%
  filter(is_vertebrate) %>%
  distinct(host_genus) %>% View()

# https://amphibiansoftheworld.amnh.org/ accessed 4/8/23
amphibian_families <- c("Family: Allophrynidae (3 sp.)",
                        "Family: Alsodidae (30 sp.)",
                        "Family: Alytidae (12 sp.)",
                        "Family: Arthroleptidae (152 sp.)",
                        "Family: Ascaphidae (2 sp.)",
                        "Family: Batrachylidae (12 sp.)",
                        "Family: Bombinatoridae (9 sp.)",
                        "Family: Brachycephalidae (79 sp.)",
                        "Family: Ceuthomantidae (4 sp.)",
                        "Family: Craugastoridae (129 sp.)",
                        "Family: Eleutherodactylidae (241 sp.)",
                        "Family: Strabomantidae (797 sp.)",
                        "Family: Brevicipitidae (37 sp.)",
                        "Family: Bufonidae (647 sp.)",
                        "Family: Calyptocephalellidae (5 sp.)",
                        "Family: Centrolenidae (164 sp.)",
                        "Family: Ceratobatrachidae (102 sp.)",
                        "Family: Ceratophryidae (12 sp.)",
                        "Family: Conrauidae (8 sp.)",
                        "Family: Cycloramphidae (37 sp.)",
                        "Family: Aromobatidae (136 sp.)",
                        "Family: Dendrobatidae (205 sp.)",
                        "Family: Dicroglossidae (218 sp.)",
                        "Family: Heleophrynidae (7 sp.)",
                        "Family: Hemiphractidae (121 sp.)",
                        "Family: Hemisotidae (9 sp.)",
                        "Family: Hylidae (1049 sp.)",
                        "Family: Hylodidae (48 sp.)",
                        "Family: Hyperoliidae (227 sp.)",
                        "Family: Leiopelmatidae (3 sp.)",
                        "Family: Leptodactylidae (234 sp.)",
                        "Family: Mantellidae (271 sp.)",
                        "Family: Megophryidae (315 sp.)",
                        "Family: Micrixalidae (24 sp.)",
                        "Family: Microhylidae (744 sp.)",
                        "Family: Limnodynastidae (44 sp.)",
                        "Family: Myobatrachidae (91 sp.)",
                        "Family: Nasikabatrachidae (2 sp.)",
                        "Family: Nyctibatrachidae (37 sp.)",
                        "Family: Odontobatrachidae (5 sp.)",
                        "Family: Odontophrynidae (55 sp.)",
                        "Family: Pelobatidae (6 sp.)",
                        "Family: Pelodytidae (4 sp.)",
                        "Family: Petropedetidae (13 sp.)",
                        "Family: Phrynobatrachidae (96 sp.)",
                        "Family: Pipidae (41 sp.)",
                        "Family: Ptychadenidae (63 sp.)",
                        "Family: Pyxicephalidae (82 sp.)",
                        "Family: Ranidae (446 sp.)",
                        "Family: Ranixalidae (18 sp.)",
                        "Family: Rhacophoridae (455 sp.)",
                        "Family: Rhinodermatidae (3 sp.)",
                        "Family: Rhinophrynidae (1 sp.)",
                        "Family: Scaphiopodidae (7 sp.)",
                        "Family: Sooglossidae (4 sp.)",
                        "Family: Telmatobiidae (61 sp.)",
                        "Family: Ambystomatidae (30 sp.)",
                        "Family: Amphiumidae (3 sp.)",
                        "Family: Cryptobranchidae (5 sp.)",
                        "Family: Hynobiidae (99 sp.)",
                        "Family: Plethodontidae (519 sp.)",
                        "Family: Proteidae (9 sp.)",
                        "Family: Rhyacotritonidae (4 sp.)",
                        "Family: Salamandridae (146 sp.)",
                        "Family: Sirenidae (7 sp.)",
                        "Family: Caeciliidae (49 sp.)",
                        "Family: Chikilidae (4 sp.)",
                        "Family: Dermophiidae (15 sp.)",
                        "Family: Grandisoniidae (24 sp.)",
                        "Family: Herpelidae (10 sp.)",
                        "Family: Ichthyophiidae (57 sp.)",
                        "Family: Rhinatrematidae (14 sp.)",
                        "Family: Scolecomorphidae (6 sp.)",
                        "Family: Siphonopidae (28 sp.)",
                        "Family: Typhlonectidae (14 sp.)")

amphibian_families <- str_split(amphibian_families, " ", simplify = T)[, 2]
amphibian_families <- amphibian_families[amphibian_families != ""]

bird_families <- fread("data/ecology/host_orders/NEW_eBird-Clements-v2022-integrated-checklist-October-2022_040823.csv")
bird_families <- unique(bird_families$FAMILY[bird_families$FAMILY != ""])
bird_families <- str_split(bird_families, " ", simplify = T)[, 1]

mammalian_families <- fread("data/ecology/host_orders/MDD_v1.11_6649species_050823.csv")$family
mammalian_families <- unique(capitalize(tolower(mammalian_families)))
mammalian_families <- mammalian_families[mammalian_families != ""]

reptile_families <- fread("data/ecology/host_orders/reptile_checklist_2023_07_050823.csv")
reptile_families <- unique(reptile_families$Family)
reptile_families <- reptile_families[reptile_families != ""]

fish_families <- fread("data/ecology/host_orders/Eschmeyer_CoF.050823.csv")
fish_families <- unique(fish_families$V1)
fish_families <- str_trim(fish_families, "both")
fish_families <- fish_families[endsWith(fish_families, "dae")]
fish_families <- fish_families[fish_families != ""]

sampled_families <- deframe(meta %>%
    filter(!is.na(host_family),
           host_family != "") %>%
    group_by(host_family) %>%
    summarise(n = n()) %>%
    filter(n >= 10) %>%
    select(host_family))

prop_sequenced <- c(sum(amphibian_families %in% sampled_families) / length(amphibian_families),
                    sum(bird_families %in% sampled_families) / length(bird_families),
                    sum(mammalian_families %in% sampled_families) / length(mammalian_families),
                    sum(reptile_families %in% sampled_families) / length(reptile_families),
                    sum(fish_families %in% sampled_families) / length(fish_families))
plot_df <- tibble(host_groups = c("Amphibians", "Birds", "Mammals", "Reptiles", "Fish"), 
       prop_sequenced = prop_sequenced) %>%
  arrange(desc(prop_sequenced))

plot_df %>%
  mutate(host_groups = factor(host_groups, plot_df$host_groups)) %>%
  ggplot() +
  geom_bar(stat = "identity", 
           aes(x = 100, y = host_groups),
           lty = "dashed",
           fill = "grey95",
           color = "black") +
  geom_bar(aes(x = prop_sequenced * 100, 
               y = host_groups,
               fill = host_groups),
           stat = "identity",
           color = "black") +
  scale_fill_manual(values = c("indianred3", "darkseagreen4", 
                               "goldenrod2", "mediumpurple3",
                               "dodgerblue4")) +
  theme_classic() +
  labs(x = "Perc. of families with >= 10 viral sequences", 
       y = "Vertebrate hosts") +
  xlim(0, 100) +
  geom_text(aes(x = prop_sequenced * 100, 
                y = host_groups,
                label = paste0(round(prop_sequenced * 100, 0), "%")),
            hjust = -0.8) +
  theme(legend.position = "none")

ggsave("results/meta_summary_out/prop_sequenced.pdf", dpi = 600, width = 6, height = 3)


count_df <- meta %>%
  # filter(!is.na(host_genus)) %>%
  # filter(host_genus != "") %>%
  # filter(!is.na(host_genus)) %>%
  mutate(host_genus = ifelse(is.na(host_genus) | host_genus == "",
                             "Missing", host_genus)) %>%
  group_by(host_genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(prop = n / sum(n))

count_filt <- count_df %>%
  head(6)

count_df %>% 
  ungroup() %>%
  filter(host_genus != "Homo") %>%
  mutate(prop = n / sum(n))

meta %>%
  filter(host_genus != "Homo" & host != "Homo sapiens") %>%
  summarise(missing = sum(collection_date == "") / n())

plot_df2 <- meta %>%
  mutate(host_genus = ifelse(is.na(host_genus) | host_genus == "",
                             "Missing", host_genus)) %>%
  mutate(host_genus = ifelse(host_genus %in% count_filt$host_genus, 
                             host_genus, 
                             "Others"))

plot_df2 %>%
  ggplot(aes(x = 1, fill = host_genus)) +
  geom_bar(aes(y = ..count..),
           color = "black") +
  coord_polar("y") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = c("steelblue3", "firebrick", "grey80",
                               "gold1", "grey35", "darkseagreen",
                               "lightpink2"))

ggsave("results/meta_summary_out/prop_by_genus.pdf", dpi = 600,
       height = 3, width = 5)



plot_df3 <- meta %>%
  filter(is_vertebrate) %>%
  filter(host != "") %>%
  filter(host_genus != "Homo", host_family != "Hominidae") %>%
  filter(country != "") %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
  mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
  rename(region = country)

world_map <- map_data("world")
world_map <- subset(world_map, region != "Antarctica") %>%
  left_join(plot_df3) %>%
  mutate(n = ifelse(is.na(n)|n == 0, 1, n))

# plot_df3 %>%
#   distinct(region) %>%
#   filter(grepl("Congo", region))
#   filter(grepl("")) %>%
#   distinct(region)
# plot_df3 %>%
#   distinct(region) %>%
#   View()
world_map %>% filter(grepl("Turk", region)) %>% distinct(region)
meta %>%
  filter(is_vertebrate) %>%
  filter(grepl("South Sudan", country))

ggplot() +
  geom_map(data = world_map, 
           map = world_map,
           aes(x = long, y = lat, map_id = region,
               fill = log10(n)),
           size = 0.25,
           color = "black") +
  scale_fill_gradient2(low = "blue", 
                       mid = "#FFF9c4",
                       high = "red",
                       midpoint = 2.5) +
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
  labs(fill = "log10(Non-human viral sequences)")

ggsave("results/meta_summary_out/world_map_surveillance_effort.pdf", 
       dpi = 600,
       width = 9, height = 5)  


# ICTV correspondence
ictv <- fread("data/metadata/ICTV_Master_Species_List_2022_MSL38.v2.060823.csv")
genome_meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")

genome_meta %>%
  filter(species != "") %>%
  distinct(species) %>%
  # filter(!(species %in% ictv$Species))
  summarise(prop = sum(species %in% ictv$Species) / n())

clique_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")

n_distinct(clique_meta$cluster)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  mutate(is_jump = T)

nrow(jump_df)
n_distinct(jump_df$clique_name)
