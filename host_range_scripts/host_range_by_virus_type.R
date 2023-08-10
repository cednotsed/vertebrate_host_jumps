rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.gt1000nt.220723.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
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

# Save family --> genome type metadata
meta_filt %>%
  distinct(family, genome_type) %>%
  fwrite("data/metadata/genome_type_metadata.csv")

genome_df <- genome_meta %>%
  group_by(family) %>%
  summarise(median_length = median(genome_length))

plot_df <- meta_filt %>%
  left_join(genome_df) %>%
  filter(is_vertebrate) %>%
  filter(!is.na(median_length)) %>% 
  group_by(family, genome_type, median_length) %>%
  summarise(n_hosts = n_distinct(host_order),
            n_genomes = n())

set.seed(66)
plot_df %>%
  ggplot(aes(x = log10(n_genomes), y = n_hosts, fill = genome_type)) +
  geom_point(
    # aes(size = log10(n_genomes)),
             color = "black",
             alpha = 1,
             size = 5,
             pch = 21) +
  geom_text_repel(aes(label = family), 
                  size = 3,
                  fontface = "italic", box.padding = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("mediumorchid3", "darkgoldenrod3", "darkolivegreen3",
                                "#2e4057", "darkorchid4", "cyan4",
                               "indianred3", "royalblue")) +
  labs(x = "No. of non-human sequences", 
       y = "No. of vertebrate orders",
       fill = "Genome type")

ggsave("results/host_range_out/host_orders_by_viral_family.png", height = 5, width = 8)
  

# plot_df %>%
#   filter(median_length < 50000) %>%
#   ggplot(aes(x = median_length, y = n_hosts, color = genome_type,
#              size = n_genomes)) +
#   geom_point()
  # geom_point(aes(size = log10(n_genomes))) +
  # geom_text_repel(aes(label = family))
# meta_filt %>%
#   group_by(host) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
