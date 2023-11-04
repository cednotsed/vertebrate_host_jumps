rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv"))

dnds_df <- fread("results/dnds_out/dnds_out.diff_hosts.genus_counts.by_gene.csv") %>%
  mutate(ka_minus_ks = ka - ks,
         kaks = ka / ks) %>%
  filter(!is.na(kaks)) %>%
  dplyr::rename(clique_name = cluster)

# Get host counts
host_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  filter(host_genus != "") %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_genus))

genome_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(n_genomes = n_distinct(accession))

genome_length <- meta %>%
  dplyr::rename(clique_name = cluster) %>% 
  group_by(clique_name) %>%
  summarise(median_length = median(genome_length))

# Gene counts
dnds_df %>%
  separate(cluster, into = c("viral_family", "_")) %>%
  group_by(viral_family) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

dnds_df %>%
  filter(grepl("Coronaviridae", clique_name)) %>%
  group_by(is_jump) %>%
  summarise(n = n())
# Coronaviridae
parsed <- dnds_df %>%
  filter(grepl("Coronaviridae", clique_name)) %>%
  mutate(gene_type = case_when(grepl("ORF1|1a|1b|polyprotein|non-struct|ns", gene_name, ignore.case = T) ~ "                  ORF1ab",
                               grepl("matrix|hema|hemma|N protein|nucleo|membrane|M protein|envelope|E protein", gene_name, ignore.case = T) |
                                 gene_name %in% c("M", "E", "N", "HE")~ "Structural",
                               grepl("S protein|spike|glyco|surface", gene_name, ignore.case = T) |
                                 gene_name %in% c("S") ~ "Entry",
                               grepl("accessory|ORF|putative|i protein|protein i|N2|3a|3b|4b|4c|5a|5b|7a|sars6|6b|7b|8a|8b|9a|9b|protein 6|6 protein|protein 3|protein 7|protein 8", gene_name, ignore.case = T) &
                                 !grepl("ORF1|1a|1b|polyprotein|non-struct|ns", gene_name, ignore.case = T) ~ "Accessory")) %>%
  group_by(anc_name, tip_name, gene_type)

parsed %>% filter(is.na(gene_type)) %>% View()
jump_count <- deframe(parsed %>% 
                        group_by(cluster, anc_name, tip_name) %>%
                        summarise(n = n()) %>%
                        nrow())

cov_plt <- parsed %>% 
  mutate(gene_type = factor(gene_type, c("Entry", "Accessory",
                                         "Structural", "ORF1ab"))) %>%
  left_join(host_counts) %>%
  filter(!is.na(gene_type)) %>%
  ggplot(aes(x = gene_type, y = log10(kaks), fill = gene_type)) +
  theme_classic() +
  labs(x = "Product type", y = "log10(ka/ks)",
       title = str_glue("Coronaviridae: No. jumps (no subsampling) = {jump_count}")) +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  scale_fill_manual(values = c("dodgerblue3", "indianred3", 
                               "goldenrod2", "mediumpurple3")) +
  coord_flip() +
  geom_pwc() +
  facet_grid(rows = vars(is_jump)) +
  theme(legend.position = "none")

ggsave("results/dnds_out/CoV_gene_dnds.png", plot = cov_plt,
       width = 8, height = 5)

# Paramyxoviridae
parsed <- dnds_df %>%
  filter(grepl("Paramyx", clique_name)) %>%
  # distinct(gene_name) %>% View()
  mutate(gene_type = case_when(grepl("nucle|matirx|matrix|NP|membrane|M protein|nuleo", gene_name, ignore.case = T)  |
                                 gene_name %in% c("M", "N") ~ "Structural",
                               grepl("phosp|polymerase|large protein|L protein", gene_name, ignore.case = T) |
                                 gene_name %in% c("L", "P") ~ "Replication-associated",
                               grepl("receptor|haem|fusion|hema|glyco|HN", gene_name, ignore.case = T) |
                                 gene_name %in% c("S", "H", "F") ~ "Entry",
                               TRUE ~ "Accessory"))

jump_count <- deframe(parsed %>% 
                        # filter(kaks < 1000) %>%
                        group_by(clique_name, anc_name, tip_name) %>%
                        summarise(n = n()) %>%
                        nrow())
para_plt <- parsed %>% 
  mutate(gene_type = factor(gene_type, c("Entry", "Accessory", 
                                         "Structural", "Replication-associated"))) %>%
  left_join(host_counts) %>%
  filter(!is.na(gene_type)) %>%
  ggplot(aes(x = gene_type, y = log10(kaks), fill = gene_type)) +
  theme_classic() +
  labs(x = "Product type", y = "log10(ka/ks)",
       title = str_glue("Paramyxoviridae: No. jumps (no subsampling) = {jump_count}")) +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  scale_fill_manual(values = c("dodgerblue3", "indianred3", 
                               "goldenrod2", "mediumpurple3")) +
  coord_flip() +
  geom_pwc() +
  facet_grid(rows = vars(is_jump)) +
  theme(legend.position = "none")

ggsave("results/dnds_out/Paramyxoviridae_gene_dnds.png", plot = para_plt,
       width = 8, height = 5)


# Filoviridae
parsed <- dnds_df %>%
  filter(grepl("Filoviridae", cluster)) %>%
  # distinct(gene_name) %>% View()
  mutate(gene_type = case_when(grepl("membrane|nucle|NP|VP40|VP|24|30|40|matrix", gene_name, ignore.case = T)  |
                                 gene_name %in% c("NP") ~ "Structural",
                               grepl("polymerase", gene_name, ignore.case = T) |
                                 gene_name %in% c("L") ~ "Replication-associated",
                               grepl("glyco|GP", gene_name, ignore.case = T) &
                                 !grepl("small|sgp", gene_name, ignore.case = T) ~ "Entry",
                               TRUE ~ "Accessory"))

# parsed %>% filter(is.na(gene_type)) %>% distinct(gene_name) %>% View()
jump_count <- deframe(parsed %>% 
                        group_by(cluster, anc_name, tip_name) %>%
                        summarise(n = n()) %>%
                        nrow())
filo_plt <- parsed %>% 
  mutate(gene_type = factor(gene_type, c("Entry", "Accessory", 
                                         "Structural", "Replication-associated"))) %>%
  left_join(host_counts) %>%
  filter(!is.na(gene_type)) %>%
  ggplot(aes(x = gene_type, y = log10(kaks), fill = gene_type)) +
  theme_classic() +
  labs(x = "Product type", y = "log10(ka/ks)",
       title = str_glue("Filoviridae: No. jumps (no subsampling) = {jump_count}")) +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  scale_fill_manual(values = c("dodgerblue3", "indianred3", 
                               "goldenrod2", "mediumpurple3")) +
  coord_flip() +
  geom_pwc() +
  theme(legend.position = "none")

ggsave("results/dnds_out/Filoviridae_gene_dnds.png", plot = filo_plt,
       width = 8, height = 5)

# Filoviridae
parsed <- dnds_df %>%
  filter(grepl("Rhabdoviridae", cluster)) %>%
  # distinct(gene_name) %>% View()
  mutate(gene_type = case_when(grepl("max|matrix|nucl|M protein", gene_name, ignore.case = T)  |
                                 gene_name %in% c("NP", "N", "M", "N protein") ~ "Structural",
                               grepl("polymerase|phos|P protein|L protein|large", gene_name, ignore.case = T) |
                                 gene_name %in% c("L", "P") ~ "Replication-associated",
                               grepl("gyl|glyco|G protein", gene_name, ignore.case = T) |
                                 gene_name %in% c("G") ~ "Entry",
                               TRUE ~ "Accessory"))

# parsed %>% filter(is.na(gene_type)) %>% distinct(gene_name) %>% View()
jump_count <- deframe(parsed %>% 
                        group_by(cluster, anc_name, tip_name) %>%
                        summarise(n = n()) %>%
                        nrow())
rhab_plt <- parsed %>% 
  mutate(gene_type = factor(gene_type, c("Entry", "Accessory", 
                                         "Structural", "Replication-associated"))) %>%
  left_join(host_counts) %>%
  filter(!is.na(gene_type)) %>%
  ggplot(aes(x = gene_type, y = log10(kaks), fill = gene_type)) +
  theme_classic() +
  labs(x = "Product type", y = "log10(ka/ks)",
       title = str_glue("Rhabdoviridae: No. jumps (no subsampling) = {jump_count}")) +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  scale_fill_manual(values = c("dodgerblue3", "indianred3", 
                               "goldenrod2", "mediumpurple3")) +
  coord_flip() +
  geom_pwc() +
  theme(legend.position = "none")

ggsave("results/dnds_out/Rhabdoviridae_gene_dnds.png", plot = rhab_plt,
       width = 8, height = 5)
ggarrange(cov_plt, para_plt, filo_plt, rhab_plt)

ggsave("results/dnds_out/merged_family_gene_dnds.png", width = 12, height = 8)
