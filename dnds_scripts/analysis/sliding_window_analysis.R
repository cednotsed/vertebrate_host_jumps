rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)
require(MSA2dist)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv"))

genome_type <- fread("data/metadata/genome_type_metadata.csv")

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

dnds_list <- list.files("results/dnds_out/all_jumps.by_gene.temp_results/", full.names = T)
dnds_df <- foreach(file_name = dnds_list, .combine = c("bind_rows")) %do% {
  fread(file_name)
}

# Remove zeroes
dnds_filt <-  dnds_df

clique_counts <- dnds_filt %>%
  left_join(host_counts) %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- dnds_filt %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

cov_df <- merged_df %>%
  filter(grepl("Coronaviridae", clique_name))

cov_split <- cov_df %>% 
  # distinct(clique_name)
  mutate(gene_type = case_when(grepl("ORF1|1a|1b|polyprotein|non_struct|ns", gene_annot, ignore.case = T) ~ "ORF1ab",
                               grepl("internal|I_protein|i_protein|protein_i", gene_annot, ignore.case = T) ~ "I protein",
                               grepl("matrix|membrane|protein_M|M_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("M") ~ "Membrane",
                               grepl("N2|N_protein|nucelo|nucleo|nucloe|envelope|E_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("E") ~ "Envelope",
                               grepl("nucelo|nucleo|nucloe", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("N") ~ "Nucleocapsid",
                               # grepl("internal|I_protein|matrix|N_protein|nucelo|nucleo|nucloe|membrane|protein_M|M_protein|envelope|E_protein", gene_annot, ignore.case = T) |
                               #   gene_annot %in% c("M", "E", "N")~ "Structural",
                               grepl("S_protein|spike|glyco|surface|hema|hemma|esterase", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("S", "HE") ~ "Entry",
                               grepl("ORF3a|3a|3b|3c|protein_3", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("orf3a") ~ "ORF3",
                               grepl("4a|4b|4c", gene_annot, ignore.case = T) ~ "ORF4",
                               grepl("5a|5b|5c", gene_annot, ignore.case = T) ~ "ORF5",
                               grepl("sars6|6a|6b|protein_6|6_protein", gene_annot, ignore.case = T) ~ "ORF6",
                               grepl("7a|7b|protein_7", gene_annot, ignore.case = T) ~ "ORF7",
                               grepl("8a|8b|protein_8", gene_annot, ignore.case = T) ~ "ORF8",
                               grepl("9a|9b|protein_9", gene_annot, ignore.case = T) ~ "ORF9",
                               grepl("32_kda|32kda", gene_annot, ignore.case = T) ~ "32kDa protein")) %>%
  filter(!is.na(gene_type)) %>%
  filter(gene_type != "ORF9")

spike <- merged_df %>%
  filter(grepl("Coronaviridae", clique_name)) %>%
  filter((grepl("S_protein|spike|spike_glyco|surface|S1_", gene_annot, ignore.case = T) |
            gene_annot %in% c("S")) &
           !grepl("hema|9.5_kDa", gene_annot))

spike %>% distinct(gene_annot) %>% View()

# Check all tips have only one spike
spike %>% distinct(tip_name) %>% nrow() == spike  %>% nrow()

# cl <- makeCluster(10)
# registerDoParallel(cl)
# 
# spike_fna <- readDNAStringSet("results/dnds_out/family_gene_alignments/CoV_spike.fna")
# n_windows <- 20
# 
# spike_res_df <- foreach(i = seq(1, length(spike_fna), 2),
#                         .combine = "bind_rows",
#                         .packages = c("foreach", "MSA2dist", "Biostrings", "tidyverse")) %dopar% {
#   # i = 1
#   print(str_glue("sequence{i}"))
#   temp_aln <- spike_fna[i:(i + 1)]
#   seq_names <- str_split(names(temp_aln), "\\|", simplify = T)[, 1]
#   
#   aln_length <- unique(width(temp_aln))
#   window_size <- floor(aln_length / (n_windows + 1))
#   window_size <- window_size - (window_size %% 3)
#   
#   idx_list <- seq(1, aln_length, window_size + 1)[1:n_windows]
#   remainder <- (idx_list[n_windows] + window_size) %% 3
#   idx_list[n_windows] <- idx_list[n_windows] - remainder
#   
#   if(as.character(temp_aln[[1]]) == as.character(temp_aln[[2]])) {
#     print(str_glue("idnetical sequence{i}"))
#     return(tibble(anc_name = seq_names[2],
#                   tip_name = seq_names[1],
#                   ka = 0,
#                   ks = 0,
#                   window = seq(n_windows)))
#   } else {
#     crumbs <- foreach(j = seq(length(idx_list))) %do% {
#       # j = 1
#       # .GlobalEnv$idx_list <- idx_list
#       # .GlobalEnv$temp_aln <- temp_aln
#       # .GlobalEnv$window_size <- window_size
#       start_pos <- idx_list[j]
#       window_aln <- subseq(temp_aln, 
#                            start = start_pos, 
#                            end = (start_pos + window_size - 1))
#       dnds <- dnastring2kaks(window_aln,
#                              model = "Li")
#       ka <- as.numeric(dnds[, "ka"])
#       ks <- as.numeric(dnds[, "ks"])
#       
#       return(tibble(anc_name = seq_names[2],
#                     tip_name = seq_names[1],
#                     ka = ka,
#                     ks = ks,
#                     window = j))
#     }
#     return(bind_rows(crumbs))
#   }
# }
# 
# stopCluster(cl)

spike_res_df <- fread("results/dnds_out/family_gene_alignments/spike_res.csv")

plot_df <- spike %>%
  distinct(clique_name, anc_name, tip_name, anc_state, tip_state, is_jump, gene_length) %>%
  left_join(spike_res_df) %>%
  arrange(clique_name, gene_length)

# Ka only
plot_df %>%
  filter(clique_name == "Coronaviridae_12") %>%
  mutate(tip_name = factor(tip_name, unique(plot_df$tip_name))) %>%
  ggplot(aes(x = window, y = tip_name, fill = log10(ka))) +
  geom_tile(na.rm = T) +
  scale_fill_gradient2(low = "blue", 
                       mid = "#FFF9c4",
                       high = "red", 
                       na.value = NA,
                       midpoint = -1.7) +
  facet_grid(rows = vars(is_jump),
             scales = "free") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Spike window", y = "Accession", fill = "log10(Ka)",
       title = "Coronaviridae_12 (i.e., SARS-CoV-2-related viruses)")

ggsave("results/dnds_out/CoV_spike_dn.png", width = 8, height = 3)

plot_df %>%
  filter(clique_name == "Coronaviridae_12") %>%
  mutate(tip_name = factor(tip_name, unique(plot_df$tip_name))) %>%
  ggplot(aes(x = window, y = tip_name, fill = log10(ks))) +
  geom_tile(na.rm = T) +
  scale_fill_gradient2(low = "blue", 
                       mid = "#FFF9c4",
                       high = "red", 
                       na.value = NA,
                       midpoint = -1.4) +
  facet_grid(rows = vars(is_jump),
             scales = "free") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Spike window", y = "Accession", fill = "log10(Ks)",
       title = "Coronaviridae_12 (i.e., SARS-CoV-2-related viruses)")

ggsave("results/dnds_out/CoV_spike_ds.png", width = 8, height = 3)

plot_df %>%
  filter(clique_name == "Coronaviridae_12") %>%
  # filter(is_jump) %>% View()
  # filter(is_jump) %>% 
  # mutate(kaks = case_when(ka == 0 & ks != 0 ~ 1e-6,
  #                         ka != 0 & ks == 0 ~ 1e6,
  #                         ka != 0 & ks != 0 ~ (ka / ks),
  #                         ka == 0 & ks == 0 ~ NA)) %>%
  # summarise(range(ka, na.rm = T))
  mutate(tip_name = factor(tip_name, unique(plot_df$tip_name))) %>%
  # filter(tip_state == "Homo sapiens") %>%
  ggplot(aes(x = window, y = tip_name, fill = ka / ks)) +
  geom_tile(na.rm = T) +
  scale_fill_gradient2(low = "blue", 
                       mid = "#FFF9c4", 
                       high = "red", 
                       na.value = NA,
                       midpoint = -0.35) +
  facet_grid(rows = vars(is_jump),
             scales = "free") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Spike window", y = "Accession", fill = "log10(Ka/ks)",
       title = "Coronaviridae_12 (i.e., SARS-CoV-2-related viruses)")

ggsave("results/dnds_out/CoV_spike_dnds.png", width = 8, height = 5)


meta %>%
  filter(grepl("Coronaviridae", cluster)) %>%
  group_by(cluster) %>%
  filter(cluster %in% spike$clique_name) %>%
  distinct(species) %>% View()
