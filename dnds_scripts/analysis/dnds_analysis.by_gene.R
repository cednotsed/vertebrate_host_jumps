rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)

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
dnds_filt <-  dnds_df %>%
  mutate(ks = ifelse(ks == 0 | ks < 0, 0, ks),
         ka = ifelse(ka == 0, 0, ka)) %>%
  filter(ka != 0 & ks != 0) %>%
  mutate(kaks = ka / ks)

clique_counts <- dnds_filt %>%
  left_join(host_counts) %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- dnds_filt %>%
  group_by(clique_name, is_jump) %>%
  summarise(min_kaks = min(kaks),
            min_dist = min(patristic_dist)) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

parsed <- dnds_df %>%
  mutate(is_polymerase = ifelse(grepl("polymerase|repli|PB|rdrp", gene_name, ignore.case = T),
                                   T, F)) %>%
  mutate(is_capsid = ifelse(grepl("capsid|envelope|spike|surface", gene_name, ignore.case = T) &
                              !grepl("encapsidation", gene_name, ignore.case = T),
                            T, F))

parsed <- dnds_df %>%
  mutate(gene_type = case_when(
                              # Polymerase
                               (grepl("rep protein|polymerase|repli|PB|rdrp|DNA pol|helicase|SPPV_035|SPPV_039|SPPV_055|SPPV_069|SPPV_071|SPPV_085|SPPV_096|SPPV_116|SPPV_119|SPPV_032|SPPV_068", gene_name, ignore.case = T) |
                                 gene_name %in% c("Rep", "REP", "rep", 
                                                  "putative rep")) |
                                 (grepl("Adenoviridae", gene_name, ignore.case = T) |
                                    gene_name %in% c("pol")) |
                                 (grepl("Rhabdoviridae", cluster) &
                                    grepl("L protein", gene_name, ignore.case = T)) |
                                 (grepl("Rhabdoviridae", cluster) &
                                    gene_name %in% c("L")) ~ "Polymerase",
                               # Surface proteins
                               (grepl("hemag|haem|putative cap|capsid|envelope|env|spike|surface|membrane", gene_name, ignore.case = T) &
                                 !grepl("encapsidation", gene_name, ignore.case = T)) |
                                 
                                 (gene_name %in% c("S protein", "large S protein", "M protein")) |
                                 
                                 (grepl("Circoviridae|Smacoviridae", cluster) &
                                  grepl("cap|coat", gene_name, ignore.case = T)) |
                                 
                                 (grepl("Adenoviridae", cluster) &
                                  grepl("PIX|PIIIa|IV-|PV|hexon|penton|fiber|PIV|IVa2|core|mu px|IX", gene_name, ignore.case = T)) |
                                 
                                 (grepl("Rhabdoviridae", cluster) &
                                  grepl("gylco|nucleoprotein|glycoprotein|matrix|M1|G protein|N protein|M protein", gene_name, ignore.case = T)) |
                                 
                                 (grepl("Rhabdoviridae", cluster) &
                                    gene_name %in% c("M", "NP", "G")) |
                                 
                                 (grepl("Coronaviridae", cluster) &
                                    grepl("HE|N2|nucleo|glycoprotein|membrane|M protein|N protein|S protein", gene_name, ignore.case = T)) |
                                 
                                 (grepl("Filoviridae", cluster) &
                                    grepl("NP|nucleo|glyco|matrix|VP40|VP35|VP24|VP30", gene_name, ignore.case = T)) |
                                 
                                 (grepl("Coronaviridae", cluster) &
                                    gene_name %in% c("S", "M", "N", 
                                                     "E", "E protein")) |
                                 
                                 (grepl("Coronaviridae", cluster) &
                                    grepl("Enelope", gene_name, ignore.case = T)) |
                                 
                                 (grepl("Circoviridae", cluster) &
                                    gene_name %in% c("CP")) |
                                 
                                 (grepl("Togaviridae|Filoviridae|Parvoviridae|Poxviridae", cluster) &
                                  gene_name %in% c("structural protein", "structural protein precursor",
                                                   "structural polyprotein", "truncated structural polyprotein",
                                                   "Structural polyprotein", "structural glycoprotein",
                                                   "virion structural protein", "structural pplyprotein",
                                                   "structural protein 2", "structural protein 1",
                                                   "Structural protein", "structural polyprotein precursor")) |
                                 
                                 (grepl("Anelloviridae|Parvoviridae|Polyomaviridae", cluster) & 
                                  grepl("vp1", gene_name, ignore.case = T)) ~ "Structural proteins",
                              # Non-structural proteins
                              (grepl("factor|oxidase|dUTPase|NTPase|reductase|phosphatase|mutase|kinase|peptidase|protease|proteinase|ligase|transcriptase|lipase|glycosylase|NS protein|nonstructural|non-structural|non structural|nsp1|nsp3|nsp7|ns1|ns2|ns3|ns4|ns5|ns6|ns7|ns8|ns9", gene_name, ignore.case = T)) |
                                
                                (grepl("Coronaviridae", cluster) &
                                   grepl("ORF|7a|3b", gene_name, ignore.case = T) &
                                   !grepl("ORF1", gene_name, ignore.case = T) ) ~ "Non-structural",
                              # Unknown proteins 
                              TRUE ~ "Unknown"))
                               # # Unknown proteins 
                               # grepl("kda|unknown|hypothetical", gene_name, ignore.case = T) ~ "Unknown"))
  # filter(is.na(gene_type))

parsed %>%
  separate(cluster, into = "viral_family", "_") %>%
  # filter(grepl("surface|spike", gene_name, ignore.case = T)) %>%
  distinct(viral_family, gene_name, tip_name) %>% View()

res_filt <- res_df %>%
  inner_join(parsed)

res_filt %>%
  ggplot(aes(x = gene_type, y = log10(kaks))) +
  geom_boxplot()

host_counts <- meta %>%
  group_by(cluster) %>%
  summarise(n_hosts = n_distinct(host_order),
            n_genomes = n_distinct(accession)) %>%
  arrange(desc(n_hosts))

merged_df <- res_filt %>%
  left_join(host_counts)

merged_df %>%
  filter(gene_type == "Non-structural") %>%
  ggplot(aes(x = n_hosts, y = kaks)) +
  geom_point()
# Hosts vs ka - ks
corr_test <- cor.test(merged_df$n_hosts, merged_df$ka_minus_ks)
r_val <- signif(corr_test$estimate, 3)
p_val <- signif(corr_test$p.value, 3)

merged_df %>%
  ggplot(aes(x = n_hosts, y = ka_minus_ks)) +
  # geom_smooth(color = "black",
  #             method = "lm",
  #             fill = "tan1") +
  geom_point(color = "steelblue4", alpha = 0.5) +
  theme_classic() +
  # xlim(0, 3.5) +
  labs(x = "No. host orders", y = "Ka - Ks",
       title = str_glue("Pearson's R = {r_val}, p = {p_val}"))

# merged_filt %>%
#   ggplot(aes(ka_minus_ks, kaks)) +
#   geom_point()
# ggsave("results/dnds_out/hosts_vs_ka_minus_ks.png", width = 5, height = 3)

# Hosts vs ka/ks
merged_filt <- merged_df %>%
  filter(!is.infinite(kaks))
corr_test <- cor.test(merged_filt$n_hosts, merged_filt$kaks)
r_val <- signif(corr_test$estimate, 3)
p_val <- signif(corr_test$p.value, 3)

merged_filt %>%
  ggplot(aes(x = n_hosts, y = kaks)) +
  # geom_smooth(color = "black",
  #             # method = "lm",
  #             fill = "tan1") +
  geom_point(color = "steelblue4", alpha = 0.5) +
  theme_classic() +
  # xlim(0, 3.5) +
  labs(x = "No. host orders", y = "ka / ks",
       title = str_glue("Pearson's R = {r_val}, p = {p_val}"))

# ggsave("results/dnds_out/hosts_vs_kaks.png", width = 5, height = 3)

# corr_test <- cor.test(log10(merged_filt$n_genomes), merged_filt$kaks)
# r_val <- signif(corr_test$estimate, 3)
# p_val <- signif(corr_test$p.value, 3)
# 
# merged_filt %>%
#   ggplot(aes(x = log10(n_genomes), y = kaks)) +
#   geom_smooth(color = "black",
#               method = "lm",
#               fill = "tan1") +
#   geom_point(color = "steelblue4", alpha = 0.5) +
#   theme_classic() +
#   # xlim(0, 3.5) +
#   labs(x = "Log10(No. genomes)", y = "Min. patristic distance",
#        title = str_glue("Pearson's R = {r_val}, p = {p_val}"))
# 
# ggsave("results/mutational_load_out/kaks_vs_sampling.png", width = 5, height = 3)

# Without genome correction
lr1 <- glm(n_hosts ~ kaks,
           merged_filt,
           family = poisson())

lr2 <- glm(n_hosts ~ log10(kaks),
           merged_filt,
           family = poisson())

comp <- compareGLM(lr1, lr2)
best_model_idx <- deframe(comp$Fit.criteria %>%
                            rownames_to_column("model") %>%
                            arrange(AIC) %>%
                            select(model) %>%
                            head(1))

best_model_formula <- deframe(comp$Models %>%
                                as_tibble() %>%
                                rownames_to_column("model") %>%
                                filter(model == best_model_idx) %>%
                                select(Formula))

best_model <- get(paste0("lr", best_model_idx))
obs_effect <- coef(best_model)[[2]]

summary(best_model)

# Permutate uncorrected model
perm_no_correction <- permute_coeff(best_model_formula, 2, merged_filt)
p_val <- perm_no_correction %>%
  summarise(p_value = sum(coefficient <= obs_effect) / n())

perm_no_correction %>%
  ggplot(aes(x = coefficient)) +
  geom_histogram(fill = "steelblue4",
                 bins = 50) +
  geom_vline(xintercept = obs_effect,
             lty = "dashed",
             color = "red") +
  theme_classic() +
  labs(x = "Slope estimate", y = "No. permutations",
       title = str_glue("Perm. test p-value = {p_val}"))

ggsave("results/dnds_out/perm_test_slope.no_correction.png", width = 5, height = 3)

