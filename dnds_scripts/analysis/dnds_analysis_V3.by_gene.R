rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)
require(rstatix)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")

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

dnds_df <- foreach(file_name = dnds_list, .combine = "bind_rows") %do% {
  if(file.size(file_name) > 246) {
    fread(file_name)
  }
}

overall_dnds <- fread("results/dnds_out/all_jumps.dnds.diff_hosts.genus_counts.csv") %>%
  mutate(kaks = ka / ks) %>%
  select(anc_name, tip_name, overall_kaks = kaks)

type_counts <- dnds_df %>%
  group_by(clique_name) %>%
  summarise(n = n_distinct(is_jump)) %>%
  filter(n == 2)

# # Randomly select lineage
# iter_df <- dnds_df %>%
#   distinct(anc_name, tip_state, clique_name, is_jump)
# 
# set.seed(66)
# iter_filt <- foreach(i = seq(nrow(iter_df)), .combine = "c") %do% {
#   # i = 1
#   row <- iter_df[i, ]
#   temp_filt <- dnds_df %>%
#     filter(anc_name == row$anc_name,
#            tip_state == row$tip_state,
#            clique_name == row$clique_name,
#            is_jump == row$is_jump) %>%
#     distinct(tip_name) %>%
#     sample_n(1, replace = F)
# 
#   return(temp_filt$tip_name)
# }

dnds_filt <- dnds_df %>%
  # filter(tip_name %in% iter_filt) %>%
  mutate(ka = ifelse(ka <= 0, 0, ka),
         ks = ifelse(ks <= 0, 0, ks)) %>% 
  filter(ka != 0 & ks != 0) %>%
  # mutate(ka = ka + 1,
         # ks = ks + 1) %>%
  mutate(log_kaks = log(ka / ks))

clique_counts <- dnds_filt %>%
  left_join(host_counts) %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- dnds_filt %>%
  left_join(genome_counts)

cov_df <- merged_df %>%
  filter(grepl("Coronaviridae", clique_name)) %>%
  mutate(gene_type = case_when(grepl("replicase|1a|1b|polyprotein|polymerase", gene_annot, ignore.case = T) ~ "Replication-associated",
                               grepl("truncated_E|matrix|N_protein|nucelo|nucleo|nucloe|membran|enelope|protein_M|M_protein|envelope", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("M", "E", "N", "E_protein")~ "Structural",
                               grepl("S1_glyco|S_protein|spike|surface|hema|hemma|esterase", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("S", "HE", "HE_protein") ~ "Entry", 
                               # grepl("ORF2|ORF3|ORF_3|ORF4|ORF5|ORF_5|ORF6|ORF7|ORF_7|ORF8|ORF_X|ORF10|protein_10|protein_14|13|small_virion|unknown|hypothe|nonstructural|putative|non_struct|kda|internal|I_protein|NS|accessory|i_protein|protein_i|3c|N2|3a|3b|4a|4b|4c|5a|5b|7a|sars6|6b|7b|8a|8b|9a|9b|protein_6|6_protein|protein_3|protein_7|protein_8", gene_annot, ignore.case = T)|
                               #   gene_annot == "I" ~ "Accessory"
                               TRUE ~ "Accessory"))

# jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv")
# jump_df %>%
#   filter(is_jump) %>%
#   group_by(family) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

# cov_df %>% distinct(gene_annot, gene_type) %>% 
#   arrange(gene_type) %>% View()

para_df <- merged_df %>%
  filter(grepl("Paramyx", clique_name)) %>%
  mutate(gene_type = case_when(grepl("nucle|matirx|matrix|NP|membrane|M_protein|nuleo", gene_annot, ignore.case = T)  |
                                 gene_annot %in% c("M", "N", "N_protein") ~ "Structural",
                               grepl("Phosp|transcriptase|rdrp|polymerase|large_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("L", "P", "P_protein", "L_protein") ~ "Replication-associated",
                               grepl("receptor|haem|fusion|hema|glyco|HN", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("S", "H", "F", "F_protein", "H_protein") ~ "Entry",
                               TRUE ~ "Accessory"))

rhab_df <- merged_df %>%
  filter(grepl("Rhabdoviridae", clique_name)) %>%
  mutate(gene_type = case_when(grepl("max|matrix|nucl|M_protein|M1|M2", gene_annot, ignore.case = T)  |
                                 gene_annot %in% c("NP", "N", "M", "N_protein", "MP") ~ "Structural",
                               grepl("polymerase|phos|P_protein|L_protein|large", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("L", "P", "LP", "PP") ~ "Replication-associated",
                               grepl("gyl|glyco|G_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("G", "GP") ~ "Entry",
                               TRUE ~ "Accessory"))

pneumo_df <- merged_df %>%
  filter(grepl("Pneumoviridae", clique_name)) %>%
  mutate(gene_type = case_when(grepl("max|nucl|M_protein|M1", gene_annot, ignore.case = T)  |
                                 gene_annot %in% c("NP", "N", "M", "N_protein",
                                                   "MP", "matrix_protein", "matrix") ~ "Structural",
                               grepl("polymerase|phos|P_protein|large|M2|matrix_protein_2", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("L", "P", "LP", "PP") ~ "Replication-associated",
                               grepl("attachment|fusion|gyl|glyco|G_protein", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("G", "GP", "F", "glycoprotein") ~ "Entry",
                               TRUE ~ "Accessory"))

circo_df <- merged_df %>% 
  filter(grepl("Circoviridae", clique_name)) %>%
  mutate(gene_type = case_when(grepl("BAQ95397.1|Caspid|cap|coat|CP|ORFC2", gene_annot, ignore.case = T)  |
                                 gene_annot %in% c("cap") ~ "Structural",
                               grepl("rep|BAQ95396.1", gene_annot, ignore.case = T) |
                                 gene_annot %in% c("Rep") ~ "Replication-associated",
                               TRUE ~ "Accessory"))

# cov_df %>%
#   filter(!is_jump) %>%
#   filter(log10(ka) < 0.5, log10(ks) < 0.5) %>%
#   ggplot(aes(x = log10(ka), y = log10(ks), color = gene_type)) +
#   geom_point() +
#   geom_smooth(method = "lm")

# Write source data
combined_df <- bind_rows(cov_df, para_df, rhab_df, circo_df) %>%
  mutate(gene_type = ifelse(gene_type == "Accessory", "Auxiliary", gene_type))

# fwrite(combined_df, "results/dnds_out/parsed_dnds.by_gene.csv")

dat_list <- list(Coronaviridae = cov_df, Paramyxoviridae = para_df,
                 Rhabdoviridae = rhab_df, Circoviridae = circo_df)

for(i in seq(length(dat_list))) {
  dat <- dat_list[[i]]
  fam_name <- names(dat_list)[i]
  
  corr_lr <- lm(log_kaks ~ clique_name,
                data = dat)

  dat <- dat %>% 
    add_column(corr_effect = corr_lr$residuals) %>%
    mutate(gene_type = factor(gene_type, c("Structural", "Entry", "Replication-associated",
                                           "Accessory")))
  
  #Get model effects
  lr <- lm(corr_effect ~ interaction(gene_type, is_jump),
           data = dat)
  
  lr_sum <- summary(lr)$coefficients
  
  rownames(lr_sum) <- gsub("interaction\\(gene_type\\, is_jump\\)", "", 
                           rownames(lr_sum))
  
  # Calculate difference in estimates
  deg_freedom <- nrow(dat) - nrow(lr_sum)
  
  model_res <- foreach(gene_type = unique(dat$gene_type),
                       .combine = "bind_rows") %do% {
                         if(gene_type != "Structural") {
                           jump_est <- lr_sum[str_glue("{gene_type}.TRUE"),]
                           nonjump_est <- lr_sum[str_glue("{gene_type}.FALSE"),]
                         } else {
                           jump_est <- lr_sum[str_glue("{gene_type}.TRUE"),]
                           nonjump_est <- lr_sum[str_glue("(Intercept)"),]
                         }
                         effect <- jump_est[["Estimate"]] - nonjump_est[["Estimate"]]
                         se_effect <- sqrt(jump_est[["Std. Error"]]^2 + nonjump_est[["Std. Error"]]^2)
                         t <- effect / se_effect
                         p_val <- ifelse(t > 0, 
                                         pt(t, df = deg_freedom, lower.tail = F),
                                         pt(t, df = deg_freedom, lower.tail = T))
                         
                         return(tibble(gene_type = gene_type,
                                       effect = effect, 
                                       std_error = se_effect, 
                                       t_val = t, 
                                       p_val = p_val))
                       }
  
  parsed_models <- model_res %>%
    mutate(effect = signif(effect, 3),
           p_val = signif(p_val, 1)) %>%
    mutate(deg_freedom = deg_freedom) %>%
    arrange(effect) 
  
  # Save results
  sink(str_glue("results/dnds_out/family_plots/{fam_name}.effects.csv"))
  print(parsed_models)
  sink()
  
  # Plot effects
  pal <- setNames(c("dodgerblue3", "indianred3", "goldenrod3", "mediumpurple3"),
                  c("Replication-associated", "Accessory", "Structural", "Entry"))
  
  effects_df <- lr_sum %>%
    as_tibble() %>%
    select(Estimate, std_error = `Std. Error`) %>%
    add_column(gene = rownames(lr_sum)) %>%
    mutate(gene = gsub("interaction\\(gene_type\\, is_jump\\)", "",
                       gene)) %>%
    separate(gene, c("gene_type", "is_jump"), "\\.") %>%
    mutate(is_jump = ifelse(gene_type == "(Intercept)", F, is_jump)) %>%
    mutate(gene_type = ifelse(gene_type == "(Intercept)", "Structural", gene_type))
  
  ycoord <- (min(effects_df$Estimate) + max(effects_df$Estimate)) / 2
  
  effects_plt <- effects_df %>%
    ggplot(aes(x = is_jump, y = Estimate, fill = gene_type)) +
    facet_grid(.~factor(gene_type, unique(parsed_models$gene_type))) +
    scale_fill_manual(values = pal) +
    geom_point(size = 3,
               pch = 21,
               color = "black") +
    geom_errorbar(aes(ymin = Estimate - std_error, ymax = Estimate + std_error),
                  width = 0) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Is host jump?", y = "Parameter estimate") +
    geom_text(aes(x = 1.5, 
                  y = ycoord, 
                  label = str_glue("Effect={effect}\np={p_val}")),
              size = 3,
              data = parsed_models) +
    labs(title = fam_name)
  # ylim(-0.5, 1.5)
  effects_plt
  
  ggsave(str_glue("results/dnds_out/family_plots/{fam_name}.effects.pdf"), 
         plot = effects_plt,
         dpi = 600,
         width = 5, height = 3)
  
  dat_plt <- dat %>%
    mutate(gene_type = factor(gene_type, parsed_models$gene_type)) %>%
    ggplot(aes(x = is_jump, y = corr_effect, fill = gene_type)) +
    facet_grid(.~factor(gene_type, unique(parsed_models$gene_type))) +
    geom_violin() +
    geom_boxplot(width = 0.1,
                 outlier.shape = NA,
                 alpha = 0.8) +
    theme_classic() +
    scale_fill_manual(values = pal) +
    theme(legend.position = "none") + 
    geom_text(aes(x = 1.5, y = 2.5, 
                  label = str_glue("Effect={effect}\np={p_val}")),
              size = 3,
              data = parsed_models) +
    labs(title = fam_name)
  # ylim(-2, 3)
  
  resid_plt <- tibble(x = lr$residuals) %>%
    ggplot(aes(x)) +
    geom_density() +
    labs(x = "Log10(dn/ds)", y = "Density") +
    theme_classic()
  
  ggsave(str_glue("results/dnds_out/family_plots/{fam_name}.raw.pdf"), 
         plot = dat_plt,
         dpi = 600,
         width = 5, height = 3)
  
  ggsave(str_glue("results/dnds_out/family_plots/{fam_name}.residuals.pdf"), 
         plot = resid_plt,
         dpi = 600,
         width = 5, height = 3)
}

# jump_df <- fread("results/dnds_out/all_jumps.dnds.diff_hosts.genus_counts.csv")
# jump_df %>%
#   separate(clique_name, c("family"), "_") %>%
#   left_join(genome_type) %>%
#   filter(is_jump) %>%
#   group_by(family, genome_type) %>%
#   summarise(n = n_distinct(tip_name)) %>%
#   arrange(desc(n))
# meta  %>%
#   left_join(genome_type) %>%
#   filter(!is_segmented) %>%
#   group_by(genome_type) %>%
#   summarise(median_length = median(genome_length))
