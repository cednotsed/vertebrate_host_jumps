rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)

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


dnds_list <- list.files("results/dnds_out/all_jumps.temp_results/", full.names = T)
dnds_df <- foreach(file_name = dnds_list, .combine = c("bind_rows")) %do% {
  fread(file_name)
}

# dnds_df <- dnds_df %>% mutate(kaks = ka / ks)
# fwrite(dnds_df, "results/dnds_out/all_jumps.dnds.diff_hosts.genus_counts.csv")

# Remove zeroes
dnds_filt <-  dnds_df %>%
  mutate(ks = ifelse(ks == 0 | ks < 0, 0, ks),
                   ka = ifelse(ka == 0, 0, ka)) %>%
  filter(!(ka == 0 | ks == 0)) %>%
  filter(ka != 0 & ks != 0) %>%
  mutate(kaks = ka / ks)

clique_counts <- dnds_filt %>%
  left_join(host_counts) %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

# Get min. kaks per clique
iter_df <- dnds_filt %>%
  distinct(clique_name, is_jump)

min_df <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  row <- iter_df[i, ]
  dnds_filt %>%
    filter(clique_name == row$clique_name,
           is_jump == row$is_jump) %>%
    arrange(kaks) %>%
    head(1)
}

# Merge data
merged_df <- min_df %>%
  select(clique_name, is_jump, patristic_dist, min_kaks = kaks) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

# merged_df <- dnds_filt %>%
#   select(clique_name, is_jump, kaks) %>%
#   group_by(clique_name, is_jump) %>%
#   summarise(min_kaks = min(kaks)) %>%
#   separate(clique_name, c("family"), "_", remove = F) %>%
#   left_join(genome_type) %>%
#   left_join(genome_counts) %>%
#   left_join(host_counts) %>%
#   left_join(clique_counts)

merged_df %>%
  fwrite("results/dnds_out/parsed_min_dnds.csv")

# Jump/non-jump counts
type_counts <- merged_df %>%
  group_by(clique_name) %>%
  summarise(n_types = n_distinct(is_jump)) %>%
  filter(n_types == 2)

merged_filt <- merged_df %>%
  filter(clique_name %in% type_counts$clique_name)

# Visualise overall dnds
wilcox <- wilcox.test(log10(merged_filt$min_kaks) ~ merged_filt$is_jump,
                      alternative = "less")

merged_filt %>%
  ggplot(aes(x = is_jump, y = log10(min_kaks), fill = is_jump)) +
  geom_boxplot(position = position_nudge(x = 0.05, y = 0),
               width = 0.1,
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.1, y = 0), alpha = 1) +
  theme_classic() +
  coord_flip() +
  geom_pwc() +
  labs(x = "Is host jump?", y = "log10(min. dn/ds)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("steelblue3", "indianred3"))

ggsave("results/dnds_out/overall_jump_v_non_jump_dnds.pdf", width = 5, height = 3.5)

# Correcting for family and n_genomes
logreg <- glm(is_jump ~ log(n_genomes) + family + log(min_kaks),
              data = merged_filt,
              family = "binomial")

summary(logreg)
exp(coef(logreg)[["log(min_kaks)"]])
nrow(merged_filt) - nrow(coef(summary(logreg)))

pnorm(0.87062 / 0.17973, lower.tail = F)
# Visualise by host range

# merged_df %>%
#   ggplot(aes(x = n_hosts, y = log10(min_kaks), fill = is_jump, color = is_jump)) +
#   geom_point() +
#   theme_classic() +
#   labs(x = "No. host genera", y = "log10(min. dn/ds)",
#        color = "Is host jump?") +
#   ylim(-2.4, 0.1) +
#   scale_fill_manual(values = c("steelblue3", "indianred3")) +
  # facet_grid(rows = vars(is_jump)) +
  # geom_smooth(method = "lm")

merged_filt %>%
  ggplot(aes(x = factor(n_hosts), y = log10(min_kaks), fill = is_jump)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  labs(x = "No. host genera", y = "log10(min. dn/ds)",
       color = "Is host jump?") +
  ylim(-2.4, 0.1) +
  scale_fill_manual(values = c("steelblue3", "indianred3"))

ggsave("results/dnds_out/dnds_v_host_range.pdf", width = 6, height = 3.5)

merged_df %>%
  ggplot(aes(x = n_hosts, y = log10(min_kaks), fill = is_jump, color = is_jump)) +
  geom_point(color = "black",
             pch = 21) +
  theme_classic() +
  labs(x = "No. host genera", y = "log10(min. dn/ds)",
       color = "Is host jump?") +
  # ylim(-2.4, 0.1) +
  scale_fill_manual(values = c("steelblue3", "indianred3")) +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  geom_smooth(method = "lm")

ggsave("results/dnds_out/dnds_v_host_range.scatter.pdf", width = 6, height = 3.5)

# Test for difference in slopes
formula_string <- "n_hosts ~ log(n_genomes) + family + log(min_kaks)"
jump_mod <- glm(as.formula(formula_string),
                data = merged_df %>% filter(is_jump),
                family = poisson())

nonjump_mod <- glm(as.formula(formula_string),
                   data = merged_df %>% filter(!is_jump),
                   family = poisson())

nrow(merged_df %>% filter(is_jump)) - nrow(coef(summary(nonjump_mod)))
nrow(merged_df %>% filter(!is_jump)) - nrow(coef(summary(nonjump_mod)))

sink("results/dnds_out/jump_poisson_model.txt")
summary(jump_mod)
sink()

sink("results/dnds_out/nonjump_poisson_model.txt")
summary(nonjump_mod)
sink()


obs_slope_diff <- coef(jump_mod)[["log(min_kaks)"]] - coef(nonjump_mod)[["log(min_kaks)"]]
obs_slope_diff

# Permutation test for difference in slopes
set.seed(66)

perm_morsels <- foreach(i = seq(1000), .combine = "c") %do% {
  perm_temp <- merged_df %>%
    group_by(clique_name) %>%
    mutate(is_jump = sample(is_jump, replace = F))
  
  lr1_jump <- glm(as.formula(formula_string),
                  data = perm_temp %>% filter(is_jump),
                  family = poisson())
  
  lr1_non <- glm(as.formula(formula_string),
                 data = perm_temp %>% filter(!is_jump),
                 family = poisson())
  
  slope_diff <- coef(lr1_jump)[["log(min_kaks)"]] - coef(lr1_non)[["log(min_kaks)"]]
  
  return(slope_diff)
}

p_val <- tibble(slope_diff = perm_morsels) %>%
  summarise(p = sum(slope_diff <= obs_slope_diff) / n())

tibble(slope_diff = perm_morsels) %>%
  ggplot(aes(x = slope_diff)) +
  geom_histogram(fill = "steelblue3") +
  geom_vline(xintercept = obs_slope_diff,
             lty = "dashed",
             color = "indianred") +
  theme_classic() +
  labs(x = "Diff. in slope estimate", y = "No. permutations",
       title = str_glue("Permutation test: p = {p_val}"))

ggsave("results/dnds_out/permutation.diff_hosts.genus_counts.png", width = 8, height = 5)

# lr <- lm(log(min_kaks) ~ is_jump + log(patristic_dist), 
#    data = merged_df)
# summary(lr)
# hist(lr$residuals)
# 
# merged_df %>%
#   ggplot(aes(log(patristic_dist), log(min_kaks), color = is_jump)) +
#   geom_point()
# summary(glm(n_hosts ~ log10(min_kaks), 
#            data = plot_df %>% filter(is_jump),
#            family = poisson()))
# 
# lr1_jump <- glm(n_hosts ~ min_kaks,
#                 data = plot_df %>% filter(is_jump),
#                 family = poisson())
# 
# lr2_jump <- glm(n_hosts ~ log(min_kaks),
#                 data = plot_df %>% filter(is_jump),
#                 family = poisson())
# 
# lr3_jump <- glm(n_hosts ~ log10(min_kaks),
#                 data = plot_df %>% filter(is_jump),
#                 family = poisson())
# 
# comp <- compareGLM(lr1_jump, lr2_jump, lr3_jump)
# 
# best_model_idx <- deframe(comp$Fit.criteria %>%
#                             rownames_to_column("model") %>%
#                             arrange(AIC) %>%
#                             select(model) %>%
#                             head(1))
# 
# best_model_formula <- deframe(comp$Models %>%
#                                 as_tibble() %>%
#                                 rownames_to_column("model") %>%
#                                 filter(model == best_model_idx) %>%
#                                 select(Formula))
# 
# best_model <- get(str_glue("lr{best_model_idx}_jump"))
# obs_effect_jump <- coef(best_model)[[2]]
# 
# print(summary(best_model))
# 
# # Test non jump
# lr1_non <- glm(n_hosts ~ min_kaks,
#                data = plot_df %>% filter(!is_jump),
#                family = poisson())
# 
# lr2_non <- glm(n_hosts ~ log10(min_kaks),
#                data = plot_df %>% filter(!is_jump),
#                family = poisson())
# 
# lr3_non <- glm(n_hosts ~ log(min_kaks),
#                data = plot_df %>% filter(!is_jump),
#                family = poisson())
# 
# comp <- compareGLM(lr1_non, lr2_non, lr3_non)
# 
# best_model_idx <- deframe(comp$Fit.criteria %>%
#                             rownames_to_column("model") %>%
#                             arrange(AIC) %>%
#                             select(model) %>%
#                             head(1))
# 
# best_model_formula <- deframe(comp$Models %>%
#                                 as_tibble() %>%
#                                 rownames_to_column("model") %>%
#                                 filter(model == best_model_idx) %>%
#                                 select(Formula))
# 
# best_model <- get(str_glue("lr{best_model_idx}_non"))
# obs_effect_nonjump <- coef(best_model)[[2]]
# 
# print(summary(best_model))
# 
# obs_slope_diff <- obs_effect_nonjump - obs_effect_jump
# host_counts <- meta %>%
#   dplyr::rename(clique_name = cluster) %>%
#   group_by(clique_name) %>%
#   summarise(n_hosts = n_distinct(host_genus),
#             n_genomes = n_distinct(accession)) %>%
#   arrange(desc(n_hosts))
# 
# jump <- res_filt %>%
#   filter(is_jump) %>%
#   select(clique_name, anc_name, tip_name, anc_state, tip_state, jump_kaks = kaks)
#   
# non_jump <- res_filt %>%
#   filter(!is_jump) %>%
#   select(clique_name, anc_state, non_jump_kaks = kaks) %>%
#   distinct()
# 
# merged_df <- jump %>%
#   left_join(non_jump) %>%
#   mutate(ratio = jump_kaks / non_jump_kaks) %>%
#   left_join(host_counts) 
# 
# res_filt %>%
#   # ggplot(aes(ka_minus_ks,  kaks)) +
#   # geom_point()
#   left_join(host_counts) %>%
#   # separate(clique_name, into = c("viral_family", "_")) %>%
#   # group_by(clique_name, is_jump) %>%
#   # summarise(median_range = max(n_hosts), 
#   #           median_kaks = median(ka / ks),
#   #           var_kaks = var(ka - ks)) %>%
#   ggplot(aes(x = median_range, y = log10(median_kaks), color = is_jump)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# merged_df %>%
#   ggplot(aes(x = n_hosts, y = log10(ratio))) +
#   geom_point() +
#   geom_smooth()
#   geom
# # Hosts vs ka - ks
# corr_test <- cor.test(merged_df$n_hosts, merged_df$ka_minus_ks)
# r_val <- signif(corr_test$estimate, 3)
# p_val <- signif(corr_test$p.value, 3)
# 
# merged %>%
#   left_join(host_counts) %>%
#   arrange(n_hosts) %>%
#   mutate(host_range = case_when(n_hosts <= 2 ~ "<3",
#                                 n_hosts > 2 & n_hosts <= 3 ~ "3-5",
#                                 n_hosts > 3 ~ ">5")) %>%
#   mutate(host_range = factor(host_range, c("<3", "3-5", ">5"))) %>%
#   # group_by(host_range) %>%
#   # summarise(n_distinct(clique_name))
#   ggplot(aes(x = host_range, y = log10(kaks), color = is_jump)) +
#   see::geom_violinhalf() +
#   # geom_smooth(method = "lm",
#   #             fill = "tan1") +
#   # geom_point(alpha = 0.5) +
#   geom_boxplot() +
#   theme_classic() +
#   # xlim(0, 3.5) +
#   labs(x = "No. host orders", y = "log10(Ks/Ka)"
#        # title = str_glue("Pearson's R = {r_val}, p = {p_val}")
#        )
# 
# # merged_filt %>%
# #   ggplot(aes(ka_minus_ks, kaks)) +
# #   geom_point()
# # ggsave("results/dnds_out/hosts_vs_ka_minus_ks.png", width = 5, height = 3)
# 
# # Hosts vs ka/ks
# merged_filt <- merged_df %>%
#   filter(!is.infinite(kaks))
# corr_test <- cor.test(merged_filt$n_hosts, merged_filt$kaks)
# r_val <- signif(corr_test$estimate, 3)
# p_val <- signif(corr_test$p.value, 3)
# 
# merged_filt %>%
#   ggplot(aes(x = n_hosts, y = log10(kaks), color = is_jump)) +
#   # geom_smooth(color = "black",
#   #             # method = "lm",
#   #             fill = "tan1") +
#   geom_point(alpha = 0.5) +
#   theme_classic() +
#   # xlim(0, 3.5) +
#   labs(x = "No. host orders", y = "ka / ks",
#        title = str_glue("Pearson's R = {r_val}, p = {p_val}"))
# 
# ggsave("results/dnds_out/hosts_vs_kaks.png", width = 5, height = 3)
# 
# # corr_test <- cor.test(log10(merged_filt$n_genomes), merged_filt$kaks)
# # r_val <- signif(corr_test$estimate, 3)
# # p_val <- signif(corr_test$p.value, 3)
# # 
# # merged_filt %>%
# #   ggplot(aes(x = log10(n_genomes), y = kaks)) +
# #   geom_smooth(color = "black",
# #               method = "lm",
# #               fill = "tan1") +
# #   geom_point(color = "steelblue4", alpha = 0.5) +
# #   theme_classic() +
# #   # xlim(0, 3.5) +
# #   labs(x = "Log10(No. genomes)", y = "Min. patristic distance",
# #        title = str_glue("Pearson's R = {r_val}, p = {p_val}"))
# # 
# # ggsave("results/mutational_load_out/kaks_vs_sampling.png", width = 5, height = 3)
# 
# # Without genome correction
# lr1 <- glm(n_hosts ~ kaks,
#            merged_filt,
#            family = poisson())
# 
# lr2 <- glm(n_hosts ~ log10(kaks),
#            merged_filt,
#            family = poisson())
# 
# comp <- compareGLM(lr1, lr2)
# best_model_idx <- deframe(comp$Fit.criteria %>%
#                             rownames_to_column("model") %>%
#                             arrange(AIC) %>%
#                             select(model) %>%
#                             head(1))
# 
# best_model_formula <- deframe(comp$Models %>%
#                                 as_tibble() %>%
#                                 rownames_to_column("model") %>%
#                                 filter(model == best_model_idx) %>%
#                                 select(Formula))
# 
# best_model <- get(paste0("lr", best_model_idx))
# obs_effect <- coef(best_model)[[2]]
# 
# summary(best_model)
# 
# # Permutate uncorrected model
# perm_no_correction <- permute_coeff(best_model_formula, 2, merged_filt)
# p_val <- perm_no_correction %>%
#   summarise(p_value = sum(coefficient <= obs_effect) / n())
# 
# perm_no_correction %>%
#   ggplot(aes(x = coefficient)) +
#   geom_histogram(fill = "steelblue4",
#                  bins = 50) +
#   geom_vline(xintercept = obs_effect,
#              lty = "dashed",
#              color = "red") +
#   theme_classic() +
#   labs(x = "Slope estimate", y = "No. permutations",
#        title = str_glue("Perm. test p-value = {p_val}"))
# 
# ggsave("results/dnds_out/perm_test_slope.no_correction.png", width = 5, height = 3)
# 
