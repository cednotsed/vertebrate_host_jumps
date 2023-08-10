rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)

# Permutation test function
permute_coeff <- function(best_model_formula, coef_idx, dat) {
  set.seed(66)
  perm_morsels <- foreach(i = seq(1000)) %do% {
    temp_df <- dat %>%
      ungroup() %>%
      mutate(n_hosts = sample(dat$n_hosts, nrow(dat), replace = F))
    
    temp_mod <- glm(as.formula(best_model_formula),
                    temp_df,
                    family = poisson())
    
    effect <- coef(temp_mod)[[coef_idx]]
    
    tibble(coefficient = effect)
  }
  
  return(bind_rows(perm_morsels))
}

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv"))

genome_type <- fread("data/metadata/genome_type_metadata.csv")

dnds_df <- fread("results/dnds_out/dnds_out.diff_hosts.genus_counts.csv") %>%
  dplyr::rename(clique_name = cluster) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type)

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

res_filt <- dnds_df %>% 
  # filter(!is.na(ka)) %>%
  mutate(ks = ifelse(ks == 0 | ks < 0, 0.000001, ks),
         ka = ifelse(ka == 0, 0.000001, ka)) %>%
  mutate(ka_minus_ks = ka - ks,
         kaks = ka / ks) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(genome_length) %>%
  filter(cds_length / median_length > 0.5)

clique_counts <- res_filt %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup() %>%
  filter(n_cliques >= 5)

plot_df <- res_filt %>% 
  group_by(clique_name, is_jump, n_hosts) %>%
  summarise(min_kaks = min(kaks),
            median_cds = median(cds_length)) %>%
  left_join(clique_counts) %>%
  filter(n_cliques >= 5)

plot_df %>%
  ggplot(aes(x = factor(n_hosts), y = log10(min_kaks))) +
  geom_boxplot(aes(fill = is_jump)) +
  theme_classic() +
  labs(x = "No. host genera", y = "log10(min. dn/ds)",
       fill = "Is host jump?") +
  geom_text(aes(x = n_hosts, y = -4.5, label = str_glue("n={n_cliques}")),
            data = clique_counts,
            color = "grey39",
            size = 3)

ggsave("results/dnds_out/dnds_visualisation.png", width = 5, height = 5)


lr1_jump <- glm(n_hosts ~ median_cds + min_kaks,
                data = plot_df %>% filter(is_jump),
                family = poisson())

lr2_jump <- glm(n_hosts ~ median_cds + log(min_kaks),
                data = plot_df %>% filter(is_jump),
                family = poisson())

lr3_jump <- glm(n_hosts ~ log(median_cds) + log(min_kaks),
                data = plot_df %>% filter(is_jump),
                family = poisson())

lr4_jump <- glm(n_hosts ~ log(median_cds) + min_kaks,
                data = plot_df %>% filter(is_jump),
                family = poisson())

comp <- compareGLM(lr1_jump, lr2_jump, lr3_jump, lr4_jump)

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

best_model <- get(str_glue("lr{best_model_idx}_jump"))
obs_effect_jump <- coef(best_model)[[2]]

print(summary(best_model))

# Test non jump
lr1_non <- glm(n_hosts ~ min_kaks,
               data = plot_df %>% filter(!is_jump),
               family = poisson())

lr2_non <- glm(n_hosts ~ log10(min_kaks),
               data = plot_df %>% filter(!is_jump),
               family = poisson())

lr3_non <- glm(n_hosts ~ log(min_kaks),
               data = plot_df %>% filter(!is_jump),
               family = poisson())

comp <- compareGLM(lr1_non, lr2_non, lr3_non)

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

best_model <- get(str_glue("lr{best_model_idx}_non"))
obs_effect_nonjump <- coef(best_model)[[2]]

print(summary(best_model))

obs_slope_diff <- obs_effect_nonjump - obs_effect_jump
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
