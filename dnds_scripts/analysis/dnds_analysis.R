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

res_df <- fread("results/source_sink_analysis/putative_host_jumps.csv") %>%
  mutate(prop_traversed = n_traverses / total_depth) 

dnds_df <- fread("results/dnds_out/dnds_out.csv")

res_filt <- res_df %>%
  left_join(dnds_df) %>% 
  filter(!is.na(ka)) %>%
  mutate(ks = ifelse(is.na(ks), 0, ks)) %>%
  mutate(ka_minus_ks = ka - ks,
         kaks = ka / ks)

host_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_order),
            n_genomes = n_distinct(accession)) %>%
  arrange(desc(n_hosts))

merged_df <- res_filt %>%
  left_join(host_counts)

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
ggsave("results/dnds_out/hosts_vs_ka_minus_ks.png", width = 5, height = 3)

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

ggsave("results/dnds_out/hosts_vs_kaks.png", width = 5, height = 3)

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

