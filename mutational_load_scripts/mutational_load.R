rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)

# Permutation test function
permute_coeff <- function(best_model_formula, coef_idx) {
  set.seed(66)
  perm_morsels <- foreach(i = seq(1000)) %do% {
    temp_df <- merged_df %>%
      ungroup() %>%
      mutate(n_hosts = sample(merged_df$n_hosts, nrow(merged_df), replace = F))
    
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
  # filter(prop_traversed < 0.2)

# Get distinct host jumps with lowest mutation load
# cl <- makeCluster(12)
# registerDoParallel(cl)

jump_morsels <- foreach(cluster_name = unique(res_df$clique_name),
                        .packages = c("tidyverse",
                                      "foreach")) %do% {
                                        print(cluster_name)
                                        res_filt <- res_df %>%
                                          filter(clique_name == cluster_name)
                                        
                                        host_jumps <- res_filt %>%
                                          distinct(anc_state, tip_state)
                                        
                                        crumbs <- foreach(i = seq(nrow(host_jumps))) %do% {
                                          anc <- host_jumps[i, ]$anc_state
                                          tip <- host_jumps[i, ]$tip_state
                                          
                                          res_filt %>%
                                            filter(anc_state == anc, tip_state == tip) %>%
                                            arrange(patristic_dist) %>%
                                            head(1)
                                        }
                                        
                                        return(bind_rows(crumbs))
                                      }

# stopCluster(cl)

res <- bind_rows(jump_morsels) %>%
  separate(clique_name, c("viral_family"), "_", remove = F)

host_counts <- meta %>%
  dplyr::rename(clique_name = cluster) %>%
  group_by(clique_name) %>%
  summarise(n_hosts = n_distinct(host_order),
            n_genomes = n_distinct(accession)) %>%
  arrange(desc(n_hosts))

merged_df <- bind_rows(jump_morsels) %>%
  separate(clique_name, c("viral_family"), "_", remove = F) %>%
  left_join(host_counts) %>%
  group_by(clique_name, viral_family, n_genomes, n_hosts) %>%
  summarise(min_dist = min(patristic_dist),
            min_nodes = min(n_traverses))

# Visualise effects
merged_df %>%
  ggplot(aes(y = n_hosts, x = log10(n_genomes))) +
  geom_smooth(method = "lm",
              color = "black",
              fill = "tan1") +
  geom_point(color = "steelblue4", alpha = 0.5) +
  theme_classic() +
  labs(x = "log10(No. of viral genomes)", y = "No. of hosts")

ggsave("results/mutational_load_out/hosts_vs_genomes.png", width = 5, height = 3)

# Hosts vs mutation
corr_test <- cor.test(merged_df$n_hosts, log10(merged_df$min_dist))
r_val <- signif(corr_test$estimate, 3)
p_val <- signif(corr_test$p.value, 3)

merged_df %>%
  ggplot(aes(y = n_hosts, x = log10(min_dist))) +
  geom_smooth(color = "black",
              method = "lm",
              fill = "tan1") +
  geom_point(color = "steelblue4", alpha = 0.5) +
  theme_classic() +
  # xlim(0, 3.5) +
  labs(x = "log10(Min. patristic distance)", y = "No. of hosts",
       title = str_glue("Pearson's R = {r_val}, p = {p_val}"))

ggsave("results/mutational_load_out/hosts_vs_mutation.png", width = 5, height = 3)

corr_test <- cor.test(log10(merged_df$n_genomes), merged_df$min_dist)
r_val <- signif(corr_test$estimate, 3)
p_val <- signif(corr_test$p.value, 3)
p_val
merged_df %>%
  ggplot(aes(x = log10(n_genomes), y = min_dist)) +
  geom_smooth(color = "black",
              method = "lm",
              fill = "tan1") +
  geom_point(color = "steelblue4", alpha = 0.5) +
  theme_classic() +
  # xlim(0, 3.5) +
  labs(x = "Log10(No. genomes)", y = "Min. patristic distance",
       title = str_glue("Pearson's R = {r_val}, p = {p_val}"))

ggsave("results/mutational_load_out/mutation_vs_sampling.png", width = 5, height = 3)

# Without genome correction
lr1 <- glm(n_hosts ~ min_dist,
           merged_df,
           family = poisson())

lr2 <- glm(n_hosts ~ log10(min_dist),
           merged_df,
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
perm_no_correction <- permute_coeff(best_model_formula, 2)
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

ggsave("results/mutational_load_out/perm_test_slope.no_correction.png", width = 5, height = 3)

# With genome correction
lr1 <- glm(n_hosts ~ n_genomes + min_dist,
           merged_df,
           family = poisson())

lr2 <- glm(n_hosts ~ log10(n_genomes) + min_dist,
           merged_df,
           family = poisson())

lr3 <- glm(n_hosts ~ n_genomes + log10(min_dist),
           merged_df,
           family = poisson())

lr4 <- glm(n_hosts ~ log10(n_genomes) + log10(min_dist),
           merged_df,
           family = poisson())

comp <- compareGLM(lr1, lr2, lr3, lr4)

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
summary(best_model)
obs_effect <- coef(best_model)[[3]]

# Permutate uncorrected model
perm_no_correction <- permute_coeff(best_model_formula, 3)
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

ggsave("results/mutational_load_out/perm_test_slope.corrected.png", width = 5, height = 3)

# 
# merged_df2 <- res_df %>%
#   group_by(clique_name) %>%
#   summarise(median_dist = median(patristic_dist)) %>%
#   left_join(host_counts)
# 
# lr1 <- glm(n_hosts ~ n_genomes + median_dist,
#            merged_df2,
#            family = poisson())
# 
# lr2 <- glm(n_hosts ~ log10(n_genomes) + median_dist,
#            merged_df2,
#            family = poisson())
# 
# lr3 <- glm(n_hosts ~ n_genomes + log10(median_dist),
#            merged_df2,
#            family = poisson())
# 
# lr4 <- glm(n_hosts ~ log10(n_genomes) + log10(median_dist),
#            merged_df2,
#            family = poisson())
# lr5 <- glm(n_hosts ~  median_dist + log10(n_genomes),
#            merged_df2,
#            family = poisson())
# 
# lr6 <- glm(n_hosts ~  log10(median_dist),
#            merged_df2,
#            family = poisson())
# lr7 <- glm(n_hosts ~  median_dist,
#            merged_df2,
#            family = poisson())
# 
# comp_res <- compareGLM(lr1, lr2, lr3, lr4, lr5, lr6, lr7)
# comp_res$Fit.criteria %>%
#   rownames_to_column("model") %>%
#   arrange(AIC) %>%
#   head(1)
# 
# summary(lr2)
# 
# # anova(lr8)
# # merged_df
# # 
# # corr_lm <- lm(n_hosts ~ log10(n_genomes),
# #               data = merged_df2)
# # 
# # corr_df <- merged_df2 %>%
# #   mutate(corr_range = corr_lm$residuals)
# # test <- lm(patristic_dist ~ corr_range,
# #            corr_df)
# # summary(test)
# # 
# # corr_df %>%
# #   ggplot(aes(x = patristic_dist, y = corr_range)) +
# #   geom_point()
# # 
# # require(mgcv)
# # require(gratia)
# # gmm <- gam(n_hosts ~ s(log10(n_genomes)) + s(log10(min_dist), k = 3),
# #            data = merged_df,
# #            select = T,
# #            family = poisson())
# # summary(gmm)
# # draw(gmm, residuals = T)
