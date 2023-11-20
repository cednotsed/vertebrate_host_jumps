rm(list = ls())
setwd("c:/git_repos/viral_sharing/")
require(tidyverse)
require(data.table)
require(Hmisc)
require(mgcv)
require(gratia)
require(ggpubr)
require(foreach)
require(doParallel)

# SRA sequencing effort
meta <- bind_rows(fread("data/metadata/sra_metadata/final_animal_metagenomes_metadata.csv"),
                  fread("data/metadata/sra_metadata/final_human_metagenomes_metadata.csv")) %>%
  mutate(host_genus = capitalize(host_genus)) %>%
  mutate(host_class = capitalize(host_class))

effort_df <- meta %>%
  group_by(host_genus) %>%
  summarise(effort = log10(sum(bases)),
            n_libraries = log10(n_distinct(run)))

# Virion data
virion <- fread("data/VIRION.v0.2.1_240922/Virion.csv.gz") %>%
  filter(DetectionMethod %in% c("PCR/Sequencing", "Isolation/Observation")) %>%
  mutate(host_genus = capitalize(HostGenus),
         host_family = capitalize(HostFamily),
         host_order = capitalize(HostOrder),
         host_class = capitalize(HostClass)) %>%
  filter(host_class == "Mammalia")

host_tax <- fread("data/VIRION.v0.2.1_240922/TaxonomyHost.csv")
virus_tax <- fread("data/VIRION.v0.2.1_240922/TaxonomyVirus.csv")

diversity_df <- virion %>%
  group_by(host_genus, host_class, host_order, host_family) %>%
  summarise(n_species = n_distinct(Virus)) %>%
  arrange(desc(n_species)) %>%
  filter(!(host_genus %in% c("Homo", ""))) %>%
  left_join(effort_df) %>%
  mutate(n_libraries = ifelse(is.na(n_libraries), 0, n_libraries),
         effort = ifelse(is.na(effort), 0, effort)) %>%
  ungroup()

# Correct for clique representation and sequencing effort
# Visualise confounders
diversity_df %>%
  select(-host_genus, -host_class, -host_order, -host_family) %>%
  pivot_longer(!n_species, names_to = "confounders", values_to = "value") %>%
  mutate(confounders = case_when(confounders == "effort" ~ "log10(sum bases sequenced)",
                                 confounders == "n_libraries" ~ "log10(SRA libraries)")) %>%
  ggplot(aes(x = value, y = log10(n_species))) +
  facet_grid(cols = vars(confounders),
             scale = "free") +
  geom_point() +
  geom_smooth(method = "lm")

ggsave("results/modelling_out/viral_diversity/virion_data/confounders_linear_visualisation.png", width = 6, height = 4)

# Get residuals
diversity_lm <- lm(log10(n_species) ~ n_libraries,
                   data = diversity_df)

corr_diversity_df <- diversity_df %>%
  mutate(corr_diversity = diversity_lm$residuals)

traits_df <- fread("data/combine_data_soria_2021/trait_data_imputed.csv",
                   stringsAsFactors = T) %>%
  rename_all(~tolower(gsub(" |-", "_", .x))) %>%
  rename_all(~gsub("-_", "", .x)) %>%
  dplyr::rename(host_genus = genus) %>%
  filter(host_genus %in% corr_diversity_df$host_genus)

# Filter traits
prop_na <- colSums(is.na(traits_df)) / nrow(traits_df)
to_keep <- names(prop_na)[prop_na < 0.20]
to_keep <- to_keep[grepl("_d|_g|_n|_mm|km|dphy", to_keep)]
to_keep <- to_keep[!(to_keep %in% c("habitat_breadth_n", 
                                    "terrestrial_non_volant",
                                    "det_diet_breadth_n",
                                    "det_nect"))]

to_keep
traits_filt <- traits_df %>%
  mutate(host_genus = as.character(host_genus)) %>%
  select(all_of(to_keep)) %>%
  group_by(host_genus) %>%
  summarise_if(is.numeric, median, na.rm = T) %>%
  filter(host_genus %in% corr_diversity_df$host_genus)

apply(traits_filt %>% select(-host_genus), 2, function(x) {sum(is.na(x))})
# Run PCA
traits_X <- traits_filt %>% 
  select(-host_genus)

pca <- prcomp(traits_X, retx = T, scale. = T, center = T)

pca_df <- as.data.frame(pca$x) %>%
  mutate(host_genus = traits_filt$host_genus)

var_explained <- pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100
plot(var_explained)

merged_df <- corr_diversity_df %>%
  inner_join(pca_df)

# filter(!(host_genus %in% c("Homo", "Sus", "Equus", "Mus", "Felis", "Canis")))

predictors <- c("PC1", "PC2", "PC3", "PC4", "PC5")

f_string <- "corr_diversity ~"
for (i in seq(length(predictors))) {
  predictor <- predictors[i]
  if (i == 1) {
    f_string <- str_glue("{f_string}s({predictor})") 
  } else {
    f_string <- str_glue("{f_string} + s({predictor})") 
  }
}

f_string
gmm <- gam(as.formula(f_string), 
           data = merged_df)

res <- summary(gmm)
obs_dev_explained <- res$dev.expl
res

# draw(gmm, residuals = T)

# Permutation test
perm_iters <- 1000

set.seed(66)
perm_morsels <- foreach(i = seq(perm_iters), 
                        .combine = "c") %do% {
                          temp_dat <- merged_df %>%
                            mutate(corr_diversity = sample(merged_df$corr_diversity, replace = F))
                          
                          gam_temp <- gam(as.formula(f_string), 
                                          data = temp_dat)
                          
                          dev_explained <- summary(gam_temp)$dev.expl
                          
                          return(dev_explained)
                        }

p <- sum(perm_morsels > obs_dev_explained) / length(perm_morsels)

tibble(x = perm_morsels) %>%
  ggplot(aes(x = x)) +
  geom_histogram(fill = "cadetblue") +
  geom_vline(xintercept = obs_dev_explained,
             lty = "dashed",
             col = "red") +
  geom_text(aes(x = obs_dev_explained + 0.05, 
                y = 75,
                label = str_glue("Observed = {round(obs_dev_explained, 2)}")),
            color = "red") + 
  theme_classic2() +
  labs(x = "Deviance explained", y = "No. of permutations", title = str_glue("Permutation test p = {p}, n_genera = {n_distinct(merged_df$host_genus)}"))

ggsave("results/modelling_out/viral_diversity/virion_data/permutation_test.PCA.png", width = 8, height = 5)


# # merged_df %>%
# #   pivot_longer(!c(host_genus, n_cliques, log_sum_bases, n_cliques_rep, corr_diversity), 
# #                names_to = "variables", values_to = "value") %>%
# #   ggplot(aes(x = corr_diversity, y = value)) +
# #   facet_grid(rows = vars(variables), scales = "free") +
# #   geom_point() +
# #   geom_smooth() 
# #   ggplot(aes(log10(corr_diversity)))
# # Plot all partial plots
# draw(gmm, residuals = T)
# 
# # Extract significant predictors
# gmm_sum <- summary(gmm)
# smooth_table <- as.data.frame(gmm_sum$s.table) %>%
#   rename_all(~tolower(gsub("-", "_", .x))) %>%
#   rownames_to_column("smooth") %>%
#   mutate(smooth = gsub("s\\(|\\)", "", smooth)) %>%
#   filter(p_value < 0.05)
# 
# # Plot partial residual plots for significant smooth terms
# to_plot <- unique(smooth_table$smooth)
# 
# plot_partial <- function(i) {
#   term <- to_plot[i]
#   print(term)
#   
#   # Get partial residuals
#   partial_resids <- deframe(partial_residuals(gmm, select = str_glue("s({term})")))
#   resid_df <- tibble(resid = partial_resids, pred_value = gmm$model[, term])
#   
#   # Get smooth estimates
#   smm <- smooth_estimates(gmm, smooth = str_glue("s({term})"))
#   
#   # Get test stats
#   f_val <- round(smooth_table[smooth_table$smooth == term, "f"], 2)
#   p_val <- signif(smooth_table[smooth_table$smooth == term, "p_value"], 3)
#   
#   plt <- smm %>%
#     add_confint() %>%
#     ggplot(aes(y = est, x = !!sym(term))) +
#     geom_point(aes(x = pred_value, y = partial_resids), 
#                color = "steelblue4",
#                alpha = 0.5,
#                size = 3,
#                data = resid_df) +
#     geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
#                 alpha = 0.2, 
#                 fill = "orange") +
#     geom_line(colour = "black", 
#               linewidth = 1,
#               alpha = 0.8) +
#     theme_bw() +
#     labs(y = "Effect on corrected viral diversity",
#          x = term,
#          title = str_glue("F = {f_val}, p = {p_val}"))
#   plt
# }
# 
# smooth_table
# ggarrange(plot_partial(1), plot_partial(2))
# # plot_partial(5), plot_partial(6))
# 
# ggsave("results/modelling_out/viral_diversity_GAM_partial_effects.png", dpi = 300, 
#        height = 6, width = 10)
# 
# merged_df %>%
#   ggplot(aes(x = generation_length_d, y = corr_diversity)) +
#   geom_point() +
#   geom_smooth()
# 
# merged_df %>%
#   arrange(desc(adult_mass_g))
# summary(gam(corr_diversity ~ s(adult_mass_g) + s(brain_mass_g) + s(interbirth_interval_d) + s(generation_length_d), 
#             data = merged_df))
# 
# meta %>% 
#   distinct(host_genus, host_class, host_order, host_family) %>%
#   filter(host_genus %in% tolower(merged_df$host_genus)) %>%
#   group_by(host_class) %>%
#   mutate(host_class = ifelse(host_class %in% c("", "lepidosauria"), "reptilia", host_class)) %>%
#   mutate(host_class = capitalize(host_class)) %>%
#   summarise(n = n()) %>%
#   mutate(host_class = capitalize(host_class)) %>%
#   arrange(desc(n))
# 
# meta %>% 
#   distinct(host_genus, host_class, host_order, host_family) %>%
#   filter(host_genus %in% tolower(merged_df$host_genus)) %>%
#   filter(host_class == "")


# # Boostrap cross-validation
# k <- 1000
# set.seed(66)
# cv_morsels <- foreach(i = seq(k),
#                       .combine = "c") %do% {
#                         train_fold <- merged_df %>%
#                           sample_n(nrow(merged_df), replace = T)
#                         test_fold <- merged_df %>%
#                           filter(!(host_genus %in% train_fold$host_genus))
#                         
#                         gam_temp <- gam(as.formula(f_string), 
#                                         data = train_fold)
#                         
#                         y <- test_fold$corr_diversity
#                         y_pred <- predict(gam_temp, test_fold)
#                         abs_error <- abs(y - y_pred) / length(y)
#                         mae <- mean(abs_error, na.rm = T)
#                         mae
#                       }
# 
# median_response <- median(abs((gmm$y)))
# 
# # Plot MAE relative to y_mean
# tibble(x = cv_morsels) %>%
#   ggplot(aes(x = x)) + 
#   geom_histogram() +
#   geom_vline(xintercept = median_response,
#              color = "indianred",
#              lty = "dashed") +
#   geom_text(aes(x = median_response + 0.15,
#                 y = 500,
#                 label = str_glue("Median response value = {round(median_response, 2)}")),
#             color = "indianred") +
#   theme_classic2() +
#   labs(x = "Mean absolute error", y = "No. of bootstraps")
# 
# ggsave("results/modelling_out/viral_diversity/bootstrap_MAE.PCA.png", width = 8, height = 5)
