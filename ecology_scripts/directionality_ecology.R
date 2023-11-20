rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Hmisc)
require(mgcv)
require(gratia)
require(ggpubr)
require(foreach)
require(doParallel)

genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv") %>%
  rename(viral_clique = cluster) %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv") %>%
              distinct(host, host_genus, host_class))

# Sequencing effort
effort_df <- genome_meta %>%
  filter(host_genus != "") %>%
  group_by(host_genus, host_class) %>%
  summarise(n_genomes = n_distinct(accession)) %>%
  arrange(desc(n_genomes)) %>%
  filter(!is.na(host_class))

good_alns <- fread("results/qc_out/good_alignments.csv")
jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  distinct(anc_genus, tip_genus, anc_name, clique_name) %>%
  filter(tip_genus != anc_genus) %>%
  group_by(anc_genus, tip_genus) %>%
  summarise(n_jumps = n()) %>%
  filter(anc_genus != "",
         tip_genus != "")

traits_df <- fread("data/ecology/combine_data_soria_2021/trait_data_imputed.csv",
                   stringsAsFactors = T) %>%
  rename_all(~tolower(gsub(" |-", "_", .x))) %>%
  rename_all(~gsub("-_", "", .x)) %>%
  dplyr::rename(host_genus = genus)

traits_parsed <- traits_df %>%
  group_by(host_genus) %>%
  summarise(across(c(ends_with("_g"), 
                     ends_with("_d"), 
                     ends_with("_mm"),
                     contains("dispersal")), 
                   ~ median(.x, na.rm = T)))

merged_df <- jump_df %>%
  left_join(traits_parsed %>% select(tip_genus = host_genus,
                                     tip_longevity = max_longevity_d,
                                     tip_mass = adult_mass_g,
                                     tip_gen = generation_length_d,
                                     tip_dispersal = dispersal_km)) %>%
  left_join(traits_parsed %>% select(anc_genus = host_genus,
                                     anc_longevity = max_longevity_d,
                                     anc_mass = adult_mass_g,
                                     anc_gen = generation_length_d,
                                     anc_dispersal = dispersal_km)) %>%
  left_join(effort_df %>% select(anc_genus = host_genus, 
                                   anc_genomes = n_genomes)) %>%
  left_join(effort_df %>% select(tip_genus = host_genus, 
                                   tip_genomes = n_genomes)) %>%
  mutate(rel_longevity = anc_longevity / tip_longevity,
         rel_mass = anc_mass / tip_mass,
         rel_gen = anc_gen / tip_gen,
         rel_dispersal = anc_dispersal / tip_dispersal,
         mean_effort = 0.5 * (anc_genomes + tip_genomes)) %>%
  mutate()
merged_df %>% View()


final <- merged_df %>%
  filter(!is.na(rel_longevity),
         !is.na(rel_mass),
         !is.na(rel_gen),
         !is.na(rel_dispersal))

final %>%
  ggplot(aes(x = log(rel_longevity), y = n_jumps)) +
  geom_point()

lr <- glm(n_jumps ~ log(mean_effort) + log(rel_mass) + log(rel_gen) + log(rel_longevity) + log(rel_dispersal),
    data = final)

summary(lr)

sum(is.na(merged_df$rel_mass))
# Correct for clique representation and sequencing effort
# Visualise confounders
diversity_df %>%
  select(-host_genus, -host_class) %>%
  pivot_longer(!n_cliques, names_to = "confounders", values_to = "value") %>%
  mutate(confounders = case_when(confounders == "effort" ~ "log10(sum bases sequenced)")) %>%
  ggplot(aes(x = value, y = log10(n_cliques))) +
  facet_grid(cols = vars(confounders),
             scale = "free") +
  geom_point() +
  geom_smooth(method = "lm")

# ggsave("results/modelling_out/viral_diversity/genome_data/confounders_linear_visualisation.png", width = 6, height = 4)

# Get residuals
diversity_lm <- lm(log10(n_cliques) ~ n_genomes,
                   data = diversity_df)

corr_diversity_df <- diversity_df %>%
  mutate(corr_diversity = diversity_lm$residuals)

parsed_df <- traits_df %>%
  group_by(host_genus) %>%
  summarise(across(c(ends_with("_g"), 
                     ends_with("_d"), 
                     ends_with("_mm")), 
                   ~ median(.x, na.rm = T))) %>%
  inner_join(diversity_df)

pr <- glm(n_cliques ~ log(n_genomes) + adult_mass_g + brain_mass_g + max_longevity_d + maturity_d + age_first_reproduction_d + gestation_length_d + generation_length_d,
          data = parsed_df,
          family = "poisson")

mat <- parsed_df %>%
  column_to_rownames("host_genus") %>%
  select(-host_class, -n_cliques, -n_genomes) %>%
  mutate(across(everything(), log))

missing <- apply(mat, 2, function(x) sum(is.na(x)))
mat <- mat[, names(missing[missing == 0])]
pca <- prcomp(mat, retx = T, center = T, scale. = T)

pca$sdev^2 / sum(pca$sdev ^2)

as.data.frame(pca$rotation) %>%
  rownames_to_column("host_genus") %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = host_genus))
cor.test(parsed_df$adult_mass_g, parsed_df$interbirth_interval_d)

parsed_df %>%
  bind_cols(pca$x) %>%
  ggplot(aes(x = n_cliques, y = PC2)) +
  geom_point()
as.data.frame(pca$x) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point()
parsed_df %>%
  ggplot(aes(x = log10(brain_mass_g))) +
  geom_histogram()
summary(pr)
require(mgcv)

lr <- glm(n_cliques ~ log(n_genomes) + log(adult_mass_g) + log(brain_mass_g),
          data = parsed_df,
          family = poisson())

summary(lr)    
# + max_longevity_d + maturity_d + age_first_reproduction_d + gestation_length_d + generation_length_d,
parsed_df %>%
  filter(is.na(brain_mass_g))
parsed_df$brain_mass_g
diversity_df %>%
  left_join(traits_df)
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
  select(all_of(to_keep)) %>%
  group_by(host_genus) %>%
  summarise_if(is.numeric, median, na.rm = T)

merged_df <- corr_diversity_df %>%
  left_join(traits_filt)

# filter(!(host_genus %in% c("Homo", "Sus", "Equus", "Mus", "Felis", "Canis")))

predictors <- colnames(merged_df)
predictors <- predictors[!(predictors %in% c("host_genus", "host_family", "host_order",
                                             "librarystrategy", "n_cliques", "norm_n_cliques",
                                             "effort", "n_cliques_rep", "corr_diversity",
                                             "host_class"))]

f_string <- "corr_diversity ~"
for (i in seq(length(predictors))) {
  predictor <- predictors[i]
  if (i == 1) {
    f_string <- str_glue("{f_string}s({predictor}, k = 5)") 
  } else {
    f_string <- str_glue("{f_string} + s({predictor}, k = 5)") 
  }
}

# f_string <- "corr_diversity ~ "
# for (i in seq(length(predictors))) {
#   predictor <- predictors[i]
#   if (i == 1) {
#     f_string <- str_glue("{f_string}{predictor}") 
#   } else {
#     f_string <- str_glue("{f_string} + {predictor}") 
#   }
# }

f_string
gmm <- gam(as.formula(f_string), 
           data = merged_df,
           select = T)

res <- summary(gmm)
obs_dev_explained <- res$dev.expl

# draw(gmm, residuals = T)

# Permutation test
perm_iters <- 1000

set.seed(66)
perm_morsels <- foreach(i = seq(perm_iters), 
                        .combine = "c") %do% {
                          temp_dat <- merged_df %>%
                            mutate(corr_diversity = sample(merged_df$corr_diversity, replace = F))
                          
                          gam_temp <- gam(as.formula(f_string), 
                                          data = temp_dat, 
                                          select =T)
                          
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

ggsave("results/modelling_out/viral_diversity/genome_data/permutation_test.png", width = 8, height = 5)

# # Boostrap cross-validation
# k <- 1000
# set.seed(66)
# 
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
# ggsave("results/modelling_out/viral_diversity/genome_data/bootstrap_MAE.png", width = 8, height = 5)



# merged_df %>%
#   pivot_longer(!c(host_genus, n_cliques, log_sum_bases, n_cliques_rep, corr_diversity), 
#                names_to = "variables", values_to = "value") %>%
#   ggplot(aes(x = corr_diversity, y = value)) +
#   facet_grid(rows = vars(variables), scales = "free") +
#   geom_point() +
#   geom_smooth() 
#   ggplot(aes(log10(corr_diversity)))
# Plot all partial plots
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
