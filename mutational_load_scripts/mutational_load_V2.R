rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(see)

genome_type <- fread("data/metadata/genome_type_metadata.csv")

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv")) %>%
  left_join(genome_type)

jump_df <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.V2.csv")

# Counts
jump_df %>%
  filter(!is_jump) %>%
  distinct(clique_name) %>%
  nrow()

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

clique_counts <- jump_df %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup()

merged_df <- jump_df %>%
  group_by(clique_name, is_jump) %>%
  summarise(min_dist = min(patristic_dist)) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  left_join(host_counts) %>%
  left_join(genome_counts) %>%
  left_join(clique_counts)

# Visualise host counts
# host_counts %>%
#   filter(clique_name %in% merged_df$clique_name) %>%
#   ggplot(aes(x = n_hosts)) +
#   geom_histogram(bins = 50,
#                  fill = "steelblue") +
#   labs(x = "No. host genera", y = "No. viral cliques") +
#   geom_vline(xintercept = 20, 
#              lty = "dashed",
#              color = "indianred") +
#   theme_classic()
# 
# ggsave("results/mutational_load_out/num_cliques_per_host_count.png", width = 8, height = 5)

# Visualise overall mutational distance
wilcox <- wilcox.test(log10(merged_df$min_dist) ~ merged_df$is_jump,
            alternative = "less")
table(merged_df$is_jump)
merged_df %>%
  ggplot(aes(x = is_jump, y = log10(min_dist), fill = is_jump)) +
  geom_boxplot(position = position_nudge(x = 0.05, y = 0),
               width = 0.1,
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.1, y = 0), alpha = 1) +
  theme_classic() +
  coord_flip() +
  geom_pwc() +
  labs(x = "Is host jump?", y = "log10(min. distance)") +
  theme(legend.position = "none",
        text = element_text(color = "black")) +
  scale_fill_manual(values = c("steelblue3", "indianred3"))

ggsave("results/mutational_load_out/overall_jump_v_non_jump_dist.pdf", width = 5, height = 3.5)

# Correcting for family and n_genomes
logreg <- glm(is_jump ~ log(n_genomes) + family + log(min_dist),
              data = merged_df,
              family = "binomial")

nrow(merged_df) - nrow(coef(summary(logreg)))
coef(summary(logreg))

exp(coef(logreg)[["log(min_dist)"]])

merged_df %>%
  ggplot(aes(x = factor(n_hosts), y = log10(min_dist), fill = is_jump)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    labs(x = "No. host genera", y = "log10(min. dist.)",
         fill = "Is host jump?") +
  scale_fill_manual(values = c("steelblue3", "indianred3"))

ggsave("results/mutational_load_out/mutational_threshold_visualisation.pdf", width = 6, height = 3.5)

# Test for relationship between sampling effort and host range
cor.test(merged_df$n_genomes, merged_df$n_hosts)

# Test for difference in slopes
formula_string <- "n_hosts ~ log(n_genomes) + family + log(min_dist)"
jump_mod <- glm(as.formula(formula_string),
                data = merged_df %>% filter(is_jump),
                family = poisson())
nrow(merged_df %>% filter(is_jump)) - nrow(coef(summary(jump_mod)))
nonjump_mod <- glm(as.formula(formula_string),
                data = merged_df %>% filter(!is_jump),
                family = poisson())
nrow(merged_df %>% filter(!is_jump)) - nrow(coef(summary(nonjump_mod)))
obs_slope_diff <- coef(jump_mod)[["log(min_dist)"]] - coef(nonjump_mod)[["log(min_dist)"]]
length(coef(jump_mod))
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

  slope_diff <- coef(lr1_jump)[["log(min_dist)"]] - coef(lr1_non)[["log(min_dist)"]]

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

ggsave("results/mutational_load_out/permutation.diff_hosts.genus_counts.png", width = 8, height = 5)


# lr1_jump <- lm(n_hosts ~ log10(n_genomes) + family + log(min_dist),
#                 data = plot_df %>% filter(is_jump))
# 
# summary(lr1_jump)
# lr1_jump <- glm(n_hosts ~ family + min_dist,
#                 data = plot_df %>% filter(is_jump),
#                 family = poisson())
# 
# lr2_jump <- glm(n_hosts ~ family + log(patristic_dist),
#            data = merged_df %>% filter(!is_jump),
#            family = poisson())
# 
# lr3_jump <- glm(n_hosts ~ family + log10(min_dist),
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
# obs_effect_jump <- coef(best_model)[["log(min_dist)"]]
# 
# print(summary(best_model))
# 
# # Test non jump
# # lr1_non <- glm(n_hosts ~ family + min_dist,
# #                 data = plot_df %>% filter(!is_jump),
# #                 family = poisson())
# # 
# # lr2_non <- glm(n_hosts ~ family + log10(min_dist),
# #                 data = plot_df %>% filter(!is_jump),
# #                 family = poisson())
# 
# lr3_non <- glm(n_hosts ~ family + log(min_dist),
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
# best_model_idx <- 3
# best_model <- get(str_glue("lr{best_model_idx}_non"))
# obs_effect_nonjump <- coef(best_model)[["log(min_dist)"]]
# 
# print(summary(best_model))
# 
# obs_slope_diff <- obs_effect_nonjump - obs_effect_jump

# nonjump_filt <- bind_rows(nonjump_morsels) %>%
#   dplyr::rename(non_anc_name = anc_name, 
#                 non_tip_name = tip_name,
#                 non_tip_state = tip_state,
#                 min_dist_non_jump = patristic_dist) %>%
#   select(-n_traverses, -total_depth)
# 
# diff_host <- bind_rows(host_morsels) %>%
#   dplyr::rename(min_dist_jump = patristic_dist) %>%
#   inner_join(nonjump_filt) %>%
#   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
#   left_join(host_counts) %>%
#   left_join(genome_counts)
# 
# fwrite(diff_host, "results/mutational_load_out/host_jump_lists/diff_host.csv")
# 
# 
# diff_host









# 
# diff_host <- jump_df %>%
#   filter(anc_state != tip_state) %>%
#   filter(anc_state != "" & tip_state != "") %>%
#   group_by(clique_name, anc_state, tip_state) %>%
#   summarise(min_dist_jump = min(patristic_dist)) %>%
#   inner_join(nonjump_filt) %>%
#   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
#   left_join(host_counts) %>%
#   left_join(genome_counts)
# 
# diff_genus <- jump_df %>%
#   filter(anc_genus != tip_genus) %>%
#   filter(anc_genus != "" & tip_genus != "") %>%
#   group_by(clique_name, anc_state, tip_state) %>%
#   summarise(min_dist_jump = min(patristic_dist)) %>%
#   inner_join(nonjump_filt) %>%
#   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
#   left_join(host_counts) %>%
#   left_join(genome_counts)
# 
# diff_family <- jump_df %>%
#   filter(anc_family != tip_family) %>%
#   filter(anc_family != "" & tip_family != "") %>%
#   group_by(clique_name, anc_state, tip_state) %>%
#   summarise(min_dist_jump = min(patristic_dist)) %>%
#   inner_join(nonjump_filt) %>%
#   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
#   left_join(host_counts) %>%
#   left_join(genome_counts)
# 
# diff_order <- jump_df %>%
#   filter(anc_order != tip_order) %>%
#   filter(anc_order != "" & tip_order != "") %>%
#   group_by(clique_name, anc_state, tip_state) %>%
#   summarise(min_dist_jump = min(patristic_dist)) %>%
#   inner_join(nonjump_filt) %>%
#   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
#   left_join(host_counts) %>%
#   left_join(genome_counts)
# 
# # Visualise differences
# diff_host %>%
#   select(-dist_ratio, -n_genomes) %>%
#   pivot_longer(!c(clique_name, anc_state, 
#                   tip_state, n_hosts,
#                   non_anc_name, non_tip_name),
#                names_to = "is_jump", values_to = "pat_dist") %>%
#   mutate(is_jump = ifelse(grepl("non", is_jump), "Non-jump", "Host jump")) %>%
#   ggplot(aes(x = n_hosts, y = log10(pat_dist), color = is_jump)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "loess") +
#   theme_classic() +
#   labs(x = "No. of host orders", y = "Log10(min. patristic dist.)",
#        color = "")
# 
# ggsave("results/mutational_load_out/jump_versus_no_jump_scatter.png", width = 8, height = 5)
# dat_list <- list(diff_order, diff_family, diff_genus, diff_host)
# 
# for (i in seq(4)) {
#   dat <- dat_list[[i]]
# 
#   lr1 <- glm(n_hosts ~ n_genomes + dist_ratio,
#              data = dat,
#              family = poisson())
#   
#   lr2 <- glm(n_hosts ~ n_genomes + log(dist_ratio),
#              data = dat,
#              family = poisson())
# 
#   lr3 <- glm(n_hosts ~ log(n_genomes) + log(dist_ratio),
#              data = dat,
#              family = poisson())
# 
#   lr4 <- glm(n_hosts ~ log(n_genomes) + dist_ratio,
#              data = dat,
#              family = poisson())
# 
#   comp <- compareGLM(lr1, lr2, lr3, lr4)
# 
#   best_model_idx <- deframe(comp$Fit.criteria %>%
#                               rownames_to_column("model") %>%
#                               arrange(AIC) %>%
#                               select(model) %>%
#                               head(1))
# 
#   best_model_formula <- deframe(comp$Models %>%
#                                   as_tibble() %>%
#                                   rownames_to_column("model") %>%
#                                   filter(model == best_model_idx) %>%
#                                   select(Formula))
# 
#   best_model <- get(paste0("lr", best_model_idx))
#   obs_effect <- coef(best_model)[[2]]
# 
#   print(summary(best_model))
#   mod_list[[i]] <- best_model
# }
# 
# 
# permute_coeff(dat_list[[1]], mod_list[[1]], 3)
# permute_coeff(dat_list[[2]], mod_list[[2]], 3)
# permute_coeff(dat_list[[3]], mod_list[[3]], 3)
# permute_coeff(dat_list[[4]], mod_list[[4]], 3)
# 
# signif(coef(mod_list[[1]])[[3]], 3)
# signif(coef(mod_list[[2]])[[3]], 3)
# signif(coef(mod_list[[3]])[[3]], 3)
# signif(coef(mod_list[[4]])[[3]], 3)
# 
# # Without correction
# lr1 <- glm(n_hosts ~ dist_ratio,
#            data = diff_host,
#            family = poisson())
# 
# lr2 <- glm(n_hosts ~ log(dist_ratio),
#            data = diff_host,
#            family = poisson())
# 
# 
# comp <- compareGLM(lr1, lr2)
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
# best_model <- get(paste0("lr", best_model_idx))
# obs_effect <- coef(best_model)[[2]]
# 
# print(summary(best_model))
# 
# 
# # fwrite(diff_order, "results/mutational_load_out/host_jump_lists/diff_order.csv")
# # fwrite(diff_family, "results/mutational_load_out/host_jump_lists/diff_family.csv")
# # fwrite(diff_genus, "results/mutational_load_out/host_jump_lists/diff_genus.csv")
# # fwrite(diff_host, "results/mutational_load_out/host_jump_lists/diff_host.csv")


# # Get minimum distance
# host_list <- jump_df %>%
#   distinct(clique_name, anc_state, tip_state)
# 
# host_morsels <- foreach(i = seq(nrow(host_list))) %do% {
#   jump_df %>%
#     filter(clique_name == host_list[i, ]$clique_name,
#            anc_state == host_list[i, ]$anc_state,
#            tip_state == host_list[i, ]$tip_state) %>%
#     arrange(patristic_dist) %>%
#     head(1)
# }
# 
# jump_filt <- bind_rows(host_morsels) %>%
#   mutate(is_jump = T)
# 
# # Get matching non jumps
# nonjump_df <- fread("results/mutational_load_out/non_jumps.csv")
# nonjump_subset <- nonjump_df
# 
# nonjump_morsels <- foreach(i = seq(nrow(jump_filt))) %do% {
#   jump_temp <- jump_filt %>%
#     filter(clique_name == jump_filt[i, ]$clique_name,
#            anc_state == jump_filt[i, ]$anc_state,
#            tip_state == jump_filt[i, ]$tip_state)
#   
#   # Match hosts if possible (CORRECTION NO MATCHING OF HOSTS)
#   nonjump_temp <- nonjump_subset %>%
#     filter(clique_name == jump_filt[i, ]$clique_name
#            # anc_state == jump_filt[i, ]$anc_state
#     ) %>%
#     arrange(patristic_dist) %>%
#     head(1)
#   
#   if (nrow(nonjump_temp) == 0) {
#     return(NULL)
#     # # Get other hosts if not possible
#     # nonjump_temp <- nonjump_subset %>%
#     #   filter(clique_name == jump_filt[i, ]$clique_name) %>%
#     #   arrange(patristic_dist) %>%
#     #   head(1)
#     if(nrow(nonjump_temp) > 0) {
#       # Remove tips from pool
#       nonjump_subset <- nonjump_subset %>%
#         filter(!(tip_name %in% nonjump_temp$tip_name))
#       return(nonjump_temp)
#     } else {
#       return(NULL)      
#     }
#   } else {
#     # Remove tips from pool
#     nonjump_subset <- nonjump_subset %>%
#       filter(!(tip_name %in% nonjump_temp$tip_name))
#     return(nonjump_temp)
#   }
# }
# 
# nonjump_filt <- bind_rows(nonjump_morsels) %>%
#   mutate(is_jump = F)
# 
# (nonjump_filt %>% nrow()) == (nonjump_filt %>% distinct(tip_name) %>% nrow())
