rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)

# Metadata
genome_type <- fread("data/metadata/genome_type_metadata.csv")
good_alns <- fread("results/qc_out/good_alignments.csv")
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.new.csv") %>%
  left_join(fread("data/metadata/parsed_host_metadata.csv")) %>%
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

# Jumps
jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  filter(clique_name %in% good_alns$clique_name) %>%
  mutate(is_jump = T)

# non-jumps
# nonjump_df <- fread("results/mutational_load_out/non_jumps.csv") %>%
#   filter(clique_name %in% good_alns$clique_name) %>%
#   mutate(is_jump = F)

# n_distinct(jump_df$clique_name) == n_distinct(nonjump_df$clique_name)

# Get matching non jumps
nonjump_df <- fread("results/mutational_load_out/non_jumps.csv") %>%
  filter(clique_name %in% jump_df$clique_name) %>%
  # Remove tips that are involved in host jumps
  filter(!(tip_name %in% jump_df$tip_name))

cl <- makeCluster(16)
registerDoParallel(cl)

tip_list <- deframe(nonjump_df %>%
                      distinct(tip_name))

set.seed(69)

nonjump_morsels <- foreach(tip = tip_list, .packages = c("tidyverse")) %dopar% {
  nonjump_temp <- nonjump_df %>%
    filter(tip_name == tip) %>%
    sample_n(1, replace = F)
}

stopCluster(cl)

nonjump_filt <- bind_rows(nonjump_morsels) %>%
  mutate(is_jump = F)

# Merge non-jumps and jumps
merged_df <- bind_rows(nonjump_filt, jump_df) %>%
  left_join(host_counts) %>%
  separate(clique_name, c("family"), "_", remove = F) %>%
  left_join(genome_type) %>%
  filter(clique_name %in% nonjump_df$clique_name)

# Check if all cliques have jumps and nonjumps
merged_df %>% 
  group_by(clique_name) %>% 
  summarise(n = n_distinct(is_jump)) %>%
  arrange(n)

fwrite(merged_df, "results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.V2.csv")
# Merge distances from finished run because sampling is non-deterministic
# test <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps_OLD.csv") %>%
#   select(-patristic_dist)

# jump_dist <- bind_rows(fread("results/mutational_load_out/putative_host_jumps.csv"),
#                        fread("results/mutational_load_out/non_jumps.csv")) %>%
#   select(anc_name, tip_name, n_traverses, patristic_dist)
# 
# test %>%
#   left_join(jump_dist)
# test %>%
#   left_join(jump_dist) %>% 
#   fwrite("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.csv")
# test1 <- merged_df %>% 
#   arrange(n_traverses) %>%
#   select(n_traverses)
# 
# test2 <- test %>%
#   arrange(n_traverses) %>%
#   select(n_traverses)
# sum(test1$n_traverses != test2$n_traverses)
# fwrite(merged_df, "results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.csv")

# Visualise host counts
host_counts %>%
  filter(clique_name %in% merged_df$clique_name) %>%
  ggplot(aes(x = n_hosts)) +
  geom_histogram(bins = 50,
                 fill = "steelblue") +
  labs(x = "No. host genera", y = "No. viral cliques") +
  geom_vline(xintercept = 20, 
             lty = "dashed",
             color = "indianred") +
  theme_classic()

# ggsave("results/mutational_load_out/num_cliques_per_host_count.png", width = 8, height = 5)

clique_counts <- merged_df %>%
  ungroup() %>%
  group_by(n_hosts) %>%
  summarise(n_cliques = n_distinct(clique_name)) %>%
  ungroup() %>%
  filter(n_cliques >= 5)

# plot_df <- merged_df %>%
#   group_by(clique_name, is_jump, n_hosts) %>%
#   summarise(min_dist = min(patristic_dist)) %>%
#   left_join(clique_counts) %>%
#   filter(n_cliques >= 5)
# 
# plot_df %>%
#   ggplot(aes(x = factor(n_hosts), y = log10(patristic_dist))) +
#   # geom_point(alpha = 0.5) +
#   geom_boxplot(aes(fill = is_jump)) +
#   theme_classic() +
#   labs(x = "No. host genera", y = "log10(min. dist.)",
#        fill = "Is host jump?") +
#   geom_text(aes(x = n_hosts, y = -3, label = str_glue("n={n_cliques}")),
#             data = clique_counts,
#             color = "grey39",
#             size = 3)
# 
# ggsave("results/mutational_load_out/mutational_threshold_visualisation.png", width = 5, height = 5)
# 
# 
# lr1_jump <- glm(n_hosts ~ min_dist,
#                 data = plot_df %>% filter(is_jump),
#                 family = poisson())
# 
# lr2_jump <- glm(n_hosts ~ log(min_dist),
#                 data = plot_df %>% filter(is_jump),
#                 family = poisson())
# 
# lr3_jump <- glm(n_hosts ~ log10(min_dist),
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
# lr1_non <- glm(n_hosts ~ min_dist,
#                data = plot_df %>% filter(!is_jump),
#                family = poisson())
# 
# lr2_non <- glm(n_hosts ~ log10(min_dist),
#                data = plot_df %>% filter(!is_jump),
#                family = poisson())
# 
# lr3_non <- glm(n_hosts ~ log(min_dist),
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
# 
# # Permutation test
# set.seed(66)
# 
# perm_morsels <- foreach(i = seq(1000), .combine = "c") %do% {
#   perm_temp <- plot_df %>%
#     group_by(clique_name) %>%
#     mutate(is_jump = sample(is_jump, replace = F))
#   
#   lr1_jump <- glm(n_hosts ~ log(min_dist),
#                   data = perm_temp %>% filter(is_jump),
#                   family = poisson())
#   
#   lr1_non <- glm(n_hosts ~ log(min_dist),
#                  data = perm_temp %>% filter(!is_jump),
#                  family = poisson())
#   
#   jump_coef <- coef(lr1_jump)[[2]]
#   nonjump_coef <- coef(lr1_non)[[2]]
#   
#   slope_diff <- nonjump_coef - jump_coef
#   
#   return(slope_diff)
# }
# 
# p_val <- tibble(slope_diff = perm_morsels) %>%
#   summarise(p = sum(slope_diff >= obs_slope_diff) / n())
# 
# tibble(slope_diff = perm_morsels) %>%
#   ggplot(aes(x = slope_diff)) +
#   geom_histogram(fill = "steelblue3") +
#   geom_vline(xintercept = obs_slope_diff,
#              lty = "dashed",
#              color = "indianred") +
#   theme_classic() +
#   labs(x = "Diff. in slope estimate", y = "No. permutations",
#        title = str_glue("Permutation test: p = {p_val}"))
# 
# ggsave("results/mutational_load_out/permutation.diff_hosts.genus_counts.png", width = 8, height = 5)
# # nonjump_filt <- bind_rows(nonjump_morsels) %>%
# #   dplyr::rename(non_anc_name = anc_name, 
# #                 non_tip_name = tip_name,
# #                 non_tip_state = tip_state,
# #                 min_dist_non_jump = patristic_dist) %>%
# #   select(-n_traverses, -total_depth)
# # 
# # diff_host <- bind_rows(host_morsels) %>%
# #   dplyr::rename(min_dist_jump = patristic_dist) %>%
# #   inner_join(nonjump_filt) %>%
# #   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
# #   left_join(host_counts) %>%
# #   left_join(genome_counts)
# # 
# # fwrite(diff_host, "results/mutational_load_out/host_jump_lists/diff_host.csv")
# # 
# # 
# # diff_host
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # diff_host <- jump_df %>%
# #   filter(anc_state != tip_state) %>%
# #   filter(anc_state != "" & tip_state != "") %>%
# #   group_by(clique_name, anc_state, tip_state) %>%
# #   summarise(min_dist_jump = min(patristic_dist)) %>%
# #   inner_join(nonjump_filt) %>%
# #   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
# #   left_join(host_counts) %>%
# #   left_join(genome_counts)
# # 
# # diff_genus <- jump_df %>%
# #   filter(anc_genus != tip_genus) %>%
# #   filter(anc_genus != "" & tip_genus != "") %>%
# #   group_by(clique_name, anc_state, tip_state) %>%
# #   summarise(min_dist_jump = min(patristic_dist)) %>%
# #   inner_join(nonjump_filt) %>%
# #   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
# #   left_join(host_counts) %>%
# #   left_join(genome_counts)
# # 
# # diff_family <- jump_df %>%
# #   filter(anc_family != tip_family) %>%
# #   filter(anc_family != "" & tip_family != "") %>%
# #   group_by(clique_name, anc_state, tip_state) %>%
# #   summarise(min_dist_jump = min(patristic_dist)) %>%
# #   inner_join(nonjump_filt) %>%
# #   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
# #   left_join(host_counts) %>%
# #   left_join(genome_counts)
# # 
# # diff_order <- jump_df %>%
# #   filter(anc_order != tip_order) %>%
# #   filter(anc_order != "" & tip_order != "") %>%
# #   group_by(clique_name, anc_state, tip_state) %>%
# #   summarise(min_dist_jump = min(patristic_dist)) %>%
# #   inner_join(nonjump_filt) %>%
# #   mutate(dist_ratio = min_dist_jump / min_dist_non_jump) %>%
# #   left_join(host_counts) %>%
# #   left_join(genome_counts)
# # 
# # # Visualise differences
# # diff_host %>%
# #   select(-dist_ratio, -n_genomes) %>%
# #   pivot_longer(!c(clique_name, anc_state, 
# #                   tip_state, n_hosts,
# #                   non_anc_name, non_tip_name),
# #                names_to = "is_jump", values_to = "pat_dist") %>%
# #   mutate(is_jump = ifelse(grepl("non", is_jump), "Non-jump", "Host jump")) %>%
# #   ggplot(aes(x = n_hosts, y = log10(pat_dist), color = is_jump)) +
# #   geom_point(alpha = 0.5) +
# #   geom_smooth(method = "loess") +
# #   theme_classic() +
# #   labs(x = "No. of host orders", y = "Log10(min. patristic dist.)",
# #        color = "")
# # 
# # ggsave("results/mutational_load_out/jump_versus_no_jump_scatter.png", width = 8, height = 5)
# # dat_list <- list(diff_order, diff_family, diff_genus, diff_host)
# # 
# # for (i in seq(4)) {
# #   dat <- dat_list[[i]]
# # 
# #   lr1 <- glm(n_hosts ~ n_genomes + dist_ratio,
# #              data = dat,
# #              family = poisson())
# #   
# #   lr2 <- glm(n_hosts ~ n_genomes + log(dist_ratio),
# #              data = dat,
# #              family = poisson())
# # 
# #   lr3 <- glm(n_hosts ~ log(n_genomes) + log(dist_ratio),
# #              data = dat,
# #              family = poisson())
# # 
# #   lr4 <- glm(n_hosts ~ log(n_genomes) + dist_ratio,
# #              data = dat,
# #              family = poisson())
# # 
# #   comp <- compareGLM(lr1, lr2, lr3, lr4)
# # 
# #   best_model_idx <- deframe(comp$Fit.criteria %>%
# #                               rownames_to_column("model") %>%
# #                               arrange(AIC) %>%
# #                               select(model) %>%
# #                               head(1))
# # 
# #   best_model_formula <- deframe(comp$Models %>%
# #                                   as_tibble() %>%
# #                                   rownames_to_column("model") %>%
# #                                   filter(model == best_model_idx) %>%
# #                                   select(Formula))
# # 
# #   best_model <- get(paste0("lr", best_model_idx))
# #   obs_effect <- coef(best_model)[[2]]
# # 
# #   print(summary(best_model))
# #   mod_list[[i]] <- best_model
# # }
# # 
# # 
# # permute_coeff(dat_list[[1]], mod_list[[1]], 3)
# # permute_coeff(dat_list[[2]], mod_list[[2]], 3)
# # permute_coeff(dat_list[[3]], mod_list[[3]], 3)
# # permute_coeff(dat_list[[4]], mod_list[[4]], 3)
# # 
# # signif(coef(mod_list[[1]])[[3]], 3)
# # signif(coef(mod_list[[2]])[[3]], 3)
# # signif(coef(mod_list[[3]])[[3]], 3)
# # signif(coef(mod_list[[4]])[[3]], 3)
# # 
# # # Without correction
# # lr1 <- glm(n_hosts ~ dist_ratio,
# #            data = diff_host,
# #            family = poisson())
# # 
# # lr2 <- glm(n_hosts ~ log(dist_ratio),
# #            data = diff_host,
# #            family = poisson())
# # 
# # 
# # comp <- compareGLM(lr1, lr2)
# # 
# # best_model_idx <- deframe(comp$Fit.criteria %>%
# #                             rownames_to_column("model") %>%
# #                             arrange(AIC) %>%
# #                             select(model) %>%
# #                             head(1))
# # 
# # best_model_formula <- deframe(comp$Models %>%
# #                                 as_tibble() %>%
# #                                 rownames_to_column("model") %>%
# #                                 filter(model == best_model_idx) %>%
# #                                 select(Formula))
# # 
# # best_model <- get(paste0("lr", best_model_idx))
# # obs_effect <- coef(best_model)[[2]]
# # 
# # print(summary(best_model))
# # 
# # 
# # # fwrite(diff_order, "results/mutational_load_out/host_jump_lists/diff_order.csv")
# # # fwrite(diff_family, "results/mutational_load_out/host_jump_lists/diff_family.csv")
# # # fwrite(diff_genus, "results/mutational_load_out/host_jump_lists/diff_genus.csv")
# # # fwrite(diff_host, "results/mutational_load_out/host_jump_lists/diff_host.csv")
