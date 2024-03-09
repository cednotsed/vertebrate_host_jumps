rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(rcompanion)
require(ggpubr)
require(see)
require(Biostrings)
require(MSA2dist)

for (virus in c("IBV", "MERS", "SC2")) {
  # virus <- "MERS"
  print(virus)
  file_names <- list.files("results/dnds_out/family_gene_alignments_V2/", "*.aln",
                           full.names = T)
  file_names <- file_names[!grepl("log", file_names)]
  file_path <- file_names[grepl(virus, file_names)]
  
  fna_filt <- readDNAStringSet(file_path)
  
  jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
    distinct(tip_name, is_jump, patristic_dist)
  
  annot_df <- fread("results/dnds_out/family_gene_alignments/master_annotation_list.csv") %>%
    filter(virus_name == virus)
  
  cl <- makeCluster(12)
  registerDoParallel(cl)
  
  morsels <- foreach(i = seq(1, length(fna_filt) - 1, 2), 
                     .packages = c("Biostrings", "tidyverse", "foreach",
                                   "MSA2dist")) %dopar% {
                                     # print(i)
   fna_temp <- fna_filt[i:(i + 1)]
   seq_names <- str_split(names(fna_temp), "\\|", simplify = T)[, 1]
   
   crumbs <- foreach(j = seq(nrow(annot_df)), 
                     .combine = "bind_rows", 
                     .packages = c("Biostrings", "tidyverse")) %do% {
     # j = 3523
     start_pos <- annot_df[j, ]$start_pos
     end_pos <- annot_df[j, ]$end_pos
     domain <- annot_df[j, ]$domain
     
     window_aln <- subseq(fna_temp, start_pos, end_pos)
     dnds <- dnastring2kaks(window_aln,
                            model = "Li")
     ka <- as.numeric(dnds[, "ka"])
     ks <- as.numeric(dnds[, "ks"])
     
     return(tibble(anc_name = seq_names[2],
                   tip_name = seq_names[1],
                   ka = ka,
                   ks = ks,
                   domain = domain))
     }
     return(crumbs)
  }
  
  stopCluster(cl)
  
  plot_df <- bind_rows(morsels) %>%
    mutate(ka = ifelse(ka <= 0, 1e-6, ka),
           ks = ifelse(ks <= 0, 1e-6, ks)) %>%
    filter(ka != 1e-6 & ks != 1e-6) %>%
    # mutate(log_kaks = ifelse(ka == 0 | ks == 0, NA, log10(ka / ks))) %>%
    mutate(log_kaks = log10((ka) / (ks))) %>%
    # mutate(log_kaks = log10(ks)) %>%
    mutate(virus_name = virus) %>%
    arrange(ka) %>%
    left_join(jump_df)
  # bind_rows(morsels) %>%
  #   filter(domain == "CT") %>%
  #   filter(ka != 0 & ks!= 0)
  # plot_df %>%
  #   ggplot(aes(log(ka + 1), log(ks + 1))) +
  #   geom_point()
  # count_df <- plot_df %>%
  #   filter(!is.na(log_kaks)) %>%
  #   group_by(domain, is_jump) %>%
  #   summarise(n = n_distinct(tip_name))
  
  fwrite(plot_df, str_glue("results/dnds_out/by_domain_V2/{virus}.csv"))
  
  count_df <- plot_df %>%
    mutate(domain = factor(domain, unique(annot_df$domain)),
           is_jump = factor(is_jump)) %>%
    filter(!is.na(log_kaks)) %>%
    group_by(domain, is_jump) %>%
    summarise(n = n_distinct(tip_name)) %>%
    ungroup() %>%
    complete(domain, is_jump) %>%
    mutate(n = replace_na(n, 0))
  
  plot_df %>%
    mutate(domain = factor(domain, unique(annot_df$domain)),
           is_jump = factor(is_jump),
           tip_name = factor(tip_name, unique(plot_df$tip_name))) %>%
    ggplot(aes(x = is_jump, y = log_kaks, fill = is_jump)) +
    geom_boxplot(position = position_dodge(preserve = "single")) +
    geom_pwc(y.position = 1.5) +
    facet_wrap(. ~ domain, 
               nrow = 1,
               strip.position = "bottom") + 
    scale_fill_manual(values = c("steelblue3", "indianred3")) +
    labs(title = virus, 
         x = "Domain",
         y = "Log10(dN/dS)") +
    geom_text(data = count_df,
              aes(x = is_jump, 
                  y = 1.5,
                  label = str_glue("n={n}"), 
                  color = is_jump),
              position = position_dodge(width = 1)) +
    scale_x_discrete(drop = F) + 
    scale_color_manual(values = c("steelblue3", "indianred3")) +
    theme_classic() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(-2, 2)
  
  ggsave(str_glue("results/dnds_out/by_domain_V2/{virus}_spike.pdf"),
         dpi = 600,
         height = 4, width = 8)
}

# 
# plot_df %>%
#   left_join(jump_df) %>%
#   mutate(domain = factor(domain, unique(annot_df$domain)),
#          tip_name = factor(tip_name, unique(plot_df$tip_name))) %>%
#   ggplot(aes(x = domain, y = tip_name, fill = log_kaks)) +
#   geom_tile() +
#   facet_grid(is_jump~.,
#              scale = "free") +
#   scale_fill_gradient2(low = "blue",
#                        mid = "#FFF9c4",
#                        high = "red",
#                        na.value = NA,
#                        midpoint = -0.3) +
#   theme_classic() +
#   labs(x = "Spike domain", 
#        y = "log10(dn/ds)") +
#   theme(legend.position = "top",
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank())

# jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv")
# jump_df %>%
#   filter(clique_name == "Coronaviridae_18") %>%
#   filter(!is_jump)
# ggsave("results/dnds_out/family_plots/CoV2_dnds.heatmap.pdf", dpi = 600, width = 5, height = 7)

# plot_df %>%
#   left_join(jump_df) %>%
#   mutate(is_rbd = ifelse(domain == "NTD", T, F)) %>%
#   ggplot(aes(x = is_rbd, y = log_kaks, fill = is_rbd)) +
#   facet_grid(is_jump ~ .) +
#   geom_boxplot() + 
#   theme_classic() +
#   geom_pwc(na.rm = T) +
#   theme(legend.position = "none") +
#   labs(x = "Is RBD?", 
#        y = "log10(dn/ds)") +
#   coord_flip()
# 
# # Get statistics
# test <- plot_df %>%
#   left_join(jump_df) %>%
#   mutate(is_rbd = ifelse(domain == "RBD", T, F)) %>%
#   filter(!is_jump) %>%
#   filter(!is.na(log_kaks))
# 
# wilcox.test(test$log_kaks ~ test$is_rbd)
# 
# ggsave("results/dnds_out/family_plots/CoV2_dnds.boxplot.pdf", dpi = 600, width = 5, height = 7)
#   
# plot_parsed <- bind_rows(morsels) %>%
#   mutate(ka = ifelse(ka <= 0, NA, ka),
#          ks = ifelse(ks <= 0, NA, ks)) %>%
#   mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
#   left_join(jump_df) %>%
#   mutate(log_kaks = log10(ka / ks))
# 
# plot_parsed %>%
#   mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
#   left_join(jump_df) %>%
#   ggplot(aes(x = domain, y = tip_name, fill = log_kaks)) +
#   geom_tile() +
#   facet_grid(is_jump~.,
#              scales = "free") +
#   scale_fill_gradient2(low = "blue",
#                        mid = "white",
#                        high = "red",
#                        na.value = NA,
#                        midpoint = 0) +
#   theme_classic() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks = element_blank())
# 
# ggsave("results/dnds_out/family_plots/CoV2_dnds.pdf", dpi = 600, width = 8, height = 4)
#  
# plot_df %>%
#   mutate(log_kaks = ifelse(ka == 1e-6 & ks == 1e-6, NA, log10(ka / ks))) %>%
#   ggplot(aes(x = domain, y = log_kaks)) +
#   geom_violinhalf()
#   # mutate(ka = ifelse(ka <= 0, NA, ka),
#          # ks = ifelse(ks <= 0, NA, ks)) %>%
#   mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
#   left_join(jump_df) %>%
#   ggplot(aes(x = domain, y = tip_name, fill = log_kaks)) +
#   geom_tile() +
#   facet_grid(is_jump~.,
#              scales = "free") +
#   scale_fill_gradient2(low = "blue", 
#                        mid = "#FFF9c4",
#                        high = "red", 
#                        na.value = NA,
#                        midpoint = 0.02) +
#   theme_classic() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks = element_blank())
# 
# ggsave("results/dnds_out/family_plots/CoV2_dn.pdf", dpi = 600, width = 4, height = 4)
# 
# bind_rows(morsels) %>%
#   mutate(ka = ifelse(ka <= 0, NA, ka),
#          ks = ifelse(ks <= 0, NA, ks)) %>%
#   mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
#   left_join(jump_df) %>%
#   ggplot(aes(x = domain, y = tip_name, fill = ks)) +
#   geom_tile() +
#   facet_grid(is_jump~.,
#              scales = "free") +
#   scale_fill_gradient2(low = "blue", 
#                        mid = "#FFF9c4",
#                        high = "red", 
#                        na.value = NA,
#                        midpoint = 0.02) +
#   theme_classic() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks = element_blank())
# 
# ggsave("results/dnds_out/family_plots/CoV2_ds.pdf", dpi = 600, width = 4, height = 4)
