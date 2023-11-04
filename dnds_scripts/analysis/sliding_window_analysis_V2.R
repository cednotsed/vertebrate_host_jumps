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

fna <- readDNAStringSet("results/dnds_out/family_gene_alignments/CoV2_spike.aln")

ref <- "Wuhan-Hu-1 spike"
which(names(fna) == ref)

# Trim gaps to wuhan-hu-1
fna_mat <- as.matrix(fna)

to_keep <- fna_mat[ref, ] != "-"

fna_filt <- fna_mat[2:nrow(fna_mat), to_keep]
fna_filt <- DNAStringSet(apply(fna_filt, 1, paste0, collapse = ""))

window <- 150

# cl <- makeCluster(16)
# registerDoParallel(cl)

# morsels <- foreach(i = seq(1, length(fna_filt) - 1, 2), 
#                    .packages = c("Biostrings", "tidyverse", "foreach",
#                                  "MSA2dist")) %dopar% {
#   # print(i)
#   fna_temp <- fna_filt[i:(i + 1)]
#   seq_names <- str_split(names(fna_temp), "\\|", simplify = T)[, 1]
#   
#   crumbs <- foreach(j = seq(1, 3822 - window + 1, window), 
#                     .combine = "bind_rows", 
#                     .packages = c("Biostrings", "tidyverse")) %do% {
#     # j = 3523
#     window_aln <- subseq(fna_temp, j, j + window - 1)
#     dnds <- dnastring2kaks(window_aln,
#                            model = "Li")
#     ka <- as.numeric(dnds[, "ka"])
#     ks <- as.numeric(dnds[, "ks"])
#     
#     return(tibble(anc_name = seq_names[2],
#                   tip_name = seq_names[1],
#                   ka = ka,
#                   ks = ks,
#                   start_pos = j))
#   }
#   return(crumbs)
# }
# 
# stopCluster(cl)

jump_df <- fread("results/dnds_out/all_jumps.per_gene_dnds.diff_hosts.genus_counts.csv") %>%
  distinct(tip_name, is_jump)

annot_df <- tibble(domain = c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"),
                   start_pos = c(1, 319, 542, 912, 1214),
                   end_pos = c(318, 541, 911, 1212, 1274)) %>%
  mutate(start_pos = start_pos * 3 - 2,
         end_pos = end_pos * 3)

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
  mutate(log_kaks = ifelse(ka == 1e-6 & ks == 1e-6, NA, log10(ka / ks))) %>%
  arrange(ka)

plot_df %>%
  left_join(jump_df) %>%
  mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT")),
         tip_name = factor(tip_name, unique(plot_df$tip_name))) %>%
  ggplot(aes(x = domain, y = tip_name, fill = log_kaks)) +
  geom_tile() +
  facet_grid(is_jump~.,
             scale = "free") +
  scale_fill_gradient2(low = "blue",
                       mid = "#FFF9c4",
                       high = "red",
                       na.value = NA,
                       midpoint = 0) +
  theme_classic() +
  labs(x = "Spike domain", 
       y = "log10(dn/ds)") +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

ggsave("results/dnds_out/family_plots/CoV2_dnds.heatmap.pdf", dpi = 600, width = 5, height = 7)

plot_df %>%
  left_join(jump_df) %>%
  mutate(is_rbd = ifelse(domain == "RBD", T, F)) %>%
  ggplot(aes(x = is_rbd, y = log_kaks, fill = is_rbd)) +
  facet_grid(is_jump ~ .) +
  geom_boxplot() + 
  theme_classic() +
  geom_pwc(na.rm = T) +
  theme(legend.position = "none") +
  labs(x = "Is RBD?", 
       y = "log10(dn/ds)") +
  coord_flip()

# Get statistics
test <- plot_df %>%
  left_join(jump_df) %>%
  mutate(is_rbd = ifelse(domain == "RBD", T, F)) %>%
  filter(!is_jump) %>%
  filter(!is.na(log_kaks))

wilcox.test(test$log_kaks ~ test$is_rbd)

ggsave("results/dnds_out/family_plots/CoV2_dnds.boxplot.pdf", dpi = 600, width = 5, height = 7)
  
plot_parsed <- bind_rows(morsels) %>%
  mutate(ka = ifelse(ka <= 0, NA, ka),
         ks = ifelse(ks <= 0, NA, ks)) %>%
  mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
  left_join(jump_df) %>%
  mutate(log_kaks = log10(ka / ks))

plot_parsed %>%
  mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
  left_join(jump_df) %>%
  ggplot(aes(x = domain, y = tip_name, fill = log_kaks)) +
  geom_tile() +
  facet_grid(is_jump~.,
             scales = "free") +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       na.value = NA,
                       midpoint = 0) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())

ggsave("results/dnds_out/family_plots/CoV2_dnds.pdf", dpi = 600, width = 8, height = 4)
 
plot_df %>%
  mutate(log_kaks = ifelse(ka == 1e-6 & ks == 1e-6, NA, log10(ka / ks))) %>%
  ggplot(aes(x = domain, y = log_kaks)) +
  geom_violinhalf()
  # mutate(ka = ifelse(ka <= 0, NA, ka),
         # ks = ifelse(ks <= 0, NA, ks)) %>%
  mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
  left_join(jump_df) %>%
  ggplot(aes(x = domain, y = tip_name, fill = log_kaks)) +
  geom_tile() +
  facet_grid(is_jump~.,
             scales = "free") +
  scale_fill_gradient2(low = "blue", 
                       mid = "#FFF9c4",
                       high = "red", 
                       na.value = NA,
                       midpoint = 0.02) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())

ggsave("results/dnds_out/family_plots/CoV2_dn.pdf", dpi = 600, width = 4, height = 4)

bind_rows(morsels) %>%
  mutate(ka = ifelse(ka <= 0, NA, ka),
         ks = ifelse(ks <= 0, NA, ks)) %>%
  mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
  left_join(jump_df) %>%
  ggplot(aes(x = domain, y = tip_name, fill = ks)) +
  geom_tile() +
  facet_grid(is_jump~.,
             scales = "free") +
  scale_fill_gradient2(low = "blue", 
                       mid = "#FFF9c4",
                       high = "red", 
                       na.value = NA,
                       midpoint = 0.02) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())

ggsave("results/dnds_out/family_plots/CoV2_ds.pdf", dpi = 600, width = 4, height = 4)
