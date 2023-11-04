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

cl <- makeCluster(16)
registerDoParallel(cl)

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
  filter(!(ka == 0 & ks == 0)) %>%
  filter(ka == 0 | ks == 0) %>%
  left_join(jump_df) %>%
  mutate(type = ifelse(ka == 0, "Only syn", "Only non-syn")) %>%
  mutate(domain = factor(domain, c("NTD", "RBD", "FP", "HR1-HR2", "TM-CT"))) %>%
  mutate(is_rbd = ifelse(domain == "RBD", T, F)) %>%
  group_by(is_rbd, is_jump, type) %>%
  summarise(n = n())
  
# plot_df %>%
#   group_by(is_rbd, is_jump) %>%
#   summarise(ratio = sum(type == "Only non-syn") / n()) %>%
#   ggplot(aes(x = is_rbd, y = ratio, fill = is_rbd)) +
#   facet_grid(.~ is_jump) +
#   geom_bar(stat = "identity") 

plot_df %>%
  ggplot(aes(x = is_rbd, y = n, fill = type)) +
  facet_grid(.~ is_jump) +
  geom_bar(stat = "identity", position = "dodge",
           color = "black") +
  theme_classic() +
  labs(x = "Is RBD?", y = "Frequency", fill = "Mutations")

chi_df <- plot_df %>%
  ungroup() %>%
  filter(!is_jump) %>%
  select(-is_jump) %>%
  pivot_wider(id_cols = is_rbd, names_from = type, values_from = n) %>%
  column_to_rownames("is_rbd")

fisher.test(chi_df)

ggsave("results/dnds_out/family_plots/CoV2_dnds.boxplot.zeroes.png", dpi = 600, width = 6, height = 4)
