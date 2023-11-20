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

file_names <- list.files("results/dnds_out/family_gene_alignments/", "*.aln",
                         full.names = T)
file_names <- file_names[!grepl("log", file_names)]

fna <- readDNAStringSet(file_names[2])

# Trim gaps to reference
ref_name <- names(fna)[1]
fna_mat <- as.matrix(fna)
to_keep <- fna_mat[ref_name, ] != "-"

# Remove reference
fna_filt <- fna_mat[2:nrow(fna_mat), to_keep]

# Convert back to DNAss
fna_filt <- DNAStringSet(apply(fna_filt, 1, paste0, collapse = ""))
aln_length <- unique(width(fna_filt))

# window <- 150
# 
# cl <- makeCluster(12)
# registerDoParallel(cl)
# 
# morsels <- foreach(i = seq(1, length(fna_filt) - 1, 2),
#                    .packages = c("Biostrings", "tidyverse", "foreach",
#                                  "MSA2dist")) %dopar% {
#   # print(i)
#   fna_temp <- fna_filt[i:(i + 1)]
#   seq_names <- str_split(names(fna_temp), "\\|", simplify = T)[, 1]
# 
#   crumbs <- foreach(j = seq(1, aln_length - window + 1, window),
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
# # bind_rows(morsels) %>%
# #   ggplot(aes(x = start_pos, y = tip_name, fill = log10(ka))) + 
# #   geom_tile()
# stopCluster(cl)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  distinct(tip_name, is_jump)