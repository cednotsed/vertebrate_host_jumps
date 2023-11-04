rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(Hmisc)
require(doParallel)

pan_dir <- "results/dnds_out/panaroo_out"
fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

runs <- list.files(pan_dir)
runs <- runs[grepl("seqid", runs)]
runs <- runs[grepl("relax", runs)]
runs
list.files(runs)
# cl <- makeCluster(5)
# registerDoParallel(cl)

morsels <- foreach(run = runs,
                   .packages = c("tidyverse", "data.table",
                                 "Biostrings", "foreach")) %do% {
  print(run)
  # run = runs[1]
  run_path <- str_glue("{pan_dir}/{run}")
  cliques <- list.files(run_path)
  
  seqid <- gsub("seqid_", "", paste0(str_split(run, "\\.")[[1]][1:2], 
                                     collapse = "."))
  protid <- gsub("prot_", "", paste0(str_split(run, "\\.")[[1]][3:4], 
                                     collapse = "."))
  seqid
  protid
  
  crumbs <- foreach(clique = cliques) %do% {
    print(clique)
    # clique = cliques[1]
    clique_path <- str_glue("{run_path}/{clique}")
    
    core_path <- str_glue("{clique_path}/core_gene_alignment.aln")
    sum_path <- str_glue("{clique_path}/summary_statistics.txt")
    
    if(file.exists(sum_path)) {
      # Count genes
      sum_df <- fread(sum_path) 
      total_genes <- deframe(sum_df %>% 
        filter(V1 == "Total genes") %>%
        select(V3))
      
      n_core <- deframe(sum_df %>% 
                         filter(V1 %in% c("Core genes", "Soft core genes")) %>%
                         summarise(n_core = sum(V3)))
      
      prop_core <- n_core / total_genes
      prop_core
      
      if (n_core > 0) {
        core_aln <- readDNAStringSet(core_path,
                                     format = "fasta")
        # Parse accession names
        names(core_aln) <- gsub("_R_", "", names(core_aln))
        
        # Get prop gaps
        full_fna <- fna[names(core_aln)]
        core_mat <- as.matrix(core_aln)
        gapped_pos <- apply(core_mat, 2, function(x) {(sum(x == "-") / length(x)) > 0.1})
        n_gapped_pos <- sum(gapped_pos)
        n_ungapped_pos <- sum(!gapped_pos)
        
        prop_ungapped <- n_ungapped_pos / width(core_aln)[1]
        
        # Get prop genome recovered
        prop_recovered <- n_ungapped_pos / median(width(full_fna))
        
      } else {
        prop_recovered <- 0
        prop_ungapped <- 0
      }
      
      return(tibble(cluster = clique,
                    seqid = seqid,
                    prot_id = protid,
                    total_genes = total_genes,
                    n_core = n_core,
                    prop_core_genes = prop_core,
                    prop_recovered = prop_recovered,
                    prop_ungapped = prop_ungapped))
    }
  }
  bind_rows(crumbs)
}

# stopCluster(cl)

bind_rows(morsels) %>%
  group_by(seqid, prot_id) %>%
  summarise(n = n(),
            n_good = sum(prop_ungapped > 0.8),
            median_core = median(prop_core_genes),
            median_ungapped = median(prop_ungapped),
            median_recovered = median(prop_recovered)) %>%
  ggplot(aes(x = seqid, y = n_good, size = median_ungapped,
             color = median_recovered, shape = prot_id)) +
  geom_point(alpha = 0.5) +
  # facet_grid(rows = vars(prot_id)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0.87) +
  labs(x = "Nuc. identity threshold",
       y = "No. good alignments",
       size = "Med. prop. ungapped positions",
       color = "Med. prop. genome recovered") +
  theme_bw()

bind_rows(morsels) %>% 
  # filter(prop_ungapped < 0.5) %>%
  filter(seqid == 0.6, prot_id == 0.7) %>%
  View()

bind_rows(morsels) %>%
  group_by(seqid, prot_id) %>%
  summarise(n = n(),
            n_gappy = sum(prop_ungapped < 0.8),
            median_core = median(prop_core_genes),
            median_ungapped = median(prop_ungapped),
            median_recovered = median(prop_recovered))

score_res <- bind_rows(morsels) %>%
  mutate(aln_score = (prop_ungapped + prop_recovered) / 2)

score_morsels <- foreach(clique_name = unique(score_res$cluster)) %do% {
  score_res %>%
    filter(cluster == clique_name) %>%
    arrange(desc(prop_ungapped)) %>%
    head(1)
}

bind_rows(score_morsels) %>%
  distinct(cluster)
  filter(prop_ungapped >= 0.8) %>%
  distinct(cluster) %>% View()


# 
# score_res <- group_by(cluster) %>%
#   summarise(max_score = max(aln_score)) %>%
#   left_join(bind_rows(morsels) %>% mutate(aln_score = (prop_ungapped + prop_recovered) / 2))
#   
