rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(treeio)
require(phytools)
require(castor)
require(Biostrings)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)
require(foreach)
require(doParallel)

# Get tree metadata
meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
host_meta <- meta %>%
  distinct(host, host_species, host_genus, host_family, host_class, host_order)

like_threshold <- 2

# Create result directory
res_dir <- str_glue("results/ancestral_reconstruction_out/varying_thresholds/ancestral_reconstructions.like{like_threshold}")

tree_dir <- "data/trees/source_sink_mini_trees/without_outgroup.masked"
aln_list <- list.files("data/alignments/source_sink_mini_trees/without_outgroup.masked",
                       "\\.aln",
                       full.names = T)
aln_list <- aln_list[grepl("masked_to", aln_list)]

tree_paths <- list.files(tree_dir, ".treefile")

# For running only unfinished runs
# done <- gsub(".csv", "", list.files(str_glue("{res_dir}/temp_results")))
# 
# tree_paths <- deframe(tibble(tree_paths = tree_paths) %>%
#     separate(tree_paths, c("clique_name"), "\\.", remove = F) %>%
#     filter(!(clique_name %in% done)) %>%
#     select(tree_paths))

tree_paths <- tree_paths[!grepl("Circoviridae_2\\.", tree_paths)]

length(tree_paths)

# Count host jumps and plot tree
root_df <- fread("results/ancestral_reconstruction_out/final_source_sink_roots.csv")

foreach(tree_path = tree_paths,
        .packages = c("tidyverse", "ape", "castor",
                      "treeio", "phytools", "Biostrings",
                      "data.table", "foreach", "randomcoloR",
                      "ggtree")) %do% {
  # tree_path = tree_paths[1]
                        
  # Read tree
  tree <- ape::read.tree(str_glue("{tree_dir}/{tree_path}"))

  ## Parse node labels
  # tree$node.label <- str_split(tree$node.label, "/", simplify = T)[, 2]
  
  # Parse tip names
  tip_labels <- gsub("_R_", "", tree$tip.label)
  tree$tip.label <- tip_labels
  
  clique_name <- str_split(tree_path, "\\.")[[1]][1]
  
  # Add node names
  tree <- makeNodeLabel(tree)
  
  # Root tree
  root_acc <- deframe(root_df %>%
                       filter(cluster == clique_name) %>%
                       select(root_tip))
  
  rooted <- root(tree, root_acc, resolve.root = T)
  
  # Fix polytomies
  rooted$edge.length[rooted$edge.length == 0] <- 1e-11
  
  if (!is.binary(rooted)) {
    rooted <- fix.poly(rooted, type = "resolve", random = F)
  }
  
  # Get alignment length to scale branch lengths
  aln_path <- aln_list[grepl(paste0(clique_name, "\\."), aln_list)]
  aln_len <- unique(width(readDNAStringSet(aln_path)))
  
  # Match metadata to tips
  tips <- rooted$tip.label
  
  meta.match <- meta[match(tips, meta$accession), ] %>%
    # Recode missing hosts as unknown
    mutate(host = ifelse(host == ""|is.na(host), "Unknown", host))
  
  tip_states <- meta.match$host

  all(tips == meta.match$accession) # Test
  x <- setNames(tip_states, tips)
  
  # Ancestral reconstruction
  fitER <- ape::ace(x, rooted,
                    model="ER",
                    type="discrete",
                    method = "ML")
  
  ancstats <- as.data.frame(fitER$lik.anc)
  ancstats$node <- (1:rooted$Nnode) + Ntip(rooted)
  
  # Plot ancestral reconstructions
  dd <- data.frame(Accession = tips,
                   host = meta.match$host,
                   country = meta.match$country) %>%
    mutate(country = ifelse(country == "", NA, country)) %>%
    mutate(label = str_glue("{Accession}|{country}|{host}"))
  
  rooted_labelled <- rooted
  rooted_labelled$tip.label <- dd$label
  
  write.tree(rooted_labelled, str_glue("results/ancestral_reconstruction_out/varying_thresholds/ancestral_reconstructions.like2/labelled_trees/{clique_name}.labelled.tree"))
}


