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

like_threshold <- 10

# Create result directory
res_dir <- str_glue("results/ancestral_reconstruction_out/varying_thresholds/ancestral_reconstructions.same_host.like{like_threshold}")
dir.create(res_dir)
dir.create(str_glue("{res_dir}/temp_results"))

tree_dir <- "data/trees/source_sink_mini_trees/without_outgroup.masked"
aln_list <- list.files("data/alignments/source_sink_mini_trees/without_outgroup.masked",
                       "\\.aln",
                       full.names = T)
aln_list <- aln_list[grepl("masked_to", aln_list)]

tree_paths <- list.files(tree_dir, ".treefile")

# # For running only unfinished runs
# done <- gsub(".csv", "", list.files(str_glue("{res_dir}/temp_results")))
# 
# tree_paths <- deframe(tibble(tree_paths = tree_paths) %>%
#                         separate(tree_paths, c("clique_name"), "\\.", remove = F) %>%
#                         filter(!(clique_name %in% done)) %>%
#                         select(tree_paths))
tree_paths <- tree_paths[!grepl("Circoviridae_2\\.", tree_paths)]

length(tree_paths)

# Count host jumps and plot tree
cl <- makeCluster(12)
registerDoParallel(cl)

root_df <- fread("results/source_sink_analysis/final_source_sink_roots.csv")

foreach(tree_path = tree_paths,
        .packages = c("tidyverse", "ape", "castor",
                      "treeio", "phytools", "Biostrings",
                      "data.table", "foreach", "randomcoloR",
                      "ggtree")) %dopar% {
                        
  # Reinitialise variables
  like_threshold <- 10
  root_df <- fread("results/source_sink_analysis/final_source_sink_roots.csv")
  res_dir <- str_glue("results/ancestral_reconstruction_out/varying_thresholds/ancestral_reconstructions.same_host.like{like_threshold}")
  tree_dir <- "data/trees/source_sink_mini_trees/without_outgroup.masked"
  
  # Analysis function
  count_transitions <- function(aln_len, tr, tip_name, ancstats, x) {
    ## FOR TESTING ##
    # tr <- rooted
    # tip_name <- tip_filt[6]
    #################
    scaled_tr <- tr
    # scaled_tr$edge.length <- scaled_tr$edge.length * aln_len
    
    tr <- as_tibble(tr)
    tip_state <- x[[tip_name]]
    
    # Get ancestral depth of tip (min depth = 1)
    anc_depth <- nrow(ancestor(tr, tip_name))
    anc_depth
    
    if (anc_depth > 2) {
      traverse_limit <- floor(anc_depth)
      
      # Traverse to limit node
      node_name <- tip_name
      
      traverse_counter <- 1
      
      continue <- T
      
      count_temp <- tibble()
      
      while (traverse_counter <= traverse_limit & continue) {
        # Get parent and corresponding stats
        parent_name <- parent(tr, node_name)$node
        node_stats <- ancstats %>%
          filter(node == parent_name) %>%
          select(-node)
        
        if (sum(!is.na(node_stats)) >= 2) {
          # Check if likelihood is high
          max_like <- max(node_stats, na.rm = T)
          fold_diff <- max_like / max(node_stats[, node_stats != max_like])
          
          if (fold_diff > 2) {
            anc_state <- colnames(node_stats)[node_stats == max_like]
            
            # Skip if ancestral state is Unknown
            if (anc_state != "Unknown") {
              if (anc_state == tip_state) {
                patristic_dist <- get_pairwise_distances(scaled_tr, 
                                                         parent_name, 
                                                         tip_name, 
                                                         as_edge_counts = F)
                
                node_name <- scaled_tr$node.label[parent_name - Ntip(scaled_tr)]
                
                count_temp <- bind_rows(count_temp, 
                                        tibble(anc_name = node_name,
                                               tip_name = tip_name, 
                                               anc_state = anc_state,
                                               tip_state = tip_state,
                                               n_traverses = traverse_counter,
                                               total_depth = anc_depth,
                                               patristic_dist = patristic_dist))
              } else {
                # Stop when host state changes
                continue <- F
              }
            } else {
              # Stop if ancestral state is unknown
              continue <- F
            }
          } else { 
            # Stop when ancestral state is ambiguous
            continue <- F
          }
        } else {
          # Stop if no. ancestral states too low
          continue <- F
        }
        
        # Get next node
        node_name <- parent_name
        traverse_counter <- traverse_counter + 1
      }
    } else {
      # Return null if tip is too shallow
      return(NULL)
    }
    
    # Return if dataframe is not empty
    if (nrow(count_temp) > 0) {
      return(count_temp)
    } else {
      return(NULL)
    }
  }
  
  # tree_path = tree_paths[20]
  print(tree_path)
  tree <- ape::read.tree(str_glue("{tree_dir}/{tree_path}"))
  
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
  
  # Count host jumps
  # Only perform traversing for tips with known hosts
  tip_filt <- names(x[x != "Unknown"])
  
  crumbs <- foreach(tip_name = tip_filt) %do% {
    count_transitions(aln_len, rooted, tip_name, ancstats, x)
  }
  
  count_df <- bind_rows(crumbs) %>%
    mutate(clique_name = clique_name)
  
  # Write temp results
  fwrite(count_df,
         str_glue("{res_dir}/temp_results/{clique_name}.csv"))
  
  return(NULL)
}


