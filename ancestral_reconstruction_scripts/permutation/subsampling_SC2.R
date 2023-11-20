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

# Create directories
res_dir <- str_glue("results/ancestral_reconstruction_out/subsampling_analysis")
dir.create(res_dir)
dir.create(str_glue("{res_dir}/temp_results"))

tree_dir <- "data/trees/source_sink_mini_trees/without_outgroup.masked"
aln_list <- list.files("data/alignments/source_sink_mini_trees/without_outgroup.masked",
                       "\\.aln",
                       full.names = T)
aln_list <- aln_list[grepl("masked_to", aln_list)]

tree_paths <- list.files(tree_dir, ".treefile")

tree_path <- tree_paths[grepl("Coronaviridae_12", tree_paths)]

root_df <- fread("results/source_sink_analysis/final_source_sink_roots.csv")

# Read tree
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

# Count number of human tips
human_tips <- deframe(meta %>%
  filter(accession %in% rooted$tip.label) %>%
  filter(host_genus == "Homo") %>%
  select(accession))

length(human_tips)

# Get alignment length to scale branch lengths
aln_path <- aln_list[grepl(paste0(clique_name, "\\."), aln_list)]
aln_len <- unique(width(readDNAStringSet(aln_path)))

# Host checker function
is_diff_host <- function(host1, host2) {
  # host1 <- "Homo sapiens"
  # host2 <- "Alouatta"
  # host1 <- "Anas"
  # host2 <- "Anatidae"
  
  taxa <- bind_rows(host_meta[host_meta$host == host1, ],
                    host_meta[host_meta$host == host2, ]) %>%
    as_tibble()
  
  if(nrow(taxa) == 2) {
    ranks <- c("host_species", "host_genus", "host_family", "host_order", "host_class")
    
    # Ignore missing rank information
    to_ignore <- apply(taxa[, ranks], 2, function(x) {any(x == "")})
    rank_filt <- ranks[!to_ignore]
    
    # If any rank is different, hosts are different
    proceed <- any(taxa[1, rank_filt] != taxa[2, rank_filt])
  } else {
    proceed <- F
    print(str_glue("Error: {host1}|{host2}"))
  }
  
  return(proceed)
}

# Analysis function
count_transitions <- function(aln_len, tr, tip_name, ancstats, x) {
  ## FOR TESTING ##
  # tr <- rooted_filt
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
    
    while (traverse_counter <= traverse_limit) {
      # Get parent and corresponding stats
      parent_name <- parent(tr, node_name)$node
      node_stats <- ancstats %>%
        filter(node == parent_name) %>%
        select(-node)
      
      if (sum(!is.na(node_stats)) >= 2) {
        # Check if likelihood is high
        max_like <- max(node_stats, na.rm = T)
        
        # Break ties
        n_max <- length(node_stats[, node_stats == max_like])
        
        if(n_max == 1) {
          fold_diff <- max_like / max(node_stats[, node_stats != max_like])
          
          if (fold_diff > like_threshold) {
            anc_state <- colnames(node_stats)[node_stats == max_like]
            
            # Skip if ancestral state is Unknown
            if (anc_state != "Unknown") {
              if (anc_state != tip_state & 
                  is_diff_host(anc_state, tip_state)) {
                patristic_dist <- get_pairwise_distances(scaled_tr, 
                                                         parent_name, 
                                                         tip_name, 
                                                         as_edge_counts = F)
                
                node_name <- scaled_tr$node.label[parent_name - Ntip(scaled_tr)]
                
                return(tibble(anc_name = node_name,
                              tip_name = tip_name, 
                              anc_state = anc_state,
                              tip_state = tip_state,
                              n_traverses = traverse_counter,
                              total_depth = anc_depth,
                              patristic_dist = patristic_dist))
              }
            }
          }
        }
      }
      # Get next node
      node_name <- parent_name
      traverse_counter <- traverse_counter + 1
    }
  }
}

like_threshold <- 2

foreach(n = c(seq(0, 900, 100), 950, 990)) %do% {
  n = 950
  foreach(i = seq(10)) %do% {
    # Randomly remove human tips
    to_remove <- sample(human_tips, n, replace = F)
    
    rooted_filt <- drop.tip(rooted, to_remove)
    
    # Match metadata to tips
    tips <- rooted_filt$tip.label
  
    meta.match <- meta[match(tips, meta$accession), ] %>%
      # Recode missing hosts as unknown
      mutate(host = ifelse(host == ""|is.na(host), "Unknown", host))
    
    tip_states <- meta.match$host
  
    all(tips == meta.match$accession) # Test
    
    x <- setNames(tip_states, tips)
    
    # Ancestral reconstruction
    fitER <- ape::ace(x, rooted_filt,
                      model="ER",
                      type="discrete",
                      method = "ML")
    
    ancstats <- as.data.frame(fitER$lik.anc)
    ancstats$node <- (1:rooted_filt$Nnode) + Ntip(rooted_filt)
    
    # Count host jumps
    # Only perform traversing for tips with known hosts
    tip_filt <- names(x[x != "Unknown"])
    
    crumbs <- foreach(tip_name = tip_filt) %do% {
      count_transitions(aln_len, rooted_filt, tip_name, ancstats, x)
    }
    
    count_df <- bind_rows(crumbs) %>% 
      mutate(n_human = 1000 - n,
             iter = i) %>%
      distinct(anc_name, tip_state, .keep_all = T)
    
    # Write temp results
    fwrite(count_df,
           str_glue("{res_dir}/temp_results/Coronaviridae_12.human{1000-n}.iter{i}.csv"))
  }
}



