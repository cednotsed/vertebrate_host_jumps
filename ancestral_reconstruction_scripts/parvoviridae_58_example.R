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
tree_dir <- "data/trees/source_sink_mini_trees/without_outgroup.masked"
aln_list <- list.files("data/alignments/source_sink_mini_trees/without_outgroup.masked",
                       "\\.aln",
                       full.names = T)
aln_list <- aln_list[grepl("masked_to", aln_list)]

root_df <- fread("results/source_sink_analysis/final_source_sink_roots.csv")

tree_paths <- list.files(tree_dir, ".treefile")
tree_path = tree_paths[grepl("Parvoviridae_58", tree_paths)]

# Read tree
tree <- ape::read.tree(str_glue("{tree_dir}/{tree_path}"))

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
  # tr <- rooted
  # tip_name <- tip_filt[3]
  #################
  scaled_tr <- tr
  # scaled_tr$edge.length <- scaled_tr$edge.length * aln_len
  
  tr <- as_tibble(tr)
  tip_state <- x[[tip_name]]
  
  # Get ancestral depth of tip (min depth = 1)
  anc_depth <- nrow(ancestor(tr, tip_name))
  
  
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
                
                node_name <- scaled_tr$node.label[parent_name - Ntip(tree)]
                
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
  mutate(country = ifelse(country == "", NA, country))

cols <- setNames(distinctColorPalette(length(unique(x))), sort(unique(x)))
pies <- nodepie(ancstats, 
                cols = 1:(ncol(ancstats) - 1),
                color = cols)

ggtree(rooted, 
       # branch.length = "none",
       ladderize = T,
       size = 0.001,
       color = "darkslategrey") %<+% dd + 
  geom_tippoint(aes(color = host),
                size = 3) + 
  scale_color_manual(values = cols) +
  theme(legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),) +
  geom_inset(pies, 
             width=0.03, 
             height=0.03) +
  geom_nodelab(size = 5, vjust = 1) +
  labs(color = "Host") +
  geom_tiplab()
  # xlim(0, 20)

# Count host jumps
# Only perform traversing for tips with known hosts
tip_filt <- names(x[x != "Unknown"])

crumbs <- foreach(tip_name = tip_filt) %do% {
  count_transitions(aln_len, rooted, tip_name, ancstats, x)
}

count_df <- bind_rows(crumbs) %>%
  mutate(clique_name = clique_name)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(clique_name == "Parvoviridae_58")

dnds_df <- fread("results/dnds_out/all_jumps.dnds.diff_hosts.genus_counts.csv") %>%
  filter(clique_name == "Parvoviridae_58")

dnds_df %>%
  mutate(kaks = ka / ks) %>% View()
  group_by(is_jump) %>%
  summarise(min_kaks = min(kaks)) %>% View()
jump_df %>%
  View()
  group_by(is_jump) %>%
  summarise(min_dist =)
meta %>%
  filter(cluster == "Parvoviridae_58")
# fwrite(jump_df %>%
#          select(tip_name, anc_name, tip_state, anc_state, is_jump),
#        "results/ancestral_reconstruction_out/example/parvoviridae_58.csv")


