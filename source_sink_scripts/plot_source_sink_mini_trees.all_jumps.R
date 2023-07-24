rm(list = ls())
setwd("C:/git_repos/viral_sharing/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(treeio)
require(castor)
require(foreach)
require(phangorn)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)

meta <- fread(str_glue("results/species_classification_out/final_cluster_metadata.2.140423.csv")) %>%
  separate(collection_date, c("Y"), "-") %>%
  mutate(host_genus = case_when(host_genus == "Hesperosciurus" ~ "Sciurus",
                                host_genus == "Neosciurus" ~ "Sciurus",
                                host_genus == "Pholidota" ~ "Manis",
                                TRUE ~ host_genus))

tree_dir <- "data/trees/source_sink_mini_trees/with_outgroup_buffer.all_jumps/"
tree_paths <- list.files(tree_dir, ".treefile")
# tree_paths <- tree_paths[grepl("Hepeviridae_17",
#                                tree_paths)]
tree_paths

tree_path <- tree_paths[1]
tree_path

# Get external outgrps
root_df <- fread("results/source_sink_analysis/source_sink_results.all_jumps.curated.csv") %>%
  filter(root_resolved == "T")

root_morsels <- foreach(tree_path = tree_paths) %do% {
  tree <- read.tree(str_glue("{tree_dir}{tree_path}"))
  clique_name <- str_split(tree_path, "\\.")[[1]][1]
  
  # Parse tips
  # tree <- tree %>%
  #     as_tibble() %>% 
  #     mutate(label = gsub("_R_", "", label)) %>%
  #     as.treedata()
  
  tree$tip.label <- gsub("_R_", "", tree$tip.label)
  
  # Root and drop external cliques
  if (clique_name %in% root_df$cluster) {
    root_names <- deframe(root_df %>%
      filter(cluster == clique_name) %>%
        select(root_cliques))

    root_names <- str_split(root_names, ";", simplify = T)[1, ]

    root_acc_filt <- deframe(meta %>%
                              filter(cluster %in% root_names) %>%
                              filter(accession %in% as.phylo(tree)$tip.label) %>%
                              select(accession))
    to_remove <- deframe(meta %>%
                           filter(accession %in% as.phylo(tree)$tip.label) %>%
                           filter(cluster != clique_name) %>%
                           select(accession))

    if(clique_name == "Rhabdoviridae_12") {
      exceptions <- c("JQ685907.1", "JQ685941.1")
      tree <- drop.tip(tree, exceptions)
    }

    tree <- root(tree, root_acc_filt)
    tree <- drop.tip(tree, to_remove)
    
    write.tree(as.phylo(tree),
               str_glue("data/trees/source_sink_mini_trees/with_outgroup_buffer.all_jumps/rooted/{clique_name}.rooted.treefile"))
  }
    
  # Match metadata to tips
  tips <- as.phylo(tree)$tip.label

  meta.match <- meta[match(tips, meta$accession), ]
  
  all(tips == meta.match$accession)
  
  # Shift Homo to top level
  host_levels <- unique(meta.match$host_genus)
  host_levels <- c("Homo", host_levels[host_levels != "Homo"])
  
  # Shift clique of interest to top level
  clique_levels <- unique(meta.match$cluster) 
  clique_levels <- c(clique_name, clique_levels[clique_levels != clique_name])
  
  dd <- data.frame(Accession = tips,
                   host = meta.match$host_genus,
                   country = meta.match$country,
                   collection_date = meta.match$Y,
                   viral_clique = meta.match$cluster) %>%
    mutate(country = ifelse(country == "", NA, country)) %>%
    mutate(of_interest = ifelse(viral_clique == clique_name, T, NA)) %>%
    mutate(host = factor(host, host_levels)) %>%
    mutate(viral_clique = factor(viral_clique, clique_levels))
  
  # Plot tree
  set.seed(66)
  random_pal <- distinctColorPalette(length(unique(dd$host)))
  random_pal2 <- distinctColorPalette(length(unique(dd$host)))
  
  # Make human tips black
  random_pal2[1] <- "black"
  
  # random_pal[length(random_pal) - 1] <- "darkolivegreen4"
  
  p_linear <- ggtree(tree,
                     size = 0.001,
                     # layout = "unrooted",
                     branch.length = "none",
                     color = "darkslategrey") %<+% dd +
    geom_tippoint(aes(hjust = 0.5, 
                      shape = viral_clique,
                      color = host), alpha = 1, size = 1.5) + 
    geom_tiplab(size = 1) + 
    # geom_tiplab(size = 1) +
    scale_color_manual(values = random_pal2) +
    geom_fruit(geom = geom_tile,
               aes(fill = of_interest),
               offset = 0.2,
               width = 0.01) +
    scale_fill_discrete(na.translate = F) +
    labs(fill = "Clique of interest") +
    # labs(fill = "Country") +
    new_scale_fill() +
    # geom_fruit(geom = geom_tile,
    #            aes(fill = as.numeric(collection_date)),
    #            color = "black",
    #            offset = 0.2001,
    #            width = 0.02) +
    # scale_fill_gradient(na.value = "grey") +
    labs(color = "Viral clique", title = clique_name)
  
  if (length(tips) < 30) {
    # p_linear <- p_linear + geom_text(aes(x = branch, label=UFboot), color = "black",
    #                                  hjust = 0.5, vjust = -0.5) 
    p_linear <- p_linear + geom_nodelab()
  }
  
  save_prefix <- gsub(".treefile", ".pdf", tree_path)
  ggsave(str_glue("results/source_sink_analysis/clique_trees/with_outgroup_buffer.all_jumps/{save_prefix}"),
         plot = p_linear,
         dpi = 100,
         width = 15,
         height = 15)
}


