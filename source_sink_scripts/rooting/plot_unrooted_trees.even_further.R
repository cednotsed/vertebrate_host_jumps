rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
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

meta <- fread(str_glue("results/clique_classification_out/final_cluster_metadata.220723.new.csv"))
mash_dir <- "results/mash_out/source_sink_mini_trees/with_buffer_outgroup.even_further"
mash_paths <- list.files(mash_dir, ".tsv")
mash_paths <- mash_paths[grepl("reoviridae", mash_paths)]
# tree_paths <- tree_paths[grepl("Hepeviridae_17",
#                                tree_paths)]
# mash_path <- mash_paths[1]

morsels <- foreach(mash_path = mash_paths) %do% {
  print(mash_path)
  clique_name <- str_split(mash_path, "\\.")[[1]][1]
  
  # Get NJ tree
  cm <- read.csv(str_glue("{mash_dir}/{mash_path}"), 
                 sep = "\t", 
                 header = T,
                 row.names = 1,
                 stringsAsFactors = F)
  
  cm <- data.matrix(cm)
  tree <- bionj(cm) 
  
  write.tree(tree, 
             str_glue("data/trees/source_sink_mini_trees/with_buffer_outgroup.even_further/unrooted/{clique_name}.n{Ntip(tree)}.unrooted.nwk"))
  
  # Match metadata to tips
  tips <- tree$tip.label
  
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
  
  save_prefix <- gsub(".tsv", ".pdf", mash_path)
  ggsave(str_glue("results/source_sink_analysis/clique_trees/with_buffer_outgroup.even_further/{save_prefix}"),
         plot = p_linear,
         dpi = 100,
         width = 15,
         height = 15)
}


