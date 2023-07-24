rm(list = ls())
setwd("C:/git_repos/viral_sharing/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)

cm <- read.csv(str_glue("results/phylogenetic_out/mash_out/all_viruses.genbank_complete_excllabhost_exclvaccine.191222.unique.tsv"), sep = "\t", 
               header = T,
               row.names = 1,
               stringsAsFactors = F)

cm <- data.matrix(cm)
tree <- nj(cm)

# Get tree metadata
meta <- read.csv(str_glue("data/metadata/all_viruses.genbank_complete_excllabhost_exclvaccine.191222.unique.csv"), 
                 check.names = F, 
                 stringsAsFactors = F)

# Save annotated tree
new_labels <- tibble(accession = tree$tip.label) %>%
  left_join(meta) %>%
  mutate(new_labels = str_glue("{accession}|{genus}"))
new_labels <- new_labels$new_labels
annot_tree <- tree
annot_tree$tip.label <- new_labels
write.tree(annot_tree, "data/trees/coronaviridae_NJ_mash.annotated.tree")

# Save genbank_title annotated tree
new_labels <- tibble(accession = tree$tip.label) %>%
  left_join(meta) %>%
  mutate(new_labels = str_glue("{accession}|{genus}|{genbank_title}"))
new_labels <- new_labels$new_labels
annot_tree <- rooted
annot_tree$tip.label <- new_labels
write.tree(annot_tree, "data/trees/coronaviridae_NJ_mash.genbank_title_annotated.tree")

# No. of genera for annotation
genera_counts <- table(meta$host_genus)
genera_filt <- names(genera_counts[genera_counts > 5])
genera_filt <- genera_filt[genera_filt != ""]
genera_filt <- genera_filt[genera_filt != "Chiroptera"]

# Match metadata to tips
meta.match <- meta[match(rooted$tip.label, meta$accession), ] %>%
  separate(accession, into = c("isolate_name"), sep = "_", remove = F) %>%
  mutate(annotation = case_when(grepl("Rhino", host) ~ "Rhinolophus",
                                grepl("Camelus", host) ~ "Camelus",
                                TRUE ~ as.character(NA))) %>%
  mutate(genus = ifelse(genus == "", NA, genus),
         new_isolate = ifelse(species == "novel" |accession == "MW719567.1", 
                              isolate_name, NA),
         is_new_isolate = ifelse(species == "novel", T, NA),
         host_genus_filt = ifelse(host_genus %in% genera_filt, 
                                  host_genus, NA),
         common_name = ifelse(common_name == "", NA, common_name),
         tip_labels = ifelse(species == "novel", isolate_name, genbank_title),
         accession_labels = ifelse(species == "novel", "novel", accession))

all(rooted$tip.label == meta.match$accession)

dd <- data.frame(Accession = meta.match$accession,
                 new_isolate = meta.match$new_isolate,
                 is_new_isolate = meta.match$is_new_isolate,
                 host = meta.match$host,
                 host_genus_filt = meta.match$host_genus_filt,
                 common_name = meta.match$common_name,
                 annotation = meta.match$annotation,
                 genus = meta.match$genus,
                 species = meta.match$genbank_title,
                 tip_labels = meta.match$tip_labels,
                 accession_labels = meta.match$accession_labels)

# Plot tree
col_pal <- c("darkcyan", "seagreen4", "darkgoldenrod3", 
                       "#2e4057", "#d1495b", "cornflowerblue", 
                       "maroon4", "darkorchid4", "slateblue4",
                       "black", "darkgrey")
                       
set.seed(66)
random_pal <- distinctColorPalette(length(unique(dd$host_genus_filt)))
# random_pal[length(random_pal) - 1] <- "darkolivegreen4"

p_circle <- ggtree(rooted, 
                   layout= "fan",
                   size = 0.1,
                   branch.length = "none",
                   color = "darkslategrey",
                   options(ignore.negative.edge = TRUE)) %<+% dd +
  geom_tippoint(aes(hjust = 0.5, color = genus), alpha = 1, size = 1) +
  scale_color_discrete(na.translate = F) +
  labs(color = "Viral genus", fill = "Host genus") +
  new_scale_color() +
  geom_tippoint(aes(color = is_new_isolate), alpha = 1, size = 2) +
  scale_color_manual(values = c("black"), na.translate = F, guide = "none") +
  geom_fruit(geom = geom_tile, 
             aes(fill = host_genus_filt),
             offset = 0.13,
             width = 10) +
  scale_fill_manual(values = random_pal, 
                    na.translate = F) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, 
             aes(fill = common_name),
             offset = 0.2,
             width = 10) +
  scale_fill_manual(values = col_pal, 
                    na.translate = F) +
  labs(fill = "Host")
p_circle
ggsave("results/phylogenetic_out/cov_NJ_mash.n2127.circle.pdf",
       plot = p_circle,
       dpi = 300,
       width = 8,
       height = 8)

# p_linear <- ggtree(rooted, 
#                    size = 0.001,
#                    branch.length = "none",
#                    color = "darkslategrey",
#                    options(ignore.negative.edge = TRUE)) %<+% dd +
#   geom_tippoint(aes(hjust = 0.5, color = genus), alpha = 1, size = 1) +
#   scale_color_discrete(na.translate = F) +
#   geom_tiplab(aes(label = tip_labels), size = 2) +
#   labs(color = "Viral genus", fill = "Host genus") +
#   new_scale_color() +
#   geom_tippoint(aes(color = is_new_isolate), alpha = 1, size = 1) +
#   scale_color_manual(values = c("black"), na.translate = F, guide = "none") +
#   geom_fruit(geom = geom_tile, 
#              aes(fill = host_genus_filt),
#              color = "grey",
#              offset = 0.15,
#              width = 10) +
#   scale_fill_manual(values = random_pal, 
#                     na.translate = F) +
#   new_scale_fill() +
#   geom_fruit(geom = geom_tile, 
#              aes(fill = common_name),
#              color = "grey",
#              offset = 0.25,
#              width = 10) +
#   scale_fill_manual(values = col_pal, 
#                     na.translate = F) +
#   labs(fill = "Host")
# 
# ggsave("results/phylogenetic_out/cov_NJ_mash.n2127.linear.pdf",
#        plot = p_linear,
#        dpi = 300,
#        width = 20,
#        height = 20)


