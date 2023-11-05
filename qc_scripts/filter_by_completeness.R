rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(data.table)
require(tidyverse)
require(foreach)
require(Biostrings)

# Merge quality summary
qual_df <- fread("results/checkv_out/quality_summary.tsv") %>%
  # bind_rows(fread("results/checkv_out/missing/quality_summary.tsv")) %>%
  dplyr::rename(accession = contig_id) %>%
  distinct()

meta <- fread("data/metadata/all_viruses.220723.filt.csv")

merged <- meta %>%
  left_join(qual_df)

# Check for missing genomes
merged %>% 
  filter(is.na(contig_length)) %>%
  nrow()

# Check for duplicates
sum(duplicated(merged$accession))
# # Visualise contamination and completeness
# merged %>%
#   ggplot(aes(x = completeness)) +
#   geom_histogram(bins = 100)
# 
# merged %>%
#   ggplot(aes(x = contamination)) +
#   geom_histogram(bins = 100)

# Filter genomes
merged_filt <- merged %>%
  filter(contamination < 5) %>%
  filter((completeness > 95 & !is_segmented) | is_segmented) %>%
  # Parse collection dates
  separate(collection_date, into = c("Y", "M", "D"), sep = "-", remove = F) %>%
  filter(Y != "") %>%
  filter(!is.na(Y)) %>%
  mutate(D = ifelse(!is.na(M) & is.na(D), "15", D)) %>%
  mutate(M = ifelse(is.na(M) & is.na(D), "06", M)) %>%
  mutate(D = ifelse(is.na(D), "01", D)) %>%
  mutate(imputed_date = as.Date(str_glue("{Y}-{M}-{D}"), "%Y-%m-%d")) %>%
  select(-Y, -M, -D) %>%
  # Filter columns
  select(accession, genbank_title, family, 
         genus, species, host, 
         isolation_source, imputed_date, country,
         genome_length = length, molecule_type, is_segmented,
         is_circular)

# Filter by genome count (AGAIN!)
genome_counts <- merged_filt %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  filter(n >= 100)

merged_filt2 <- merged_filt %>%
  filter(family %in% genome_counts$family)

# Check final family counts
plot_df <- merged_filt2 %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

plot_df %>%
  mutate(family = factor(family, unique(plot_df$family))) %>%
  ggplot(aes(x = family, y = n, fill = family)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Viral family", y = "No. genomes")

ggsave("results/qc_out/family_genome_counts.pdf", width = 8, height = 5)

# host_plts <- merged_filt2 %>%
#   group_by(family, host) %>%
#   summarise(n = n()) %>%
#   ggplot(aes(x = host, y = n)) +
#   facet_wrap(~family) +
#   geom_bar(stat = "identity")

# Write filtered metadata
fwrite(merged_filt2, "data/metadata/all_viruses.220723.filt.QCed.csv")

# Read genomes
fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.fna")

# Parse genome names
names(fna) <- str_split(names(fna), " ", simplify = T)[, 1]

fna_filt <- fna[names(fna) %in% merged_filt2$accession]
nrow(merged_filt2) == length(fna_filt)

writeXStringSet(fna_filt, "data/genomes/all_viruses.220723.filt.formatted.QCed.fna")

sum(duplicated(names(fna)))
