require(tidyverse)
require(data.table)
require(Biostrings)

meta <- fread("data/metadata/all_viruses.220723.filt.csv")

fna <- readDNAStringSet("data/genomes/all_viruses.220723.filt.formatted.fna")

res_df <- tibble(accession = names(fna), 
       parsed_length = width(fna)) %>%
  separate(accession, c("accession"), "\\ ") %>%
  left_join(meta)

res_df %>%
  filter(parsed_length == length) %>%
  View()

all(res_df$parsed_length == res_df$length)
