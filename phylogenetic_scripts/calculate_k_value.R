rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(stringr)

meta <- fread("data/metadata/all_viruses.220723.filt.QCed.csv")

family_list <- deframe(meta %>% distinct(family))
family_list

# Probability of observing a k-mer by chance
q <- 0.01

# no. of letters in alphabet
n_alphabet <- 4

k_morsels <- foreach(family_name = family_list) %do% {
  print(family_name)
  
  # Calculate median length of genomes per family
  l <- deframe(meta %>%
    filter(family == family_name) %>%
    summarise(median_length = median(genome_length)))
  
  # Calculate k value for mash
  k <- ceiling(log(l * (1 - q)/ q, base = n_alphabet))
  
  tibble(family = family_name, k = k)
}

bind_rows(k_morsels) %>%
  arrange(desc(k))


