rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Hmisc)

host_meta <- fread("data/metadata/parsed_host_metadata.csv") %>%
  rename(parsed_host = host)

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))
# mutate(host_species = replace_na(host_species, ""),
#        host_genus = replace_na(host_genus, ""))

# Parse host meta
host_list <- meta %>%
  distinct(host)

# Parse host with more than two words
word_count <- host_list %>%
  mutate(parsed_host = host) %>%
  mutate(parsed_host = gsub("  ", " ", parsed_host)) %>%
  mutate(n_words = str_count(parsed_host, "\\w+"))

short_names <- word_count %>%
  filter(n_words <= 2) %>%
  mutate(parsed_host = gsub(" sp\\.| subsp\\.| subgen\\.", "", parsed_host))

long_names <- word_count %>%
  filter(n_words > 2) %>%
  separate(parsed_host, c("host1", "host2"), " ", remove = F) %>%
  mutate(parsed_host = str_glue("{host1} {host2}")) %>%
  mutate(parsed_host = gsub(" sp\\.| subsp\\.| subgen\\.", "", parsed_host))

combined <- bind_rows(long_names, short_names) %>%
  mutate(parsed_host = ifelse(!(parsed_host %in% host_meta$parsed_host),
                              word(parsed_host, 1),
                              parsed_host)) %>%
  left_join(host_meta) %>%
  select(-host1, -host2, -n_words)

# Merge host meta with main metadata df
parsed_meta <- meta %>%
  left_join(combined) %>%
  mutate(host_species = replace_na(host_species, ""),
         host_genus = replace_na(host_genus, ""),
         host_family = replace_na(host_family, ""),
         host_order = replace_na(host_order, ""),
         host_class = replace_na(host_class, ""))

fwrite(parsed_meta, "data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.parsed_host.csv")
