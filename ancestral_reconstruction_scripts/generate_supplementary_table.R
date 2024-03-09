rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)

parsed <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  distinct(anc_name, anc_state, tip_state, clique_name) %>%
  mutate(involves_humans = ifelse(anc_state %in% c("Homo", "Homo sapiens") |
                                    tip_state %in% c("Homo", "Homo sapiens"),
                                  "Yes", "No")) %>%
  mutate(event_type = case_when(anc_state %in% c("Homo", "Homo sapiens") ~ "Anthroponotic",
                                 tip_state %in% c("Homo", "Homo sapiens") ~ "Zoonotic",
                                 TRUE ~ "Animal-to-animal"))

fwrite(parsed, "results/ancestral_reconstruction_out/host_jump_lists/distinct_jumps.csv")
  
parsed %>%
  summarise(n = n_distinct(clique_name))

parsed %>%
  filter(involves_humans == "Yes") %>%
  group_by(event_type) %>%
  summarise(n = n())

parsed2 <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  mutate(involves_humans = ifelse(anc_state %in% c("Homo", "Homo sapiens") |
                                    tip_state %in% c("Homo", "Homo sapiens"),
                                  "Yes", "No")) %>%
  mutate(event_type = case_when(anc_state %in% c("Homo", "Homo sapiens") ~ "Anthroponotic",
                                tip_state %in% c("Homo", "Homo sapiens") ~ "Zoonotic",
                                TRUE ~ "Animal-to-animal"))

fwrite(parsed2, "results/ancestral_reconstruction_out/host_jump_lists/all_jumps.csv")

nrow(parsed2)  
