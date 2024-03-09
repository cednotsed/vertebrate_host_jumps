rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

df <- fread("data/VIRION.v0.2.1_240922/Virion.csv.gz")

df %>%
  filter(HostGenus == "homo") %>%
  filter(ICTVRatified) %>%
  filter(DetectionMethod != "Not specified") %>%
  distinct(VirusGenus) %>%
  fwrite("../human_infecting_genera.txt")

  df %>% distinct(DetectionMethod)
  distinct(VirusGenus) %>% View()
