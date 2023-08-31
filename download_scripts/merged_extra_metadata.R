rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)

acc <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.full_acc_list.csv")

nrow(acc)

meta1 <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.gt1000nt.220723.csv") %>%
  select(-Sequence_Type)

meta2 <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.lt1000nt.220723.csv") %>%
  select(colnames(meta1))

merged <- meta1 %>%
  bind_rows(meta2)

fwrite(merged, "data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.csv")

