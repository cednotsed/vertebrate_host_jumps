rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(Hmisc)
require(ggpubr)
require(randomcoloR)
require(see)
require(ggrepel)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv")
host_meta <- fread("data/metadata/parsed_host_metadata.csv")
genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  mutate(host = ifelse(grepl("sp\\.", host), 
                       gsub(" sp\\.", "", host), 
                       host)) %>%
  mutate(host = ifelse(host %in% c("Homo", "Homo sapiens"), 
                       "Homo sapiens", 
                       host)) %>%
  filter(host != "") %>%
  dplyr::rename(clique_name = cluster)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/putative_host_jumps.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name)

zoo_df <- jump_df %>%
  mutate(anc_state = ifelse(anc_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", anc_state),
         tip_state = ifelse(tip_state %in% c("Homo", "Homo sapiens"), 
                            "Homo sapiens", tip_state)) %>%
  filter(anc_state == "Homo sapiens" | tip_state == "Homo sapiens") %>%
  mutate(event_type = ifelse(anc_state != "Homo sapiens", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state)) %>%
  inner_join(host_meta)

iter_df <- zoo_df %>% distinct(clique_name, anc_state, tip_state, event_type)

zoo_df %>%
  distinct(clique_name, anc_name, tip_state, .keep_all = T)

zoo_filt <- foreach(i = seq(nrow(iter_df)), .combine = "bind_rows") %do% {
  # i = 1
  row <- iter_df[i, ]
  anc <- row$anc_state
  tip <- row$tip_state
  clique <- row$clique_name
  
  temp_filt <- zoo_df %>%
    filter(clique_name == clique,
           anc_state == anc,
           tip_state == tip) %>%
    distinct(anc_name, .keep_all = T)
  
  return(temp_filt)
}

# Bootstrap zoonotic proportion
obs_df <- zoo_filt %>%
  group_by(event_type) %>%
  summarise(n_jumps = n()) %>%
  mutate(prop = n_jumps / sum(n_jumps))

clique_list <- unique(zoo_filt$clique_name)

set.seed(66)
cl <- makeCluster(12)
registerDoParallel(cl)

perm_morsels <- foreach(i = seq(1000), .packages = c("foreach", "tidyverse")) %dopar% {
  print(i)
  temp_morsels <- foreach(clique = clique_list) %do% {
    temp_df <- zoo_df %>%
      filter(clique_name == clique)
    
    # Encode hosts as numeric
    host_list <- unique(c(temp_df$anc_state, temp_df$tip_state))
    host_number <- as.numeric(as.factor(host_list))
    
    # Recode shuffled hosts
    code_list <- setNames(host_number, host_list)
    recode_list <- setNames(sample(host_list, replace = F), host_number)
    
    recoded_temp <- temp_df %>%
      mutate(anc_recoded = recode(anc_state, !!!code_list),
             tip_recoded = recode(tip_state, !!!code_list)) %>%
      mutate(anc_recoded = recode(anc_recoded, !!!recode_list),
             tip_recoded = recode(tip_recoded, !!!recode_list))
      # select(anc_state, anc_recoded, tip_state, tip_recoded) %>% View()
    recoded_temp
  }
  
  temp_df <- bind_rows(temp_morsels) %>%
    filter(anc_recoded == "Homo sapiens" | tip_recoded == "Homo sapiens") %>%
    mutate(event_type = ifelse(anc_recoded != "Homo sapiens", "Zoonotic", "Anthroponotic")) %>%
    mutate(host = ifelse(event_type == "Zoonotic", anc_recoded, tip_recoded)) %>%
    distinct(clique_name, anc_name, tip_recoded, .keep_all = T) %>%
    group_by(event_type) %>%
    summarise(n_jumps = n()) %>%
    mutate(iter = i)
  
  return(temp_df)
}

stopCluster(cl)


obs_df <- zoo_filt %>%
  summarise(obs_prop = sum(event_type == "Anthroponotic") / sum(event_type == "Zoonotic"))

plot_df <- bind_rows(perm_morsels) %>%
  pivot_wider(id_cols = iter, names_from = "event_type", values_from = "n_jumps") %>%
  mutate(prop = Anthroponotic / Zoonotic)

plot_df %>%
  ggplot(aes(x = log(prop))) +
  geom_histogram() +
  geom_vline(xintercept = log(obs_df$obs_prop))

plot_df %>%
  summarise(pval = sum(prop >= obs_df$obs_prop) / 1000)
