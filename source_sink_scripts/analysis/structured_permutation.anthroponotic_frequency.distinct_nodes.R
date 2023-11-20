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

genome_meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  dplyr::rename(clique_name = cluster)

jump_df <- fread("results/ancestral_reconstruction_out/host_jump_lists/like2.diff_hosts.genus_counts.all_jumps.V2.csv") %>%
  filter(is_jump) %>%
  distinct(clique_name, anc_name, tip_state, .keep_all = T)

jump_df %>% View()
zoo_df <- jump_df %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  mutate(host = ifelse(event_type == "Zoonotic", anc_state, tip_state))

## Permutation test
clique_list <- unique(jump_df$clique_name)

set.seed(69)
cl <- makeCluster(12)
registerDoParallel(cl)

perm_morsels <- foreach(i = seq(500), .packages = c("foreach", "tidyverse")) %dopar% {
  print(i)
  
  # Recode within cliques
  temp_morsels <- foreach(clique = clique_list) %do% {
    temp_df <- jump_df %>%
      filter(clique_name == clique)
    
    # Encode hosts as numeric
    host_list <- unique(c(temp_df$anc_state, temp_df$tip_state))
    host_number <- as.numeric(as.factor(host_list))
    
    # Generate random numeric to host mappings
    host_to_number <- setNames(host_number, host_list)
    number_to_host_list <- setNames(sample(host_list, replace = F), host_number)
    
    # Recode shuffled hosts
    recoded_temp <- temp_df %>%
      mutate(anc_recoded = recode(anc_state, !!!host_to_number),
             tip_recoded = recode(tip_state, !!!host_to_number)) %>%
      mutate(anc_recoded = recode(anc_recoded, !!!number_to_host_list),
             tip_recoded = recode(tip_recoded, !!!number_to_host_list))
    
    return(recoded_temp)
  }
  
  perm_df <- bind_rows(temp_morsels) %>%
    select(anc_recoded, tip_recoded, tip_name) %>%
    left_join(genome_meta %>% distinct(anc_recoded = host, anc_genus = host_genus)) %>%
    left_join(genome_meta %>% distinct(tip_recoded = host, tip_genus = host_genus)) %>%
    filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
    mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic"))
  
  perm_stats <- perm_df %>%
    group_by(event_type) %>%
    summarise(n_jumps = n()) %>%
    mutate(iter = i,
           total_jumps = nrow(perm_df))
  
  return(perm_stats)
}

stopCluster(cl)

obs_df <- zoo_df %>%
  summarise(obs_prop = sum(event_type == "Anthroponotic") / sum(event_type == "Zoonotic"))

plot_df <- bind_rows(perm_morsels) %>%
  pivot_wider(id_cols = iter, names_from = "event_type", values_from = "n_jumps") %>%
  mutate(prop = Anthroponotic / Zoonotic)

p_val <- sum(plot_df$prop > obs_df$obs_prop) / 1000

plot_df %>%
  ggplot(aes(x = log(prop))) +
  geom_histogram() +
  geom_vline(xintercept = log(obs_df$obs_prop),
             color = "red",
             lty = "dashed") +
  geom_text(x = 0.9, y = 40, label = str_glue("p={p_val}"),
            color = "red")

ggsave("results/source_sink_analysis/structured_permutation_test.png")
