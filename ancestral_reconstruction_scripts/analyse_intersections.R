rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)
require(Hmisc)
require(doParallel)
require(UpSetR)

threshold_list <- setNames(c(2, 5, 10), c("two", "five", "ten"))

df_morsels <- foreach(idx = seq(length(threshold_list))) %do% {
  t <- threshold_list[[idx]]
  
  merged <- fread(str_glue("results/ancestral_reconstruction_out/host_jump_lists/like{t}.diff_hosts.genus_counts.all_jumps.V2.csv")) %>%
    filter(is_jump) %>%
    select(clique_name, anc_name, anc_state, tip_state, tip_name) %>%
    distinct(clique_name, anc_state, tip_state, anc_name) # Get distinct host jump events
  
  # Get distinct host pairs
  jump_df <- merged %>% 
    distinct(clique_name, anc_state, tip_state)
  
  # For each host pair, add integer column indicating jump count
  merged <- foreach(i = seq(nrow(jump_df)), .combine = "bind_rows") %do% {
    row <- jump_df[i, ]
    temp <- merged %>%
      filter(anc_state == row$anc_state, 
             tip_state == row$tip_state,
             clique_name == row$clique_name)
    temp %>%
      select(-anc_name) %>%
      mutate(jump_index = seq(nrow(temp)))
  }
  
  # Add column for intersection
  return(merged %>% add_column(!!names(threshold_list)[idx] := 1))
}

# Get presence absence columns
all <- df_morsels[[1]] %>%
  full_join(df_morsels[[2]]) %>%
  full_join(df_morsels[[3]]) %>%
  mutate(across(all_of(c("two", "five", "ten")), 
                ~ replace_na(.x, 0))) %>%
  select(-clique_name, -anc_state, -tip_state, -jump_index) %>%
  as.data.frame()

# Plot upset
pdf("results/ancestral_reconstruction_out/upset_plot.pdf",
    height = 5, width = 4)
upset(all, sets = c("two", "five", "ten"), 
      order.by = "freq")
dev.off()


# Intersection between all
# intersect(intersect(df_morsels[[1]] %>% 
#                       select(-one), 
#                     df_morsels[[2]] %>%
#                       select(-five)), 
#           df_morsels[[3]] %>%
#             select(-ten))

# Intersection between one and five
# intersect(df_morsels[[1]] %>% 
#             select(-one), 
#           df_morsels[[2]] %>%
#             select(-five))

# Plot anthro versus zoonotic
zoo_morsels <- foreach(idx = seq(length(threshold_list)), .combine = "bind_rows") %do% {
  t_name <- names(threshold_list)[idx]
  t <- threshold_list[[idx]]
  merged_df <- fread(str_glue("results/ancestral_reconstruction_out/host_jump_lists/like{t}.diff_hosts.genus_counts.all_jumps.V2.csv")) %>%
    filter(is_jump) %>%
    distinct(clique_name, anc_state, tip_state, anc_name, .keep_all = T) %>% # Get distinct host jump events
    mutate(event_type = case_when(anc_genus == "Homo" ~ "Anthroponotic",
                                 tip_genus == "Homo"~ "Zoonotic")) %>% 
    filter(!is.na(event_type)) %>%
    group_by(event_type) %>%
    summarise(n = n()) %>%
    mutate(threshold = t_name)
}

zoo_morsels %>%
  mutate(threshold = factor(threshold, c("two", "five", "ten"))) %>%
  ggplot(aes(x = event_type, y = n, fill = event_type)) +
  geom_bar(stat = "identity",
           color = "black") +
  facet_grid(rows = vars(threshold)) +
  geom_text(aes(label = n, vjust = 0)) +
  theme_classic() +
  labs(x = "Event type", y = "No. jumps") +
  theme(legend.position = "none")

ggsave("results/ancestral_reconstruction_out/varying_thresholds_anthro_zoo_frequency.pdf",
       dpi = 600,
       width = 4,
       height = 6)
