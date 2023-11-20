rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(data.table)
require(tidyverse)
require(ape)
require(randomcoloR)
require(castor)
require(foreach)
require(adephylo)
require(doParallel)
require(easyPubMed)
require(Hmisc)

contree <- read.tree("data/supertree/olival_cytb_supertree.tree")
contree$tip.label <- gsub("_", " ", contree$tip.label)

df <- fread("data/VIRION.v0.2.1_240922/Virion.csv.gz") %>%
  filter(VirusNCBIResolved)

df_filt <- df %>%
  filter(DetectionMethod %in% c("PCR/Sequencing", "Isolation/Observation")) %>%
  filter(HostNCBIResolved) %>%
  filter(ICTVRatified)

host_count <- df_filt %>%
  group_by(Virus) %>%
  summarise(n_hosts = n_distinct(Host)) %>%
  arrange(desc(n_hosts)) %>%
  filter(n_hosts > 1)

virus_list <- host_count$Virus

morsels <- foreach(virus = virus_list) %do% {
  # taxid = virus_taxid_list[2]
  host_list <- deframe(df_filt %>%
      filter(Virus == virus) %>%
      distinct(Host) %>%
      arrange(Host))
  
  edgelist <- as_tibble(t(combn(host_list, 2, simplify = T))) %>%
    mutate(viral_species = virus)
  
  return(edgelist)
}

parsed <- bind_rows(morsels) %>%
  rename(host1 = V1, host2 = V2) %>%
  mutate(host1 = capitalize(host1), 
         host2 = capitalize(host2)) %>%
  filter(host1 %in% contree$tip.label & host2 %in% contree$tip.label)

parsed
# # Check reverse
# foreach(i = seq(nrow(parsed))) %do% {
#   row <- parsed[i, ]
#   n <- parsed %>%
#     filter(host2 == row$host1, host1 == row$host2) %>%
#     nrow()
#   print(n)
# }

sharing_df <- parsed %>%
  group_by(host1, host2) %>%
  summarise(n = n_distinct(viral_species)) %>%
  arrange(desc(n)) %>%
  ungroup()

cl <- makeCluster(12)
registerDoParallel(cl)

dist_df <- foreach(i = seq(nrow(sharing_df)), 
                        .packages = c("tidyverse", "castor"),
                        .combine = "bind_rows") %dopar% {
  row <- sharing_df[i, ]
  row %>% 
    mutate(host_dist = get_pairwise_distances(contree, row$host1, row$host2))
}

# Get disease citations
host_list <- unique(c(sharing_df$host1, sharing_df$host2))

pubmed_df <- foreach(sp = host_list,
                          .packages = c("tidyverse", "easyPubMed"),
                          .combine = "bind_rows") %dopar% {
  # sp <- "Sus scrofa"
  query_string <- str_glue("{sp}[all] AND (virus[all] OR pathogen[all] OR [disease])")
  pubmed_search <- get_pubmed_ids(query_string)

  tibble(host_species = sp, pubmed_count = pubmed_search$Count)
}

pubmed_parsed <- pubmed_df %>%
  mutate(pubmed_count = as.numeric(pubmed_count))

# fwrite(pubmed_parsed,
#        str_glue("results/ecology_out/pubmed_counts.csv"))
pubmed_parsed <- fread(str_glue("results/ecology_out/pubmed_counts.csv"))

merged_df <- dist_df %>%
  left_join(pubmed_parsed %>% select(host1 = host_species, host1_n = pubmed_count)) %>%
  left_join(pubmed_parsed %>% select(host2 = host_species, host2_n = pubmed_count)) %>%
  mutate(mean_count = 0.5 * (host1_n + host2_n))

pois_mod <- glm(n ~ log(mean_count + 1) + host_dist,
    data = merged_df,
    family = poisson())

summary(pois_mod)

merged_df %>%
  ggplot(aes(x = host_dist, y = log(n))) +
  geom_point() +
  geom_smooth(method = "lm")

perm_morsels <- foreach(i = seq(1000), .packages = c("tidyverse"), .combine = "c") %dopar% {
  temp <- merged_df %>%
    mutate(n = sample(n, nrow(merged_df), replace = F))
  
  pois_temp <- glm(n ~ log(mean_count + 1) + host_dist,
                  data = temp,
                  family = poisson())
  
  slope <- coef(summary(pois_temp))[3, 1]
  return(slope)
}

hist(unlist(perm_morsels))

obs_slope <- coef(summary(pois_mod))[3, 1]
obs_slope
p_val <- sum(unlist(perm_morsels) < obs_slope) / 1000

## WITHOUT HUMANS ##
nonhuman_df <- merged_df %>%
  filter(host1 != "Homo sapiens",
         host2 != "Homo sapiens")

nonhuman_df %>% nrow() / merged_df %>% nrow() * 100

pois_mod <- glm(n ~ log(mean_count + 1) + host_dist,
                data = nonhuman_df,
                family = poisson())

summary(pois_mod)

merged_df %>%
  ggplot(aes(x = host_dist, y = log(n))) +
  geom_point() +
  geom_smooth(method = "lm")

perm_morsels <- foreach(i = seq(1000), .packages = c("tidyverse"), .combine = "c") %dopar% {
  temp <- merged_df %>%
    mutate(n = sample(n, nrow(merged_df), replace = F))
  
  pois_temp <- glm(n ~ log(mean_count + 1) + host_dist,
                   data = temp,
                   family = poisson())
  
  slope <- coef(summary(pois_temp))[3, 1]
  return(slope)
}

stopCluster(cl)

hist(unlist(perm_morsels))

obs_slope <- coef(summary(pois_mod))[3, 1]
obs_slope
p_val <- sum(unlist(perm_morsels) < obs_slope) / 1000
p_val
