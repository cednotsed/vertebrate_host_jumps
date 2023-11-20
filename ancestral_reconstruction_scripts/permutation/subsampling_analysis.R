rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  distinct(host, host_genus)
res_dir <- str_glue("results/ancestral_reconstruction_out/subsampling_analysis/temp_results")
file_list <- list.files(res_dir, full.names = T)

df <- foreach(file = file_list, .combine = "bind_rows") %do% {
  fread(file)
}

df %>%
  left_join(meta %>% select(anc_state = host, anc_genus = host_genus)) %>%
  left_join(meta %>% select(tip_state = host, tip_genus = host_genus)) %>%
  filter(anc_genus == "Homo" | tip_genus == "Homo") %>%
  mutate(event_type = ifelse(anc_genus != "Homo", "Zoonotic", "Anthroponotic")) %>%
  group_by(n_human, iter) %>%
  summarise(n_anthro = sum(event_type == "Anthroponotic"),
            n_zoo = sum(event_type == "Zoonotic")) %>%
  mutate(ratio = n_anthro / n_zoo,
         prop_human = round(n_human / (579 + n_human), 3)) %>%
  ggplot(aes(x = factor(prop_human), y = ratio)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, 
             lty = "dashed",
             color = "red") +
  theme_bw() +
  labs(x = "Prop. human genomes", y = "Anthroponotic-zoonotic ratio")

  ggsave("results/source_sink_analysis/SC2_ratio_by_threshold.pdf",
         dpi = 600,
         width = 5, height = 5)
100/570


