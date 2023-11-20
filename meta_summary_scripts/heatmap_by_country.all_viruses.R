rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)
require(Hmisc)
require(maps)
require(foreach)

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.parsed_host.csv") %>%
  select(accession, host, host_genus, genbank_title, country, family, species)

meta %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(4)

# SC2 by country
sc2_df <- meta %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
               genbank_title,
               ignore.case = T)) %>%
  filter(country != "") %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
  mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
  rename(region = country)

# Influenza A
flu_df <- meta %>%
  filter(species == "Alphainfluenzavirus influenzae") %>%
  filter(country != "") %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
  mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
  rename(region = country)

hiv_df <- meta %>%
  filter(species == "Human immunodeficiency virus 1") %>%
  filter(country != "") %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
  mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
  rename(region = country)

hep_df <- meta %>%
  filter(species == "Hepacivirus hominis") %>%
  filter(country != "") %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
  mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
  rename(region = country)

world_map <- map_data("world")

sc2_map <- subset(world_map, region != "Antarctica") %>%
  left_join(sc2_df) %>%
  mutate(n = ifelse(is.na(n)|n == 0, 1, n))

flu_map <- subset(world_map, region != "Antarctica") %>%
  left_join(flu_df) %>%
  mutate(n = ifelse(is.na(n)|n == 0, 1, n))

hiv_map <- subset(world_map, region != "Antarctica") %>%
  left_join(hiv_df) %>%
  mutate(n = ifelse(is.na(n)|n == 0, 1, n))

hep_map <- subset(world_map, region != "Antarctica") %>%
  left_join(hep_df) %>%
  mutate(n = ifelse(is.na(n)|n == 0, 1, n))

plot_list <- list(SC2 = sc2_map, flu = flu_map, hiv = hiv_map, HCV = hep_map)

for (i in seq(length(plot_list))) {
  plot_name <- names(plot_list)[i]
  dat <- plot_list[[i]]
  ggplot() +
    geom_map(data = dat, 
             map = dat,
             aes(x = long, y = lat, map_id = region,
                 fill = log10(n)),
             size = 0.25,
             color = "black") +
    scale_fill_gradient2(low = "blue", 
                         mid = "#FFF9c4",
                         high = "red",
                         midpoint = 3.5,
                         limits = c(0, 7)) +
    
    theme_classic() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank()) +
    labs(fill = "Log10(viral sequences)")
  
  ggsave(str_glue("results/meta_summary_out/world_map_surveillance_effort.{plot_name}.pdf"), 
         dpi = 600,
         width = 9, height = 5) 
}

# Barplot
# plot_df <- meta %>%
#   filter(host_genus == host_name) %>% 
#   filter(country != "") %>%
#   group_by(country) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n)) %>%
#   mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
#   mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
#   rename(region = country)

# Barplot
bar_list <- list(SC2 = sc2_df, flu = flu_df, hiv = hiv_df, HCV = hep_df)

foreach(i = seq(length(bar_list))) %do% {
  dat <- bar_list[[i]]
  virus_name <- names(bar_list)[i]
  dat %>%
    arrange(desc(n)) %>%
    head(10) %>%
    mutate(region = factor(region, rev(unique(dat$region)))) %>%
    ggplot(aes(x = n, y = region)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18, angle = 45, hjust = 1)) +
    labs(x = "Viral sequences",
         title = virus_name)
  
  ggsave(str_glue("results/meta_summary_out/surveillance_effort_barplot.{virus_name}.pdf"), 
         dpi = 600,
         width = 4, height = 5)
}
