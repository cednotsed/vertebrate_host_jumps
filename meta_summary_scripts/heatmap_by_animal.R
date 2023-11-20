rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(ggrepel)
require(Hmisc)
require(maps)

meta <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.220723.V2.parsed_host.csv") %>%
  select(accession, host, host_genus, 
         host_family, genbank_title, country, 
         family, species) %>%
  filter(host_genus != "Homo") %>%
  filter(host_family != "Hominidae") %>%
  filter(host_genus != "")

count_df <- meta %>%
  # filter(species == "Alphainfluenzavirus influenzae") %>%
  group_by(host_genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(4)

for (host_name in unique(count_df$host_genus)) {
  # Proportion most sequenced viruses
  virus_count <- meta %>%
    filter(host_genus == host_name) %>% 
    filter(country != "") %>%
    group_by(species) %>%
    summarise(n_genomes = n()) %>%
    arrange(desc(n_genomes)) %>%
    mutate(perc_genomes = signif(n_genomes / sum(n_genomes) * 100, 2))
  
  virus_count %>%
    fwrite(str_glue("results/meta_summary_out/surveillance_effort.{host_name}.csv"))
  
  plot_df <- meta %>%
    filter(host_genus == host_name) %>% 
    filter(country != "") %>%
    group_by(country) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    mutate(country = ifelse(country == "United Kingdom", "UK", country)) %>%
    mutate(country = ifelse(grepl("Viet", country), "Vietnam", country)) %>%
    rename(region = country)
  
  # Barplot
  plot_df %>%
    arrange(desc(n)) %>%
    head(10) %>%
    mutate(region = factor(region, rev(unique(plot_df$region)))) %>%
    ggplot(aes(x = n, y = region)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18, angle = 45, hjust = 1)) +
    labs(x = "Viral sequences",
         title = host_name)
  
  ggsave(str_glue("results/meta_summary_out/surveillance_effort.{host_name}.pdf"), 
         dpi = 600,
         width = 4, height = 5)
  
  # Heatmap
  world_map <- map_data("world")
  
  plot_map <- subset(world_map, region != "Antarctica") %>%
    left_join(plot_df) %>%
    mutate(n = ifelse(is.na(n)|n == 0, 1, n))
  
  ggplot() +
    geom_map(data = plot_map, 
             map = plot_map,
             aes(x = long, y = lat, map_id = region,
                 fill = log10(n)),
             size = 0.25,
             color = "black") +
    scale_fill_gradient2(low = "blue", 
                         mid = "#FFF9c4",
                         high = "red",
                         midpoint = 3,
                         limits = c(0, 6)) +
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
    labs(fill = "Log10(viral sequences)",
         title = host_name)

  ggsave(str_glue("results/meta_summary_out/world_map_surveillance_effort.{host_name}.pdf"), 
         dpi = 600,
         width = 9, height = 5) 
}
