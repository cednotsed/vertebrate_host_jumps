setwd("c:/git_repos/vertebrate_host_jumps/")
require(easyPubMed)
require(tidyverse)
require(data.table)
require(Hmisc)
require(foreach)
require(mgcv)
require(gratia)

# Virion data
virion <- fread("data/VIRION.v0.2.1_240922/Virion.csv.gz") %>%
  filter(DetectionMethod %in% c("PCR/Sequencing", "Isolation/Observation")) %>%
  mutate(host_species = capitalize(Host),
         host_genus = capitalize(HostGenus),
         host_family = capitalize(HostFamily),
         host_order = capitalize(HostOrder),
         host_class = capitalize(HostClass)) %>%
  filter(host_class == "Mammalia")

diversity_df <- virion %>%
  group_by(host_species, host_genus, host_class, host_order, host_family) %>%
  summarise(n_species = n_distinct(Virus)) %>%
  arrange(desc(n_species)) %>%
  filter(!(host_genus %in% c("Homo", ""))) %>%
  ungroup()

# Host traits
# Population density
density_df <- fread("data/santini_2022/Dataset.140723.csv") %>% 
  rename(host_species = AcceptedName_COL) %>% 
  select(-Order.1, -Family.1, -AcceptedGenus_COL,
         -Order, -Family, -sp_id) %>% 
  relocate(host_species, .before = 1) %>% 
  group_by(host_species) %>%
  summarise(median_density = median(Density_km, na.rm = T))

density_df

# Species richness within range
richness_df <- fread("data/richness_data_tucker_2022/richness_data.140723.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>% 
  select(host_species = species, species_richness_10kmbuff, 
         human_population_density_10kmbuff, human_footprint_10kmbuff) %>%
  group_by(host_species) %>%
  summarise_if(is.numeric, median, na.rm = T)

# Population size
size_df <- fread("data/greenspoon_2023/population_data.140723.csv") %>%
  select(host_species = binomial, total_pop)

# Other traits
albery_df <- fread("data/albery_2022/host_traits.140723.csv") %>% 
  select(host_species = hHostNameFinal, hWildDomFAO, hAllZACites, 
         hMarOTerr, hDiseaseZACites, hHuntedIUCN,
         hArtfclHbttUsrIUCN, AreaHost, Area.crop,
         crop.perc, Area.pasture, urban.perc,
         Area.urban) %>%
  filter(!is.na(AreaHost)) %>%
  mutate(host_species = gsub("_", " ", host_species))

# COMBINE
traits_df <- fread("data/combine_data_soria_2021/trait_data_imputed.csv",
                   stringsAsFactors = T) %>%
  rename(host_species = iucn2020_binomial) %>%
  rename_all(~tolower(gsub(" |-", "_", .x))) %>%
  rename_all(~gsub("-_", "", .x)) %>%
  dplyr::rename(host_genus = genus) %>%
  filter(host_species %in% diversity_df$host_species)

# Filter traits
prop_na <- colSums(is.na(traits_df)) / nrow(traits_df)
to_keep <- names(prop_na)[prop_na < 0.2]
to_keep <- to_keep[grepl("_d|_g|_n|_mm|km|dphy|foraging|trophic|host_species", to_keep)]
to_keep <- to_keep[!(to_keep %in% c("family", "host_genus", "species", 
                                    "phylacine_binomial", "order"))]

to_keep
traits_filt <- traits_df %>%
  select(all_of(to_keep))

# # Get disease citations
# pubmed_morsels <- foreach(sp = unique(traits_df$host_species)) %do% {
#   # sp <- "Sus scrofa"
#   query_string <- str_glue("{sp}[all] AND (virus[all] OR pathogen[all] OR [disease])")
#   pubmed_search <- get_pubmed_ids(query_string)
#   
#   tibble(host_species = sp, pubmed_count = pubmed_search$Count)
# }
# 
# pubmed_df <- bind_rows(pubmed_morsels) %>%
#   mutate(pubmed_count = as.numeric(pubmed_count))
# 
# fwrite(pubmed_df,
#        str_glue("results/modelling_out/viral_diversity/virion_data/pubmed_counts.n{n_distinct(pubmed_df$host_species)}.csv"))

pubmed_test <- fread("results/modelling_out/viral_diversity/virion_data/pubmed_counts.n1233.csv")

# Merge data
merged_df <- diversity_df %>% 
  inner_join(pubmed_test)

merged_df %>%
  ggplot(aes(x = log10(pubmed_count + 1), y = log10(n_species + 1))) +
  geom_point()

diversity_lm <- lm(log10(n_species + 1) ~ log10(pubmed_count + 1),
                   data = merged_df)

corr_df <- tibble(host_species = merged_df$host_species, 
                  corr_diversity = diversity_lm$residuals)

final_df_traits <- diversity_df %>% 
  inner_join(corr_df) %>%
  inner_join(traits_filt) %>%
  mutate(trophic_level = factor(trophic_level),
         foraging_stratum = factor(foraging_stratum),
         trophic_level = factor(trophic_level),
         host_order = factor(host_order))
         # Area.urban = log10(Area.urban + 1),
         # AreaHost = log10(AreaHost + 1))

gmm1 <- gam(corr_diversity ~ s(adult_mass_g) + s(adult_body_length_mm) +
              s(max_longevity_d) + s(female_maturity_d) + s(age_first_reproduction_d) +
              s(gestation_length_d) + s(litter_size_n) + s(litters_per_year_n) +
              trophic_level + foraging_stratum,
              # s(median_density) + s(human_population_density_10kmbuff) + s(human_footprint_10kmbuff) +
              # s(species_richness_10kmbuff),
            data = final_df_traits,
            select = T)

final_df_density <- diversity_df %>% 
  inner_join(corr_df) %>%
  inner_join(traits_filt) %>%
  inner_join(density_df) %>%
  mutate(trophic_level = factor(trophic_level),
         foraging_stratum = factor(foraging_stratum),
         trophic_level = factor(trophic_level),
         host_order = factor(host_order))

final_df_density
gmm2 <- gam(corr_diversity ~ s(adult_mass_g) + s(adult_body_length_mm) +
      s(max_longevity_d) + s(female_maturity_d) + s(age_first_reproduction_d) +
      s(gestation_length_d) + s(litter_size_n) + s(litters_per_year_n) +
      trophic_level + foraging_stratum + 
      s(median_density),
      data = final_df_density,
      select = T)

summary(gmm2)
draw(gmm2, residuals = T)

final_df_human <- diversity_df %>% 
  inner_join(corr_df) %>%
  inner_join(traits_filt) %>%
  inner_join(density_df) %>%
  inner_join(richness_df) %>%
  # inner_join(size_df) %>%
  mutate(trophic_level = factor(trophic_level),
         foraging_stratum = factor(foraging_stratum),
         trophic_level = factor(trophic_level),
         host_order = factor(host_order))

gmm3 <- gam(corr_diversity ~ s(adult_mass_g) + s(adult_body_length_mm) +
              s(max_longevity_d) + s(female_maturity_d) + s(age_first_reproduction_d) +
              s(gestation_length_d) + s(litter_size_n) + s(litters_per_year_n) +
              trophic_level + foraging_stratum + 
              s(median_density) + s(species_richness_10kmbuff) + s(human_population_density_10kmbuff) +
              s(human_footprint_10kmbuff),
            data = final_df_human,
            select = T)

summary(gmm3)
draw(gmm3, residuals = T)

final_df_all <- diversity_df %>% 
  inner_join(corr_df) %>%
  inner_join(traits_filt) %>%
  inner_join(density_df) %>%
  inner_join(albery_df) %>%
  inner_join(richness_df) %>%
  mutate(trophic_level = factor(trophic_level),
         foraging_stratum = factor(foraging_stratum),
         trophic_level = factor(trophic_level),
         host_order = factor(host_order))

gmm4 <- gam(corr_diversity ~ s(adult_mass_g) + s(adult_body_length_mm) +
              s(max_longevity_d) + s(female_maturity_d) + s(age_first_reproduction_d) +
              s(gestation_length_d) + s(litter_size_n) + s(litters_per_year_n) +
              trophic_level + foraging_stratum + 
              s(median_density) + s(species_richness_10kmbuff) + s(human_population_density_10kmbuff) +
              s(human_footprint_10kmbuff) + s(Area.urban) + s(AreaHost),
            data = final_df_all,
            select = T)
test <- final_df_all %>%
  select_if(is.numeric)
reshape2::melt(cor(test)) %>%
  ggplot(aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
summary(gmm4)
draw(gmm4, residuals = T)

cor.test(final_df_all$human_footprint_10kmbuff, final_df_all$human_population_density_10kmbuff,
         method = "spearman")
