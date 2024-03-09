test <- fread("results/clique_classification_out/final_cluster_metadata.220723.csv") %>%
  filter(species == "Alphainfluenzavirus influenzae")
test2 <- fread("data/metadata/Orthomyxoviridae.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(species == "Alphainfluenzavirus influenzae")

host_meta <- fread("data/metadata/parsed_host_metadata.csv")

human_test2 <- test2 %>%
  left_join(host_meta) %>%
  filter(host_genus == "Homo") %>%
  distinct(country, host, isolation_source, collection_date) %>%
  sample_n(1000)
 
animal_test2 <- test2 %>%
  left_join(host_meta) %>%
  filter(host_genus != "Homo") %>%
  distinct(country, host, isolation_source, collection_date)

test2_filt <- bind_rows(animal_test2, human_test2)
test2_filt

colnames(test2)
test2 %>%
  filter(Length > 2000) %>%
  filter(Species == "Alphainfluenzavirus influenzae") %>%
  distinct(Host)

plot_df1 <- test %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(20)

plot_df2 <- test2_filt %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  head(20)

plt1 <- plot_df1 %>%
  mutate(country = factor(country, unique(plot_df1$country))) %>%
  ggplot(aes(x = country, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Original dataset")

plt2 <- plot_df2 %>%
  mutate(country = factor(country, unique(plot_df2$country))) %>%
  ggplot(aes(x = country, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Full dataset")

ggpubr::ggarrange(plt1, plt2)

plt1 <- test2_filt %>%
  ggplot(aes(x = as.Date(collection_date))) +
  geom_histogram(bins = 100) +
  xlim(as.Date("1900-01-01"), as.Date("2023-12-01")) +
  labs(x = "Collection date", title = "Full dataset")

plt2 <- test %>%
  ggplot(aes(x = as.Date(imputed_date))) +
  geom_histogram(bins = 100) +
  xlim(as.Date("1900-01-01"), 
       as.Date("2023-12-01")) +
  labs(x = "Collection date", title = "Original dataset")

require(ggpubr)
ggarrange(plt2, 
          plt1)


hist(as.Date(test2_filt$collection_date), breaks = "year")
hist(as.Date(test$imputed_date), breaks = "year", )

median(as.Date(test2$Collection_Date), na.rm = T)
median(as.Date(test$imputed_date), na.rm = T)
min(as.Date(test$imputed_date), na.rm = T)
min(as.Date(test2$Collection_Date), na.rm = T)
min(as.Date(test2$Collection_Date), na.rm = T)
max(as.Date(test2$Collection_Date), na.rm = T)
max(as.Date(test$imputed_date), na.rm = T)
test$species
hist(test2$Length)
test %>%
  filter(family == "Orthomyxoviridae") %>%
  distinct(species)

test %>%
  filter(species == "Alphainfluenzavirus influenzae") %>%
  View()

test %>%
  filter(species == "Alphainfluenzavirus influenzae") %>%
  group_by(host) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
