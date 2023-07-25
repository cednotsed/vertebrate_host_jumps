rm(list = ls())
setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)

## PLEASE RMBR TO SET SEED ##
set.seed(66)

dat <- fread("data/metadata/all_viruses.taxid10239.excl_provirus_env_lab_vax.gt1000nt.220723.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))

# Remove non-vertebrate viruses
to_keep <- c("Poxviridae", "Herpesviridae", "Picornaviridae",
             "Peribunyaviridae", "Birnaviridae", "Circoviridae",
             "Rhabdoviridae", "Hantaviridae", "Polyomaviridae",
             "Paramyxoviridae", "Bornaviridae", "Filoviridae",
             "Nairoviridae", "Spinareoviridae", "Sedoreoviridae",
             "Phenuiviridae", "Arenaviridae", "Parvoviridae",
             "Redondoviridae", "Adenoviridae", "Togaviridae",
             "Orthomyxoviridae", "Arteriviridae", "Coronaviridae",
             "Tobaniviridae", "Nanhypoviridae", "Olifoviridae",
             "Gresnaviridae", "Nanghoshaviridae", "Cremegaviridae",
             "Caliciviridae", "Asfarviridae", "Alloherpesviridae",
             "Astroviridae", "Flaviviridae", "Retroviridae",
             "Nyamiviridae", "Metaviridae", "Hepadnaviridae",
             "Genomoviridae", "Smacoviridae", "Picobirnaviridae",
             "Reoviridae", "Anelloviridae", "Papillomaviridae",
             "Adomaviridae", "Hepeviridae", "Kolmioviridae",
             "Pneumoviridae", "Amnoonviridae", "Nodaviridae",
             "Matonaviridae", "Pseudoviridae")

segmented_families <- c("Arenaviridae", "Birnaviridae", "Peribunyaviridae", 
                        "Orthomyxoviridae", "Picobirnaviridae", "Reoviridae")

family_filt <- dat %>%
  filter(family %in% to_keep) %>%
  filter(family != "")

# Remove families with less than 100 genomes
family_count <- family_filt %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  filter(n > 100)

family_filt <- family_filt %>%
  filter(family %in% family_count$family)
  
# Remove sars-cov-2 samples
sars2_ref <- family_filt %>%
  filter(accession == "MN908947.3")

human_sars2 <- family_filt %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
               genbank_title,
               ignore.case = T) &
           host == "Homo sapiens") %>%
  filter(!grepl("partial|proviral|protein|structural|polymerase|segment|capsid|replicase|NSP|gene,|mrna|cds|ORF",
                genbank_title,
                ignore.case = T)) %>%
  filter(length > 28000) %>%
  distinct(country, host, isolation_source, collection_date, .keep_all = T) %>%
  sample_n(1000, replace = F)

animal_sars2 <- family_filt %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
               genbank_title,
               ignore.case = T) &
           host != "Homo sapiens") %>%
  filter(host != "") %>%
  filter(length > 28000)

no_sars2 <- family_filt %>%
  filter(!(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
                 genbank_title,
                 ignore.case = T)))

# Get segmented and non-segmented viruses
segmented <- no_sars2 %>%
  filter(family %in% segmented_families)

non_segmented <- no_sars2 %>%
  filter(!(family %in% segmented_families)) %>%
  filter(!grepl("partial|proviral|protein|structural|polymerase|segment|capsid|replicase|NSP|gene,|mrna|cds|ORF",
                genbank_title,
                ignore.case = T))

# Filter segmented
pico <- segmented %>%
  filter(family == "Picobirnaviridae") %>% 
  filter(grepl("RdRP|polymerase|segment 2", 
               genbank_title,
               ignore.case = T) &
           grepl("complete", 
                 genbank_title, 
                 ignore.case = T))

arena <- segmented %>%
  filter(family == "Arenaviridae") %>% 
  filter(grepl("complete", 
               genbank_title, 
               ignore.case = T) &
           !grepl("partial", 
                  genbank_title,
                  ignore.case = T) &
           grepl("segment L|L gene", 
                 genbank_title,
                 ignore.case = T))

birna <- segmented %>%
  filter(family == "Birnaviridae") %>% 
  filter(grepl("complete", 
               genbank_title, 
               ignore.case = T) &
           !grepl("partial", 
                  genbank_title,
                  ignore.case = T) &
           grepl("polymerase|VP1|segment B|ORF1", 
                 genbank_title,
                 ignore.case = T))

peri <- segmented %>%
  filter(family == "Peribunyaviridae") %>% 
  filter(grepl("complete", 
               genbank_title, 
               ignore.case = T) &
           !grepl("partial", 
                  genbank_title,
                  ignore.case = T) &
           grepl("polymerase|segment L", 
                 genbank_title,
                 ignore.case = T))

# Filter Flu B sequences
ortho <- segmented %>%
  filter(family == "Orthomyxoviridae") %>% 
  filter(!grepl("partial",
                genbank_title,
                ignore.case = T)) %>%
  filter(grepl("PB 1|PB1 |polymerase basic 1", 
               genbank_title,
               ignore.case = T) &
           !grepl("PB1-F2", 
                  genbank_title,
                  ignore.case = T) &
           grepl("complete", 
                 genbank_title, 
                 ignore.case = T))

ortho_human <- ortho %>%
  filter(species == "Betainfluenzavirus influenzae") %>%
  filter(host == "Homo sapiens") %>%
  distinct(country, host, isolation_source, collection_date, .keep_all = T)

ortho_animal <- ortho %>%
  filter(species == "Betainfluenzavirus influenzae") %>%
  filter(host != "Homo sapiens")

ortho_filt <- ortho %>%
  filter(species != "Betainfluenzavirus influenzae") %>%
  bind_rows(ortho_human, ortho_animal)
  
segmented_filt <- bind_rows(pico, birna, peri, ortho_filt, arena)

# Filter other viruses non-segmented (no sars2)
ns_human_viruses <- non_segmented %>%
  filter(host == "Homo sapiens") %>%
  distinct(species, country, isolation_source, collection_date, .keep_all = T)

ns_animal_viruses <- non_segmented %>% 
  filter(host != "Homo sapiens")

# Merge all
final_merged <- bind_rows(sars2_ref, animal_sars2, human_sars2,
                          segmented_filt, ns_animal_viruses, ns_human_viruses) %>%
  mutate(is_segmented = ifelse(family %in% segmented_families, T, F))

print(str_glue("after family filtering: {nrow(family_filt)}"))
print(str_glue("Animal viruses (sars2): {nrow(animal_sars2)}"))
print(str_glue("Human viruses (sars2): {nrow(human_sars2)}"))
print(str_glue("Segmented viruses: {nrow(segmented_filt)}"))
print(str_glue("Animal viruses (non-segmented): {nrow(ns_animal_viruses)}"))
print(str_glue("Human viruses (non-segmented): {nrow(ns_human_viruses)}"))
print(str_glue("Final filtered viruses: {nrow(final_merged)}"))

# Check if script made duplicates
all(final_merged$accession == unique(final_merged$accession))

# Annotate circular genomes
segmented_family <- c("Arenaviridae", "Birnaviridae", "Peribunyaviridae",
                        "Orthomyxoviridae", "Picobirnaviridae", "Reoviridae")
circular_family <- c("Anelloviridae", "Circoviridae", "Genomoviridae",
                     "Hepadnaviridae", "Smacoviridae")
final_merged %>%
  mutate(is_segmented = ifelse(family %in% segmented_family, T, F),
         is_circular = ifelse(family %in% circular_family, T, F)) %>%
  fwrite("data/metadata/all_viruses.220723.filt.csv")

final_merged %>%
  select(accession) %>%
  fwrite("data/metadata/all_viruses.220723.filt.accessions_only.txt",
         col.names = F,
         eol = "\n")

