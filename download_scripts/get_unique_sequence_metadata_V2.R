setwd("c:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)

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

family_filt <- dat %>%
  filter(family %in% to_keep) %>%
  filter(!grepl("partial|proviral|protein|structural|polymerase|segment|capsid|replicase|NSP|gene,|mrna|cds|ORF", 
                genbank_title,
                ignore.case = T))

# Remove sars-cov-2 samples
sars2_ref <- family_filt %>%
  filter(accession == "MN908947.3")

human_sars2 <- family_filt %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
               genbank_title,
               ignore.case = T) &
           host == "Homo sapiens") %>%
  filter(length > 28000) %>%
  distinct(country, host, isolation_source, .keep_all = T)

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

# Filter other viruses (no sars2)
human_viruses <- no_sars2 %>%
  filter(host == "Homo sapiens") %>%
  distinct(species, country, isolation_source, collection_date, .keep_all = T)

animal_viruses <- no_sars2 %>% 
  filter(host != "Homo sapiens")

print(str_glue("after family filtering: {nrow(family_filt)}"))
print(str_glue("Animal viruses (sars2): {nrow(animal_sars2)}"))
print(str_glue("Human viruses (sars2): {nrow(human_sars2)}"))
print(str_glue("Animal viruses (no sars2): {nrow(animal_viruses)}"))
print(str_glue("Human viruses (no sars2): {nrow(human_viruses)}"))

final_merged <- bind_rows(sars2_ref, animal_sars2, human_sars2,
                          animal_viruses, human_viruses)
final_merged %>%
  fwrite("data/metadata/all_viruses.220723.filt.csv")

final_merged %>%
  select(accession) %>%
  fwrite("data/metadata/all_viruses.220723.filt.accessions_only.txt",
         col.names = F,
         eol = "\n")
