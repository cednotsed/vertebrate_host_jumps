segmented %>%
  filter(family == "Picobirnaviridae") %>% 
  filter(grepl("RdRP|polymerase|segment 2", 
               genbank_title,
               ignore.case = T) &
           grepl("complete", 
                 genbank_title, 
                 ignore.case = T)) %>% View()

segmented %>%
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
                 ignore.case = T)) %>% 
  filter(grepl("Influenza A", genbank_title, ignore.case = T))

segmented %>%
  filter(family == "Arenaviridae") %>% 
  filter(grepl("complete", 
               genbank_title, 
               ignore.case = T) &
          !grepl("partial", 
                 genbank_title,
                 ignore.case = T) &
          grepl("segment L|L gene", 
                genbank_title,
                ignore.case = T)) %>% View()
  
segmented %>%
  filter(family == "Birnaviridae") %>% 
  filter(grepl("complete", 
               genbank_title, 
               ignore.case = T) &
           !grepl("partial", 
                  genbank_title,
                  ignore.case = T) &
           grepl("polymerase|VP1|segment B|ORF1", 
                 genbank_title,
                 ignore.case = T)) %>% View()

segmented %>%
  filter(family == "Peribunyaviridae") %>% 
  filter(grepl("complete", 
               genbank_title, 
               ignore.case = T) &
           !grepl("partial", 
                  genbank_title,
                  ignore.case = T) &
           grepl("polymerase|segment L", 
                  genbank_title,
                  ignore.case = T)) %>% View()


ns_animal_viruses %>%
  filter(family == "Circoviridae") %>%
  # filter(host == "") %>% View()
  group_by(host) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% View()
#   
ns_animal_viruses <- non_segmented %>% 
  filter(host != "Homo sapiens")

final_merged %>%
  group_by(family) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>% View()


family_filt %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|SARS-CoV-2", 
               genbank_title,
               ignore.case = T) &
           host == "Homo sapiens") %>%
  filter(length > 28000) %>%
  filter(!grepl("partial|proviral|protein|structural|polymerase|segment|capsid|replicase|NSP|gene,|mrna|cds|ORF",
                genbank_title,
                ignore.case = T)) %>%
  distinct(country, host, isolation_source, collection_date, .keep_all = T)
