setwd("c:/git_repos/viral_sharing/")
require(data.table)
require(tidyverse)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/all_viruses.140423.filt.csv")

spec_to_remove <- deframe(meta %>% 
  distinct(species) %>%
  filter(grepl("Gokushovirinae environmental samples|Temperate fruit|Gossypium|Emaravirus|Stralarivirus|Rice ragged|Nerrivikvirus|Peanut|Cucurbit|Cacao|Pea necrotic|white spot|Bhendi yellow|Prokaryotic|Macroptilium|Camellia|Chickpea|Beet|Chrysanthemum|avocado|Tobacco|southern rice|strawberry|prunus|Shallot|Rice stripe tenuivirus|Broad bean|sugarcane|Beet curly top virus|onion|Prune dwarf virus|garlic|Blainvillea|Faba bean|banana|Digitaria streak|Hop stunt|wheat|Leek yellow stripe virus|Rice black streaked|leaf|grapevine|pepper|plum|melon|Cardamom|polar freshwater|peach|cherry|tomato|potato|citrus|cucumber|apple|mosaic|maize", species, ignore.case = T)))

meta_filt <- meta %>% 
  filter(!grepl("Cystoviridae|Partitiviridae|Marnaviridae|Cystoviridae|Becurtovirus|Sobemovirus|Begomovirus|Mitoviridae|Caudovirales|Gokushovirinae|Totiviridae|Inovirus|Botourmiaviridae|Narnaviridae|Leviviridae|Myoviridae|siphoviridae|Geminiviridae|Genomoviridae|Nanoviridae|Microviridae|Microvirus|Rice|phage|Inoviridae", species, ignore.case = T)) %>%
  filter(!grepl("Zanthoxylum|Haloquadratum|Areca catechu|Petroselinum|Ageratum|Brassica|Sophora|Medicago|Sulfolobus|Vigna|Sida acuta|Solanum", host)) %>% 
  filter(!(species %in% spec_to_remove)) %>%
  filter(host != "")

aai_df <- fread("results/checkv_out/all_viruses.140423.filt/completeness.tsv") %>%
  dplyr::rename(accession = contig_id) %>%
  right_join(meta_filt)

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

genome_filt <- aai_df %>%
  filter(family %in% to_keep) %>%
  filter(aai_completeness > 95) %>%
  select(accession, genbank_title, family, 
         genus, species, host, 
         isolation_source, collection_date, imputed_date, 
         country)

large <- deframe(genome_filt %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>% 
  filter(n > 500) %>%
  distinct(family))

others <- genome_filt %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>% 
  filter(n <= 500)

fna <- readDNAStringSet("data/genomes/all_viruses.140423.filt.fna")
names(fna) <- str_split(names(fna), " ", simplify = T)[, 1]

fna_filt <- fna[names(fna) %in% genome_filt$accession]
nrow(genome_filt) == length(fna_filt)

fwrite(genome_filt, "data/metadata/all_viruses.140423.filt.QCed.csv")
writeXStringSet(fna_filt, "data/genomes/all_viruses.140423.filt.QCed.fna")

## AFTER FORMATTING ##
fna_filt2 <- readDNAStringSet("data/genomes/all_viruses.140423.filt.QCed.formatted.fna")

foreach(fam = large) %do% {
  temp_filt <- genome_filt %>%
    filter(family == fam)
  
  temp_fna <- fna_filt2[names(fna_filt2) %in% temp_filt$accession]  
  writeXStringSet(temp_fna, str_glue("data/genomes/viral_family_subsets/{fam}.140423.filt.QCed.formatted.fna"))
}

# Others
others_filt <- genome_filt %>%
  filter(family %in% others$family)

others_fna <- fna_filt2[names(fna_filt2) %in% others_filt$accession]  
writeXStringSet(others_fna, str_glue("data/genomes/viral_family_subsets/others.140423.filt.QCed.formatted.fna"))
