rm(list = ls())
setwd("C:/git_repos/vertebrate_host_jumps/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(MSA2dist)
require(doParallel)

good_alns <- fread("results/qc_out/good_alignments.csv")

jump_df <- fread("results/mutational_load_out/host_jump_lists/diff_hosts.genus_counts.all_jumps.csv") %>%
  as_tibble() %>%
  filter(clique_name %in% good_alns$clique_name)

# Reconstructed sequences
node_dir <- "results/source_sink_analysis/ancestral_reconstructions/ancestral_sequences"
node_list <- list.files(node_dir, 
                       "\\.fna",
                       full.names = T)

# Alignments
aln_dir <- "data/alignments/source_sink_mini_trees/without_outgroup.masked"
aln_list <- list.files(aln_dir, full.names = T)
aln_list <- aln_list[!grepl("\\.log|trim_to", aln_list)]

# Annotations
gff_path <- "results/dnds_out/ncbi_annotation/all_gene_annotations.290723.n53xxx.parsed.gff3"
gff_df <- fread(gff_path)

# Ensure gffs all have product info
gff_filt <- gff_df %>%
  filter(V3 == "CDS") %>%
  mutate(has_annot = grepl("product", V9)) %>%
  filter(has_annot) %>%
  as_tibble()

# Fix annotation function
fix_annot <- function(gene_aln, gene_pos, full_aln) {
  remainder <- ncol(gene_aln) %% 3
  to_add <- 3 - remainder
  
  # Add nucleotides to complete codon
  end_pos <- max(gene_pos)
  positions_fixed <- c(gene_pos, seq(end_pos + 1, end_pos + to_add))
  
  # If adding nucleotides exceeds aln length, subtract instead
  if (max(positions_fixed) > ncol(full_aln)) {
    positions_fixed <- gene_pos[1:(length(gene_pos) - remainder)]
  }
  
  gene_temp <- full_aln[, positions_fixed]
  
  # Return only if annot starts with Met
  seq_temp <- apply(gene_temp, 1, paste0, collapse = "")
  seq_temp <- gsub("-", "N", seq_temp)
  prot_temp <- Biostrings::translate(DNAStringSet(seq_temp),
                                     if.fuzzy.codon = "solve")
  
  good_start <- all(as.matrix(prot_temp)[, 1] == "M")
  
  if (good_start) {
    return(gene_temp)
  } else {
    return(NA)
  }
}

cl <- makeCluster(12)
registerDoParallel(cl)

clique_list <- jump_df %>% 
  distinct(clique_name, is_jump)

morsels <- foreach(i = seq(nrow(clique_list)), .packages = c("data.table",
                                                         "tidyverse",
                                                         "MSA2dist",
                                                         "Biostrings",
                                                         "foreach")) %dopar% {
  # i = 7
  # clique_name <- clique_list[i, ]
  # tip_name <- "MN737436.1"
  # node_name <- "Node205"
  
  row <- clique_list[i, ]
  
  clique <- row$clique_name
  
  print(str_glue("Creating mega CDS for {clique}"))
  jump_bool <- row$is_jump
  
  jump_filt <- jump_df %>%
    filter(clique_name == clique,
           is_jump == jump_bool)
  
  # Count number of accessions in mega cds
  acc_counter <- 0
  
  mega_cds <- foreach(q = seq(nrow(jump_filt)), .combine = "cbind") %do% {
    # q = 1
    subrow <- jump_filt[q, ]
    tip_name <- subrow$tip_name
    node_name <- subrow$anc_name
    tip_state <- subrow$tip_state
    anc_state <- subrow$anc_state
    
    # print(str_glue("{clique}: {node_name} -> {tip_name}"))
  
    aln_path <- aln_list[grepl(str_glue("{clique}\\."), aln_list)]
  
    tip_aln <- readDNAStringSet(aln_path)
    # Parse accession names _R_
    names(tip_aln) <- gsub("_R_", "", names(tip_aln))
    
    node_path <- node_list[grepl(str_glue("{clique}."), node_list)]
    node_aln <- readDNAStringSet(node_path)
    
    tip_seq <- tip_aln[tip_name]
    node_seq <- node_aln[node_name]
    
    # Trim gaps from tip sequence
    aln_pair <- as.matrix(c(tip_seq, node_seq))
    to_keep <- aln_pair[1, ] != "-"
    
    aln_filt <- aln_pair[, to_keep]
    
    # Get CDS
    genes <- gff_filt %>%
      filter(V1 == tip_name)
    
    # Check if any annotation is available at all
    if(nrow(genes) > 0) {
      gene_counter <- 0
      fix_counter <- 0
      bad_counter <- 0
      
      # Iterate through genes
      annot_list <- deframe(genes %>% distinct(V9))
      
      sub_cds <- foreach(annot_string = annot_list, .combine = "cbind") %do% {
        # annot_string = annot_list[2]
        # Get gene annotation
        annot <- str_split(annot_string, ";")[[1]]
        gene_annot <- annot[grepl("product", annot)]
        gene_annot <- gsub("product=", "", gene_annot)
        
        subgene_list <- genes %>%
          filter(V9 == annot_string)
        
        # Get gene positions
        gene_positions <- foreach(k = seq(nrow(subgene_list)), .combine = "c") %do% {
          positions <- subgene_list[k, ]$V4:subgene_list[k, ]$V5
          return(positions)
        }
        
        # print(str_glue("For {tip_name}, {gene_annot} is {length(gene_positions)}nt!!"))
  
        # Check if gene positions exceed aln length
        if(max(gene_positions) > ncol(aln_filt)) {
          print(str_glue("ERROR!!for {tip_name}, {gene_annot} annot exceeds aln length"))
          bad_counter <- bad_counter + 1
          return(NULL)
        }
        
        # Check if gene length is multiple of 3
        bad_length <- length(gene_positions) %% 3 != 0
        
        gene_aln <- aln_filt[, gene_positions]
        
        if(bad_length) {
          gene_fixed <- fix_annot(gene_aln, gene_positions, aln_filt)
          
          # Add counter if gene was fixed
          if (!any(is.na(gene_fixed))) {
            fix_counter <- fix_counter + 1
            
            nt_added <- ncol(gene_fixed) - ncol(gene_aln)
            print(str_glue("for {tip_name}, {gene_annot}: added {nt_added}nt"))
          } else {
            # Otherwise return nothing
            bad_counter <- bad_counter + 1
            print(str_glue("for {tip_name}, could not fix {gene_annot}"))
            return(NULL)
          }
        } else {
          # If length is good then return original gene aln
          gene_fixed <- gene_aln
        }
        
        # Reverse complement if gene on minus strand
        strand_val <- unique(subgene_list$V7)
        
        if (length(strand_val) == 1) {
          if (strand_val == "-") {
            gene_rc <- DNAStringSet(apply(gene_fixed, 1, paste0, collapse = ""))
            gene_rc <- reverseComplement(gene_rc)
            gene_fixed <- as.matrix(gene_rc)
          } 
          
          gene_counter <- gene_counter + 1
          # All good, return the gene!
          # But first remove all ambiguous codons
          gene_fixed <- foreach(codon_start = seq(1, ncol(gene_fixed), by = 3),
                                .combine = "cbind") %do% {
            codon <- gene_fixed[, codon_start:(codon_start + 2)]
            if (any(codon %in% c("-", "N"))) {
              return(NULL)
            } else {
              return(codon)
            }
          }
          
          # Then check if gene still has positive length
          if (!is.null(gene_fixed)) {
            # gene_annot_list <- c(gene_annot_list, gene_annot)
            # print(str_glue("{ncol(gene_fixed) %% 3 == 0}"))
            return(gene_fixed)
          } else {
            return(NULL)
            print(str_glue("ERROR!for {tip_name}, {gene_annot}: all ambiguous!!"))
          }        
        } else {
          # Return nothing if strand annotation is bad
          print(str_glue("ERROR!for {tip_name}, {gene_annot}: subgenes on different strands??"))
          return(NULL)
        }
      }
      
      if(ncol(sub_cds) > 0) {
        acc_counter <- acc_counter + 1
        return(sub_cds)
      } else {
        return(NULL)
      }
    }
  }
  
  # ncol(mega_cds)
  # Concatenate all gene alignments and calculate dnds
  if(!is.null(mega_cds)) {
    cds_length <- ncol(mega_cds)
    print(str_glue("{clique} has a mega CDS of {cds_length}nt comprising {acc_counter} accessions"))
    concat_seq <- DNAStringSet(apply(mega_cds, 1, paste0, collapse = ""))
    dnds <- dnastring2kaks(concat_seq, model = "Li")
    ka <- as.numeric(dnds[, "ka"])
    ks <- as.numeric(dnds[, "ks"])
    
    return(tibble(clique_name = clique, 
                  n_accs = acc_counter,
                  is_jump = jump_bool,
                  cds_length = cds_length,
                  ka = ka,
                  ks = ks))
  } else {
    return(NULL)
  }
}

stopCluster(cl)

res <- bind_rows(morsels)

# res %>%
  # arrange(desc(n_bad))
# 
fwrite(res, "results/dnds_out/dnds_out.diff_hosts.genus_counts.all_jumps.mega.csv")
