# Crossing host boundaries: the evolutionary drivers and correlates of viral host jumps
Authors: Cedric C.S. Tan (cedriccstan@gmail.com), Lucy van Dorp, Francois Balloux

DOI: https://doi.org/10.1101/2023.09.01.555953

## Clique assignment pipeline
1. Download sequence metadata from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) (TaxID=10239, excluding provirus, environmental, lab host and vaccine strain sequences).
1. [Sequence metadata curation](download_scripts/get_unique_sequence_metadata_V2.R)
2. [Download sequences](download_scripts/download_genomes_using_acc_list.sh)
3. [Use CheckV to QC sequences and remove low quality sequences](qc_scripts)
4. [Extract sequences for each viral family](phylogenetic_scripts/extract_viral_family_genomes.R)
5. [Mash genomes within each family](phylogenetic_scripts/run_mash_multiple.sh)
6. [Create NJ trees and run InfoMap to generate final clusters](clique_classification_scripts/generate_final_clusters.R)
7. Assess performance metrics:
   * [Concordance with NCBI taxonomy](clique_classification_scripts/optimise_threshold_for_metrics.R)
   * [Concordance with ICTV taxonomy](clique_classification_scripts/optimise_threshold_for_metrics.ICTV.R)
   * [Clique monophyly](clique_classification_scripts/optimise_threshold_for_monophyly_V2.R)

## Tree rooting
1. [Extract sequences for viral cliques](ancestral_reconstruction_scripts/rooting/get_mini_tree_genomes.all_jumps.even_further.R), including 10 sequences where the minimum Mash distance to any sequence in the clique of interest is 0.3-0.5.
2. [Identify suitable tips to root viral cliques](ancestral_reconstruction_scripts/rooting/get_roots_for_clique_trees.even_further.R)

## Phylogenetic reconstruction
1. [Align](phylogenetic_scripts/align_multiple_genomes_mafft.sh) clique sequences
2. [Mask](phylogenetic_scripts/mask_alignment.R) gappy positions
3. [Build tree](phylogenetic_scripts/build_trees_with_bootstraps_myriad.sh) with ancestral sequence reconstruction toggled

## Ancestral host reconstruction and host jump inference
1. [Identifying host jumps](ancestral_reconstruction_scripts/varying_thresholds/ancestral_reconstruction.varying_threshold.R)
2. [Identifying non-host jumps](ancestral_reconstruction_scripts/varying_thresholds/ancestral_reconstruction.varying_threshold.R)
3. [Parse and filter](ancestral_reconstruction_scripts/varying_thresholds/get_jump_non_jump.V2.R) host jumps and non-host jumps

## Main analyses
### Metadata summary
[here](meta_summary_scripts) and [here](host_range_scripts)
### Anthroponotic frequency
[here](source_sink_scripts)
### Mutational distance
[here](mutational_load_scripts)
### dN/dS
[here](dnds_scripts)
