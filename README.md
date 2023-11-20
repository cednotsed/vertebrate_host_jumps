# Crossing host boundaries: the evolutionary drivers and correlates of viral host jumps
Authors: Cedric C.S. Tan, Lucy van Dorp, Francois Balloux

DOI: https://doi.org/10.1101/2023.09.01.555953

## Rough data pipeline
1. Download sequence metadata from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) (TaxID=10239, excluding provirus, environmental, lab host and vaccine strain sequences).
1. [Sequence metadata curation](download_scripts/get_unique_sequence_metadata_V2.R)
2. [Download sequences](download_scripts/download_genomes_using_acc_list.sh)
3. [Use CheckV to QC sequences and remove low quality sequences](qc_scripts)
4. 

