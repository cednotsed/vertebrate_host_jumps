#acc_list=../data/metadata/all_viruses.genbank_complete_excllabhost_exclvaccine.191222.unique.accession_only.txt
#out_dir=../data/genomes/all_viruses.genbank_complete_excllabhost_exclvaccine.191222.unique

#acc_list=../data/metadata/all_viruses.genbank_complete_excllabhost_exclvaccine.191222.sars_filt.accession_only.txt
#out_dir=../data/genomes/all_viruses.genbank_complete_excllabhost_exclvaccine.191222.sars_filt

# For missing files
acc_list=../data/metadata/missing_genomes.accessions_only.txt
out_dir=../data/genomes/all_viruses.220723.filt

# Check rows
cat $acc_list|wc -l

# Parallel download
cat $acc_list \
|xargs \
	-I {} \
	-P 5 \
	ncbi-acc-download \
		--api-key 9728b911e655022461913b0412f82bf09507 \
                --format fasta \
                --out $out_dir/{}.fna \
                {}; \
	sleep 1
