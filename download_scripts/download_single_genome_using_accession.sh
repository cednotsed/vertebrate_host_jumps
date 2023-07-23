acc_list=../data/metadata/all_viruses.220723.filt.accessions_only.txt
out_dir=../data/genomes/all_viruses.220723.filt

# Check rows
cat $acc_list|wc -l

ncbi-acc-download \
	--api-key 9728b911e655022461913b0412f82bf09507 \
        --format fasta \
        --out $out_dir/{}.fna \
        KY486145.1

