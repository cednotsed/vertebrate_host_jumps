ori_dir=../data/genomes/all_viruses.220723.filt
fna_dir=../data/genomes/missing
out_path=../data/genomes/all_viruses.220723.filt.missing.fna

while read fna 
do
	cp $ori_dir/${fna}.fna $fna_dir

done < ../data/metadata/missing_genomes.accessions_only.txt

find $fna_dir -maxdepth 1 -type f -name '*.fna' -print0 \
	| sort -zV \
	| xargs -0 cat \
	> $out_path
