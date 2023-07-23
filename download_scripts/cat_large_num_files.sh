fna_dir=../data/genomes/all_viruses.220723.filt
out_path=../data/genomes/all_viruses.220723.filt.fna

find $fna_dir -maxdepth 1 -type f -name '*.fna' -print0 \
	| sort -zV \
	| xargs -0 cat \
	> $out_path
