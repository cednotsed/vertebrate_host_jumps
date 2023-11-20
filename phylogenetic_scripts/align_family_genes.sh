threads=12

#genomes=../results/dnds_out/family_gene_alignments/SC2-spike.all_genomes.fna

in_dir=../results/dnds_out/family_gene_alignments
for genomes in $in_dir/MERS*all_genomes.fna
do
	alignment_out=$(echo $genomes|sed 's/.fna/.aln/g')

	echo $genomes
	echo $alignment_out

	fftns \
		--anysymbol \
		--adjustdirection \
		--thread ${threads} \
		${genomes} \
		1> ${alignment_out} \
		2> ${alignment_out}.log
done
