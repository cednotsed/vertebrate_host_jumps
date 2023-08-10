threads=12

genomes=../results/dnds_out/family_gene_alignments/CoV2_spike.fna

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
