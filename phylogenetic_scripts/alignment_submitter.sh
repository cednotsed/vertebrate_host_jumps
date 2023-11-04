genome_dir=/home/zcbtct0/Scratch/vertebrate_host_jumps/data/genomes/source_sink_mini_trees/full_alignments

for genomes in $genome_dir/*.fna
do
	echo $genomes
	qsub align_multiple_genomes_mafft.sh $genomes
done
