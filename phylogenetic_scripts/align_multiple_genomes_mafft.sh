#! /usr/bin/bash
#$ -pe smp 8
#$ -l h_rt=48:0:0
#$ -l mem=8G
#$ -cwd
#$ -N alignment
source /home/zcbtct0/miniconda3/etc/profile.d/conda.sh
conda activate philodendron

genomes=$1

threads=8
alignment_out=$(echo $genomes|sed 's/.fna/.aln/g'|sed "s|genomes|alignments|g")

echo $genomes
echo $alignment_out

fftns \
	--reorder \
	--anysymbol \
	--adjustdirection \
	--thread ${threads} \
	${genomes} \
	1> ${alignment_out} \
	2> ${alignment_out}.log
