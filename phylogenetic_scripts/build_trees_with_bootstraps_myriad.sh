#! /usr/bin/bash
#$ -pe smp 12
#$ -l h_rt=48:0:0
#$ -l mem=8G
#$ -cwd
#$ -N supp_tree
source /home/zcbtct0/miniconda3/etc/profile.d/conda.sh
conda activate philodendron

aln_path=$1
#aln_path=../data/alignments/source_sink_mini_trees/without_outgroup.masked/Hepeviridae_4.n113.masked_to_7133pos.aln

output=$(echo $aln_path| sed "s|alignments|trees|g"| sed "s|.aln|.tree|g")
echo $output

iqtree2 \
	-nt 12 \
	-s $aln_path \
	-m MFP \
	--prefix $output \
	-B 1000 \
	--keep-ident \
	-asr

echo $aln_path DONE!!!
