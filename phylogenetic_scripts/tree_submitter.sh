aln_dir=/home/zcbtct0/Scratch/vertebrate_host_jumps/data/alignments/source_sink_mini_trees/without_outgroup.masked

for aln_path in $aln_dir/*masked_to*aln
do
	echo $aln_path
	qsub build_trees_with_bootstraps_myriad.sh $aln_path

done
