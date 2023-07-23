db=../databases/checkv-db-v1.2/
in_path=../data/genomes/assembly_scaffolds/pooled_contigs_gt500.MHV_filt.fna
out_dir=../results/checkv_out

checkv end_to_end \
	-t 16 \
	-d $db \
	--remove_tmp \
	--restart \
	$in_path $out_dir

