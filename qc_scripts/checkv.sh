db=../databases/checkv-db-v1.5/
#in_path=../data/genomes/all_viruses.220723.filt.formatted.fna
#out_dir=../results/checkv_out

in_path=../data/genomes/all_viruses.220723.filt.missing.formatted.fna
out_dir=../results/checkv_out/missing

checkv end_to_end \
	-t 16 \
	-d $db \
	--remove_tmp \
	--restart \
	$in_path $out_dir

