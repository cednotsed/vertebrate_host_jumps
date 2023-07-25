in_dir=../results/dnds_out/agat_filt
out_basedir=../results/dnds_out/panaroo_out

for clique_path in $in_dir/Anelloviridae_4
do
	clique=$(echo $clique_path|sed "s|$in_dir/||g")
	echo $clique	
	echo $clique_path

	panaroo \
		--threads 4 \
		-i $clique_path/*.gff \
		-o $out_basedir/$clique \
		--clean-mode sensitive \
		--threshold 0.7 \
		--core_threshold 0.95 \
		-a core \
		--aligner mafft
done
