prokka_dir=../results/dnds_out/prokka_out
out_basedir=../results/dnds_out/agat_fix


for clique_path in $prokka_dir/*
do
	clique=$(echo $clique_path|sed "s|$prokka_dir/||g")
	
	mkdir $out_basedir/$clique

	for gff in $clique_path/*.gff
	do
		echo $gff
		prefix=$(echo $gff|sed "s|$clique_path/||g"|sed "s|.gff||g")
		fna=$clique_path/$prefix.fna
		out_dir=$out_basedir/$clique/$prefix
		echo $prefix
		echo $out_dir
		echo $fna

		# Run AGAT
		~/AGAT/bin/agat_sp_prokka_fix_fragmented_gene_annotations.pl \
			-gff $gff \
			--fasta $fna \
			--db ~/miniconda3/envs/prokka/db/kingdom/Viruses/sprot \
			-o $out_dir \
	  	 	--frags
	done
done

#~/AGAT/bin/agat_sp_prokka_fix_fragmented_gene_annotations.pl -gff OQ518334.1.gff --fasta OQ518334.1.fna --db ~/miniconda3/envs/prokka/db/kingdom/Viruses/sprot  -o outfolder2 --frags
