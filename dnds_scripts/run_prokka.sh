genome_dir=../data/genomes/dnds_analysis

out_basedir=../results/dnds_out/prokka_out
echo $out_basedir

# Iterate through each clique of interest
for clique_path in $genome_dir/*
do
	clique=$(echo $clique_path|sed "s|$genome_dir/||g")
	echo $clique
	mkdir $out_basedir/$clique
	
	# Iterate through each genome in clique
	for genome_path in $clique_path/*.fna
	do 
		prefix=$(echo $genome_path| sed "s|$clique_path/||g"| sed "s|.fna||g")
		echo $prefix	
		# Run prokka
		prokka \
	               	--prefix $prefix \
        	        --kingdom Viruses \
       	        	--evalue 1e-09 \
	               	--coverage 80 \
			--cpus 8 \
			--outdir $out_basedir/$clique \
			--force \
	               	$genome_path
	done
done

