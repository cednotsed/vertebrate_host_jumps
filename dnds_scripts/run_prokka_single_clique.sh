genome_dir=../data/genomes/dnds_analysis
out_basedir=../results/dnds_out/prokka_out
out_basedir=./
# Iterate through each clique of interest
clique=$1
clique=Flaviviridae_3
echo $clique
#mkdir $out_basedir/$clique

clique_path=../data/genomes/dnds_analysis/$clique

# Iterate through each genome in clique
for genome_path in $clique_path/*.fna
do 
	prefix=$(echo $genome_path| sed "s|$clique_path/||g"| sed "s|.fna||g")
	echo $prefix	
#	echo $genome_path
	# Run prokka
	prokka \
               	--prefix $prefix \
       	        --kingdom Viruses \
  	        --evalue 1e-03 \
               	--coverage 80 \
		--cpus 1 \
		--outdir $out_basedir/$clique \
		--force \
		--quiet \
	        $genome_path
done

echo "$clique is done!"
