genome_dir=../data/genomes/dnds_analysis
out_basedir=../results/dnds_out/prokka_out
echo $out_basedir

ls $genome_dir|xargs -P 16 -I '{}' sh run_prokka_single_clique.sh '{}'

