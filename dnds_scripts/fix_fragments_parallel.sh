prokka_dir=../results/dnds_out/prokka_out
out_basedir=../results/dnds_out/agat_fix

ls $prokka_dir|xargs -P 16 -I '{}' sh fix_fragments_single_clique.sh '{}'
