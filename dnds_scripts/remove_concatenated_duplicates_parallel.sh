kill_dir=../results/dnds_out/kill_lists
in_dir=../results/dnds_out/agat_fix
out_basedir=../results/dnds_out/agat_filt

ls $kill_dir|xargs -P 4 -I '{}' sh remove_concatenated_duplicates_single_clique.sh '{}'
