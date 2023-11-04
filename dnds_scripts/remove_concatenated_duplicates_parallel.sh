kill_dir=../results/dnds_out/kill_lists

ls $kill_dir|xargs -P 16 -I '{}' sh remove_concatenated_duplicates_single_clique.sh '{}'
