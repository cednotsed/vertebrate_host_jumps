kill_dir=../results/dnds_out/kill_lists
in_dir=../results/dnds_out/agat_fix
out_basedir=../results/dnds_out/agat_filt

clique=$1
clique_path=$kill_dir/$clique

clique=$(echo $clique_path|sed "s|$kill_dir/||g")

mkdir $out_basedir/$clique
echo $clique_path

# Parse gffs
for kill_list in $clique_path/*.txt
do
	prefix=$(echo $kill_list|sed "s|$clique_path/||g"|sed "s|.txt||g")
	gff=$in_dir/$clique/$prefix/$prefix.gff
	out_path=$out_basedir/$clique/$prefix.gff
	echo $kill_list
	echo $gff
	echo $prefix
	echo $out_path

	# Run AGAT
	~/AGAT/bin/agat_sp_filter_feature_from_kill_list.pl \
		-gff $gff \
		--kill_list $kill_list \
		-o $out_path
done
