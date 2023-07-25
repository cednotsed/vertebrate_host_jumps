#genome_dir=../data/genomes/viral_family_subsets
#out_dir=../results/mash_out/viral_family_subsets

genome_dir=../data/genomes/source_sink_mini_trees/with_buffer_outgroup
out_dir=../results/mash_out/source_sink_mini_trees/with_buffer_outgroup

n_threads=12

for in_fna in $genome_dir/*.fna
do
	echo $in_fna

	out_path=$(echo $in_fna|sed "s|\\.fna||g"| sed "s|$genome_dir|$out_dir|g")

	echo $out_path

	mash sketch \
		-p $n_threads \
		-i $in_fna \
		-o $out_path \
		-k 13

	mash dist \
		-p $n_threads \
		 ${out_path}.msh ${out_path}.msh \
		-t \
		> ${out_path}.tsv
done
