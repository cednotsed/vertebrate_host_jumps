genome_dir=../data/genomes/viral_family_subsets
out_dir=../results/mash_out/viral_family_subsets
n_threads=12

for in_fna in $genome_dir/Anello*.fna
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
