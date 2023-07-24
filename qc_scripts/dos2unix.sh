fna_dir=../data/genomes/viral_family_subsets

for fna in $fna_dir/*.fna
do
	dos2unix $fna
done
