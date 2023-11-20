in_dir=mini_gene_alignments
out_dir=mini_gene_alignments_sorted
#for i in $out_dir/*
#do 
#	echo $i
#	prefix=$(echo $i|sed "s|$out_dir/||g")
#	echo $prefix
#	ls mini_gene_alignments/$prefix-*|xargs -P16 -I @ sh -c "cp @ mini_gene_alignments_sorted/$prefix" x-sh @
#xargs  -I '{}' sh -c "cp ../mini_gene_alignments/'{}'* ./'{}'" x-sh '{}'
#done

for i in $in_dir/*
do
	cp $i $out_dir/$(echo $i|sed "s|$in_dir/||g"|cut -d'-' -f1)
done
