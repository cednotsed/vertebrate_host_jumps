in_path=../results/dnds_out/ncbi_annotation/sequence.chunk_6.gff3
out_path=$(echo $in_path|sed "s|.gff3|.parsed.gff3|g")

dos2unix $in_path

cat $in_path| \
	grep -v '^#'| \
	awk 'NF' \
	> $out_path
