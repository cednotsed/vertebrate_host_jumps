wkdir=../results/dnds_out/vapid_annotation/
in_fna=$wkdir/for_vapid/references.fna
in_sbt=$wkdir/for_vapid/example.sbt
in_meta=$wkdir/for_vapid/dummy_meta.csv
db=~/VAPiD/all_virus.fasta

# Format fasta
#dos2unix $in_fna

python ~/VAPiD/vapid3.py \
	$in_fna \
	$in_sbt \
	--metadata_loc $in_meta \
	--db $db

