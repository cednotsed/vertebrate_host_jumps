#in_fna=./test.fna
#out_fna=./test.formatted.fna

in_fna=../data/genomes/all_viruses.220723.filt.fna
out_fna=../data/genomes/all_viruses.220723.filt.formatted.fna

dos2unix $in_fna

# Remove empty lines
awk 'NF' $in_fna | \
# Replace weird characters with Ns
sed '/^>/! s/[^ACTG\r\n]/N/g' > $out_fna
