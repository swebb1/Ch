
## Make sure contig names are unique to each genome, e.g. chr1 and 1, or rename spike as chr1_spike etc.

## Combine reference and index with BWA
ln -s /homes/genomes/chicken/GRCg7b/GRCg7b.fa .
ln -s /homes/genomes/d.melanogaster/dm6/dm6.fa .
cat GRCg7b.fa dm6.fa > gg7b_dm6.fa 
bwa index -a bwtsw -p gg7b_dm6 gg7b_dm6.fa

## Get chromosome files
cut -f 1,2 /homes/genomes/chicken/GRCg7b/GRCg7b.fa.fai | perl -lane '$e=$F[1]+1000; print "$F[0]\t0\t$e";'  > gg7b.bed
cut -f 1,2 /homes/genomes/d.melanogaster/dm6/dm6.fa.fai| perl -lane '$e=$F[1]+1000; print "$F[0]\t0\t$e";' >  dm6.bed

## Length files for genome and ucsc converted (chr)
awk 'BEGIN {OFS="\t"}; { print $1, $3 }' gg7b.bed > gg7b.len  
grep -v "^N" gg7b.bed | awk 'BEGIN {OFS="\t"}; { print $1, "chr"$1 }' > chrmap.txt 
perl ../bin/scripts/chrMap.pl --in gg7b.len --chr chrmap.txt > gg7b.ucsc.len 

