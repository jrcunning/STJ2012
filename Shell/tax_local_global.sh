#!/bin/bash

dir=$(dirname $1)

### Get blast results (local alignment) for each sequence
formatdb -p F -i data/ITS2db_trimmed.fasta
blastall -p blastn -i $1 -d data/ITS2db_trimmed.fasta -b 1 -m 8 > $dir/blast_results.txt
# Sort by query sequence and then E-value in descending order, then get only the first hist for each query sequence
sort -k1,1 -k11,11nr $dir/blast_results.txt | sort -u -k1,1 > $dir/blast_results_filtered.txt
# cleanup 
rm -f data/ITS2db_trimmed.fasta.*

### Get top needle hit (global alignment) for each sequence
# split fasta file into separate file for each sequence to feed to get_needle
rm -rf $dir/needle; mkdir $dir/needle
split -l 2 $1 $dir/needle/
# rename files according to OTU name
for i in $dir/needle/*
do
d="$(head -1 "$i" | cut -d " " -f1 | cut -c 2-).fasta";
mv "$i" $dir/needle/$d
done
# build function to get top needle hit from a sequence in above format
get_needle() {
	needle -gapopen 10.0 -gapextend 0.5 -endopen 0.0 -endextend 0.0 -asequence $1 -bsequence data/ITS2db_trimmed.fasta -outfile $1.out -aformat score 
}
export -f get_needle
# feed each line of tab-separated sequences to the get_needle function
parallel get_needle :::  $dir/needle/*.fasta
# get top hit for each sequence and write to file
parallel "sort -nt'(' -k2 {} | tail -n 1" ::: $dir/needle/*.out | sort -k1,1 > $dir/needle_results.txt


## merge together blast results and needle results to compare
join $dir/blast_results_filtered.txt $dir/needle_results.txt > $dir/tax_results.txt