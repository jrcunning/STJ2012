#!/bin/bash

dir=$(dirname $1)

### Get blast results (local alignment) for each sequence
formatdb -p F -i $2
blastall -p blastn -i $1 -d $2 -b 1 -m 8 > $dir/blast_results.txt
# Sort by query sequence and then E-value in descending order, then get only the first hit for each query sequence
sort -k1,1 -k11,11nr $dir/blast_results.txt | sort -u -k1,1 > $dir/blast_tophits.txt
# cleanup 
rm -f $2.*
