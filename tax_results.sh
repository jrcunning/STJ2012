#!/bin/bash

DIR=dirname $1


# Assign taxonomy
assign_taxonomy.py -i $1 -r data/ITS2db.fasta -t data/id_to_taxonomy.txt -m blast -o $DIR/blast_taxonomy

# Make list of "no blast hits"
awk '/No blast hit/' $DIR/blast_taxonomy/*_tax_assignments.txt > $DIR/blast_taxonomy/no_blast_hits.txt

# Make OTU table excluding no blast hits
make_otu_table.py -i $DIR/nosingles_otus.txt \
-t $DIR/blast_taxonomy/all_rep_set_rep_set_tax_assignments.txt \
-e $DIR/blast_taxonomy/no_blast_hits.txt -o $DIR/$2.biom
rm $DIR/*.tsv  # delete old OTU .tsv if present
biom convert -i $DIR/$2.biom -o data/$2.tsv --to-tsv

# Add clade to tax table and copy to data directory
cut -d$'\t' -f4- $DIR/blast_taxonomy/*_tax_assignments.txt | cut -c-1 | paste - $DIR/blast_taxonomy/*_tax_assignments.txt > data/$2_tax_assignments.txt
