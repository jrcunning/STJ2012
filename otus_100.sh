#!/bin/bash

# Cluster entire dataset at 100% similarity
pick_otus.py -i data/fasta/combined_seqs_chimera_filtered.fasta -s 1.0 --optimal_uclust -o data/otus_100/

# Remove singletons before picking rep set and assigning taxonomy (to save time)
awk '$3 ~ /./ {print}' data/otus_100/combined_seqs_chimera_filtered_otus.txt > data/otus_100/nosingles_otus.txt

# Get rep set of 100% OTUs
pick_rep_set.py -i data/otus_100/nosingles_otus.txt \
-f data/fasta/combined_seqs_chimera_filtered.fasta \
-o data/otus_100/rep_set.fasta

# Assign taxonomy
assign_taxonomy.py -i data/otus_100/rep_set.fasta -m blast \
-r data/Unaligned_ITS2_Database_31July16.fasta \
-t data/id_to_taxonomy_31July16.txt \
-o data/otus_100/blast_taxonomy

# Make list of "no blast hits"
awk '/No blast hit/' data/otus_100/blast_taxonomy/rep_set_tax_assignments.txt > data/otus_100/blast_taxonomy/no_blast_hits.txt

# Make OTU table excluding no blast hits
make_otu_table.py -i data/otus_100/nosingles_otus.txt \
-t data/otus_100/blast_taxonomy/rep_set_tax_assignments.txt \
-e data/otus_100/blast_taxonomy/no_blast_hits.txt \
-o data/otus_100/100_otu_table.biom
rm data/100_otu_table.tsv  # Delete old OTU .tsv if present
biom convert -i data/otus_100/100_otu_table.biom -o data/100_otu_table.tsv --to-tsv

# Add clade to tax table and copy to data directory
cut -d$'\t' -f4- data/otus_100/blast_taxonomy/rep_set_tax_assignments.txt | cut -c-1 | paste - data/otus_100/blast_taxonomy/rep_set_tax_assignments.txt > data/100_tax_assignments.txt
