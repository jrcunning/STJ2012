#!/bin/bash

# Cluster entire dataset at 97% similarity
pick_otus.py -i data/fasta/combined_seqs_trimmed.fasta -s 0.97 -o data/otus_97/

# Remove singletons before picking rep set and assigning taxonomy (to save time)
awk '$3 ~ /./ {print}' data/otus_97/combined_seqs_chimera_filtered_otus.txt > data/otus_97/nosingles_otus.txt

# Get rep set of 97% OTUs
pick_rep_set.py -i data/otus_97/nosingles_otus.txt \
-f data/fasta/combined_seqs_chimera_filtered.fasta \
-o data/otus_97/97_otus_rep_set.fasta