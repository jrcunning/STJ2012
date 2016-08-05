#!/bin/bash

# Unzip sequence data archive
unzip -j data/Cunning_3967Raw10232015.zip -d data/fastq

# Generate config files for merging reads with illumina-utils
mkdir -p data/merge
cp data/fastq_list.txt data/fastq/fastq_list.txt
iu-gen-configs data/fastq/fastq_list.txt -o data/merge

# Merge paired reads for each sample (in parallel)
parallel iu-merge-pairs {} --min-overlap-size 150 --enforce-Q30-check --marker-gene-stringent ::: data/merge/*.ini

# Filter sequences - keep only those with 3 mismatches or less (in parallel)
parallel iu-filter-merged-reads {} --max-mismatches 3 ::: data/merge/*_MERGED

# Add QIIME labels
add_qiime_labels.py -m data/mapping_file.txt -i data/merge/ -c InputFileName -o data/fasta

# Chimera checking
identify_chimeric_seqs.py -i data/fasta/combined_seqs.fna --suppress_usearch61_ref -m usearch61 -o data/fasta/usearch61_chimeras
filter_fasta.py -f data/fasta/combined_seqs.fna -o data/fasta/combined_seqs_chimera_filtered.fasta -s data/fasta/usearch61_chimeras/chimeras.txt -n
