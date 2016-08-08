
# Prepare ITS2 reference database
./get_ITS2_database.sh

# Merge, QC, and trim Illumina reads
./merge_qc_reads.sh

# Cluster 97% OTUs within samples
./otus_97.sh

# Assign taxonomy and generate OTU table for 97% OTUs
./tax_results.sh data/otus_97/all_rep_set_rep_set.fasta 97_otus

# Cluster 100% OTUs
./otus_100.sh

# Assign taxonomy and generate OTU table for 97% OTUs
./tax_results.sh data/otus_100/rep_set.fasta 100_otus
