all: data/otus_97/97_otus.tsv data/otus_97/nw_tophits.tsv data/otus_97/poormatch_IDs.txt data/otus_97_bysample/97_otus_bysample.tsv data/otus_97_bysample/nw_tophits.tsv data/otus_97_bysample/poormatch_IDs.txt data/otus_100/100_otus.tsv data/otus_100/nw_tophits.tsv data/otus_100/poormatch_IDs.txt

# 12. Determine which 97%-across-sample OTUs are probably not Symbiodinium
data/otus_97/poormatch_IDs.txt: data/otus_97/nw_tophits.tsv Shell/run_blast_poormatches.sh
	bash Shell/run_blast_poormatches.sh data/otus_97/nw_tophits.tsv data/otus_97/97_otus_rep_set.fasta

# 11. Assign taxonomy to 97%-across-sample OTUs
data/otus_97/nw_tophits.tsv: data/otus_97/97_otus_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

# 10. Cluster 97%-across-sample OTUs
data/otus_97/97_otus.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_97.sh
	bash Shell/otus_97.sh

# 9. Determine which 97%-within-sample OTUs are probably not Symbiodinium
data/otus_97_bysample/poormatch_IDs.txt: data/otus_97_bysample/nw_tophits.tsv Shell/run_blast_poormatches.sh
	bash Shell/run_blast_poormatches.sh data/otus_97_bysample/nw_tophits.tsv data/otus_97_bysample/all_rep_set_rep_set.fasta

# 8. Assign taxonomy to 97%-within-sample OTUs
data/otus_97_bysample/nw_tophits.tsv: data/otus_97_bysample/all_rep_set_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

# 7. Cluster 97%-within-sample OTUs
data/otus_97_bysample/97_otus_bysample.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_97_bysample.sh
	bash Shell/otus_97_bysample.sh

# 6. Determine which 100% OTUs are probably not Symbiodinium
data/otus_100/poormatch_IDs.txt: data/otus_100/nw_tophits.tsv Shell/run_blast_poormatches.sh
	bash Shell/run_blast_poormatches.sh data/otus_100/nw_tophits.tsv data/otus_100/100_otus_rep_set.fasta

# 5. Assign taxonomy to 100% OTUs
data/otus_100/nw_tophits.tsv: data/otus_100/100_otus_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

# 4. Cluster 100% OTUs
data/otus_100/100_otus.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_100.sh
	bash Shell/otus_100.sh

# 3. Prepare ITS2 reference database
data/ITS2db_trimmed_derep.fasta: data/ITS2db_raw.fasta Shell/prep_ITS2db.sh
	bash Shell/prep_ITS2db.sh

# 2. Download ITS2 reference sequences
data/ITS2db_raw.fasta: data/accn_nos.txt Shell/fetch_seqs.sh
	bash Shell/fetch_seqs.sh

# 1. Merge paired reads for each sample.
data/fasta/combined_seqs_trimmed.fasta: data/Cunning_3967Raw10232015.zip data/fastq_list.txt Shell/merge_qc_reads.sh
	bash Shell/merge_qc_reads.sh
