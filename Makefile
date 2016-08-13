data/97_otus_tax_assignments.txt data/97_otus.tsv: data/otus_97/97_otus_rep_set.fasta data/ITS2db_trimmed.fasta
	bash tax_results.sh $< 97_otus

data/97_otus_bysample_tax_assignments.txt data/97_otus_bysample.tsv: data/otus_97_bysample/all_rep_set_rep_set.fasta data/ITS2db_trimmed.fasta
	bash tax_results.sh $< 97_otus_bysample
	
data/100_otus_tax_assignments.txt data/100_otus.tsv: data/otus_100/100_otus_rep_set.fasta data/ITS2db_trimmed.fasta
	bash tax_results.sh $< 100_otus
	
data/otus_97/97_otus_rep_set.fasta: data/fasta/combined_seqs_trimmed.fasta
	bash otus_97.sh
	
data/otus_100/100_otus_rep_set.fasta: data/fasta/combined_seqs_trimmed.fasta
	bash otus_100.sh

data/otus_97_bysample/all_rep_set_rep_set.fasta: data/fasta/combined_seqs_trimmed.fasta
	bash otus_97_bysample.sh
	
data/ITS2db_trimmed.fasta: data/accn_nos.txt
	bash get_ITS2_database.sh

#data/fasta/combined_seqs_trimmed.fasta: data/Cunning_3967Raw10232015.zip data/fastq_list.txt
#	bash merge_qc_reads.sh
