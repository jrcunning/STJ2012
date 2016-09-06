all: data/otus_97/97_otus.tsv data/otus_97/nw_tophits.tsv data/otus_97/poormatch_IDs.txt data/otus_97_bysample/97_otus_bysample.tsv data/otus_97_bysample/nw_tophits.tsv data/otus_97_bysample/poormatch_IDs.txt data/otus_100/100_otus.tsv data/otus_100/nw_tophits.tsv data/otus_100/poormatch_IDs.txt data/otus_97_byspecies/97_otus_byspecies.tsv data/otus_97_byspecies/nw_tophits.tsv data/otus_97_byspecies/poormatch_IDs.txt

data/otus_97/poormatch_IDs.txt: data/otus_97/nw_tophits.tsv Shell/run_blast_poormatches.sh
	bash Shell/run_blast_poormatches.sh data/otus_97/nw_tophits.tsv data/otus_97/97_otus_rep_set.fasta

data/otus_97/nw_tophits.tsv: data/otus_97/97_otus_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

data/otus_97/97_otus.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_97.sh
	bash Shell/otus_97.sh

data/otus_97_bysample/poormatch_IDs.txt: data/otus_97_bysample/nw_tophits.tsv Shell/run_blast_poormatches.sh
        bash Shell/run_blast_poormatches.sh data/otus_97_bysample/nw_tophits.tsv data/otus_97_bysample/all_rep_set_rep_set.fasta

data/otus_97_bysample/nw_tophits.tsv: data/otus_97_bysample/all_rep_set_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

data/otus_97_bysample/97_otus_bysample.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_97_bysample.sh
	bash Shell/otus_97_bysample.sh

data/otus_97/poormatch_IDs.txt: data/otus_97_byspecies/nw_tophits.tsv Shell/run_blast_poormatches.sh
        bash Shell/run_blast_poormatches.sh data/otus_97_byspecies/nw_tophits.tsv data/otus_97_byspecies/all_rep_set_rep_set.fasta

data/otus_97_byspecies/nw_tophits.tsv: data/otus_97_byspecies/all_rep_set_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

data/otus_97_byspecies/97_otus_byspecies.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_97_byspecies.sh
	bash Shell/otus_97_byspecies.sh

data/otus_100/poormatch_IDs.txt: data/otus_100/nw_tophits.tsv Shell/run_blast_poormatches.sh
        bash Shell/run_blast_poormatches.sh data/otus_100/nw_tophits.tsv data/otus_100/100_otus_rep_set.fasta

data/otus_100/nw_tophits.tsv: data/otus_100/100_otus_rep_set.fasta data/ITS2db_trimmed_derep.fasta R/run_nw.R
	R --vanilla < R/run_nw.R --args $< data/ITS2db_trimmed_derep.fasta

data/otus_100/100_otus.tsv: data/fasta/combined_seqs_trimmed.fasta Shell/otus_100.sh
	bash Shell/otus_100.sh
	
data/ITS2db_trimmed_derep.fasta: data/ITS2db_raw.fasta Shell/prep_ITS2db.sh
	bash Shell/prep_ITS2db.sh
	
data/ITS2db_raw.fasta: data/accn_nos.txt Shell/fetch_seqs.sh
	bash Shell/fetch_seqs.sh

data/fasta/combined_seqs_trimmed.fasta: data/Cunning_3967Raw10232015.zip data/fastq_list.txt Shell/merge_qc_reads.sh
	bash Shell/merge_qc_reads.sh
