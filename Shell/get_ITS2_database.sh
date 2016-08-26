#!/bin/bash

# use accession numbers to download sequences from NCBI in parallel
getseq() {
	n=$1
	acn=${n##*_}
	echo $1 > data/dbseqs/${acn}.fasta
	curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acn}&rettype=fasta" >> data/dbseqs/${acn}.fasta
}
export -f getseq
rm -rf data/dbseqs; mkdir data/dbseqs
parallel -a data/accn_nos.txt getseq 
cat data/dbseqs/*.fasta > data/ITS2_Database.fasta

	
# Reformat names of new database to be only accession numbers & subtypes and put sequences on one line
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' data/ITS2_Database.fasta > data/ITS2_Database_inline.fasta
awk '!/^>gi/ {print}' data/ITS2_Database_inline.fasta > data/ITS2_Database_ready.fasta

# Download sequences from supplementary material of Green et al. 2014 and add to database
curl https://peerj.com/articles/386/DataS1_PerlScripts_alignment_Bioinformatics.zip > data/green.zip
unzip -p data/green.zip DataS1* > data/greenseqs.aln
awk '/Haplotype/ {a[$1] = a[$1]"\n"$2}END{for(i in a){print ">B1."i""a[i]}}' data/greenseqs.aln | sed 's/-//g' | \
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | \
sed 's/Haplotype\(.*$\)/\1_Green2014/' > data/greenseqs.fasta
rm data/green.zip
rm data/greenseqs.aln
cat data/greenseqs.fasta >> data/ITS2_Database_ready.fasta

# Trim primers from database sequences using cutadapt
# Trim forward primers using cutadapt
#   Allow error rate of 15% (0 indels/mismatches)
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 data/ITS2_Database_ready.fasta -o data/ITS2_Database_trimF.fasta
# Trim forward primers again (may be multiple internal primer sequences), but do not discard sequences that do not contain primer
cutadapt -g GTGAATTGCAGAACTCCGTG -e 0.15 data/ITS2_Database_trimF.fasta -o data/ITS2_Database_trimF2.fasta
# Trim reverse primers using cutadapt
#   Do not remove sequences that do not have reverse primer sequence
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 data/ITS2_Database_trimF2.fasta -o data/ITS2_Database_trimF2_trimR.fasta
# Trim reverse primers again for any remaining internally
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 data/ITS2_Database_trimF2_trimR.fasta -o data/ITS2db_trimFR.fasta

# Remove Ns if present in sequences -- noting that several were deposited in GenBank using Ns as gaps...
awk '!/>/ { gsub("N","") }; { print $0 }' data/ITS2db_trimFR.fasta > data/ITS2db.fasta


# Generate id_to_taxonomy file
TAB=$'\t'  # ( because \t is not recognized by sed on mac os x )
sed -e 's/>\([A-Z]\)\(.*\)\(_.*$\)/\1\2\3'"${TAB}"'Symbiodiniaceae;Symbiodinium;Clade\1;\1\2;_;_/' -e 'tx' -e 'd' -e ':x' data/ITS2db.fasta > data/id_to_taxonomy.txt


awk '/>A/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeA_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeA_align_trimE.fasta data/cladeA_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeA_align_trimES.fasta data/cladeA_align_trimE.fasta

awk '/>B/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeB_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeB_align_trimE.fasta data/cladeB_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeB_align_trimES.fasta data/cladeB_align_trimE.fasta

awk '/>C/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeC_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeC_align_trimE.fasta data/cladeC_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeC_align_trimES.fasta data/cladeC_align_trimE.fasta

awk '/>D/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeD_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeD_align_trimE.fasta data/cladeD_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeD_align_trimES.fasta data/cladeD_align_trimE.fasta

awk '/>E/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeE_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeE_align_trimE.fasta data/cladeE_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeE_align_trimES.fasta data/cladeE_align_trimE.fasta

awk '/>F/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeF_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeF_align_trimE.fasta data/cladeF_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeF_align_trimES.fasta data/cladeF_align_trimE.fasta

awk '/>G/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeG_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeG_align_trimE.fasta data/cladeG_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeG_align_trimES.fasta data/cladeG_align_trimE.fasta

awk '/>H/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeH_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeH_align_trimE.fasta data/cladeH_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeH_align_trimES.fasta data/cladeH_align_trimE.fasta

awk '/>I/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeI_align.fasta
o-smart-trim --min-percent 90 -E -o data/cladeI_align_trimE.fasta data/cladeI_align.fasta
o-smart-trim --min-percent 90 -S -o data/cladeI_align_trimES.fasta data/cladeI_align_trimE.fasta

cat data/*_trimES.fasta | sed 's/-//g' > data/ITS2db_trimmed.fasta



# Clean up intermediate files
rm -r data/dbseqs
rm data/ITS2_Database.fasta
rm data/ITS2_Database_inline.fasta
rm data/ITS2_Database_ready.fasta
rm data/ITS2_Database_trimF.fasta
rm data/ITS2_Database_trimF2.fasta
rm data/ITS2_Database_trimF2_trimR.fasta
rm data/ITS2db_trimFR.fasta
rm data/clade*
rm data/greenseqs.fasta


### Identify identical reference sequences (post-trimming)
# Get database sequences that are not uniqueq
usearch -derep_fulllength data/ITS2db_trimmed.fasta -output data/ITS2db_trimmed_notunique.fasta -sizeout -minuniquesize 2
# Get fasta file of all the non-unique sequences
usearch -usearch_global data/ITS2db_trimmed.fasta -db data/ITS2db_trimmed_notunique.fasta -strand plus -id 1.0 \
-matched data/ITS2db_trimmed_notuniques.fasta
# Re-cluster non-unique sequences at 100% identity to get list of which ones are identical
pick_otus.py -i data/ITS2db_trimmed_notuniques.fasta -s 1.0 -o data
# list of identical sequence groups is in output file: data/ITS2db_trimmed_notuniques_otus.txt

# Clean up
rm data/ITS2db_trimmed_notunique.fasta
rm data/ITS2db_trimmed_notuniques_clusters.uc
rm data/ITS2db_trimmed_notuniques.fasta
rm data/ITS2db_trimmed_notuniques_otus.log
