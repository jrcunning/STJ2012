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
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 data/ITS2_Database_trimF2_trimR.fasta -o data/ITS2db_trimmed.fasta

# Remove Ns if present in sequences -- noting that several were deposited in GenBank using Ns as gaps...
awk '!/>/ { gsub("N","") }; { print $0 }' data/ITS2db_trimmed.fasta > data/ITS2db.fasta


# Generate id_to_taxonomy file
TAB=$'\t'  # ( because \t is not recognized by sed on mac os x )
sed -e 's/>\([A-Z]\)\(.*\)\(_[A-Z][A-Z][0-9][0-9][0-9][0-9][0-9][0-9]\)/\1\2\3'"${TAB}"'Symbiodiniaceae;Symbiodinium;Clade\1;\1\2;_;_/' -e 'tx' -e 'd' -e ':x' data/ITS2db.fasta > data/id_to_taxonomy.txt


awk '/>A/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeA_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeA_align_trimE.fasta data/cladeA_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeA_align_trimES.fasta data/cladeA_align_trimE.fasta

awk '/>B/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeB_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeB_align_trimE.fasta data/cladeB_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeB_align_trimES.fasta data/cladeB_align_trimE.fasta

awk '/>C/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeC_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeC_align_trimE.fasta data/cladeC_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeC_align_trimES.fasta data/cladeC_align_trimE.fasta

awk '/>D/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeD_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeD_align_trimE.fasta data/cladeD_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeD_align_trimES.fasta data/cladeD_align_trimE.fasta

awk '/>E/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeE_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeE_align_trimE.fasta data/cladeE_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeE_align_trimES.fasta data/cladeE_align_trimE.fasta

awk '/>F/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeF_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeF_align_trimE.fasta data/cladeF_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeF_align_trimES.fasta data/cladeF_align_trimE.fasta

awk '/>G/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeG_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeG_align_trimE.fasta data/cladeG_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeG_align_trimES.fasta data/cladeG_align_trimE.fasta

awk '/>H/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeH_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeH_align_trimE.fasta data/cladeH_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeH_align_trimES.fasta data/cladeH_align_trimE.fasta

awk '/>I/ {print; getline; print}' data/ITS2db.fasta | muscle -out data/cladeI_align.fasta
o-smart-trim --min-percent 80 -E -o data/cladeI_align_trimE.fasta data/cladeI_align.fasta
o-smart-trim --min-percent 80 -S -o data/cladeI_align_trimES.fasta data/cladeI_align_trimE.fasta

cat data/*_trimES.fasta | sed 's/-//g' > data/ITS2db_trimmed.fasta



# Clean up intermediate files
rm -r data/dbseqs
rm data/ITS2_Database.fasta
rm data/ITS2_Database_inline.fasta
rm data/ITS2_Database_ready.fasta
rm data/ITS2_Database_trimF.fasta
rm data/ITS2_Database_trimF2.fasta
rm data/ITS2_Database_trimF2_trimR.fasta
rm data/clade*
