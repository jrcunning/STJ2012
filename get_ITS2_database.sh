#!/bin/bash

awk '/>/ {print}' data/Unaligned_ITS2_Database_31July16.fasta > data/seqIDs1.txt
sed 's/\.1$//' data/seqIDs1.txt > data/seqIDs.txt  # remove trailing '.1' on some accession numbers
cat data/addseqs.txt >> data/seqIDs.txt

# Get Subtype and Accession number links from old database and use accession number to download sequences from NCBI
cat data/seqIDs.txt | \
while read n
do
	acn=${n##*_}
	echo $n
	curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acn}&rettype=fasta"
done > data/ITS2_Database.fasta
	
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
cutadapt -a AAGCATATAAGTAAGCGGAGG -e 0.15 data/ITS2_Database_trimF2_trimR.fasta -o data/ITS2db.fasta

# Generate id_to_taxonomy file
TAB=$'\t'  # ( because \t is not recognized by sed on mac os x )
sed -e 's/>\([A-Z]\)\(.*\)\(_[A-Z][A-Z][0-9][0-9][0-9][0-9][0-9][0-9]\)/\1\2\3'"${TAB}"'Symbiodiniaceae;Symbiodinium;Clade\1;\1\2;_;_/' -e 'tx' -e 'd' -e ':x' data/ITS2db.fasta > data/id_to_taxonomy.txt

# Clean up intermediate files
rm data/seqIDs1.txt
rm data/seqIDs.txt
rm data/ITS2_Database.fasta
rm data/ITS2_Database_inline.fasta
rm data/ITS2_Database_ready.fasta
rm data/ITS2_Database_trimF.fasta
rm data/ITS2_Database_trimF2.fasta
rm data/ITS2_Database_trimF2_trimR.fasta
