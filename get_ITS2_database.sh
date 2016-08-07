#!/bin/bash

# Get Subtype and Accession number links from old database and use accession number to download sequences from NCBI
cat data/Unaligned_ITS2_Database_31July16.fasta | awk '/>/ {print}' | \
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

# Clean up intermediate files
rm data/ITS2_Database.fasta
rm data/ITS2_Database_inline.fasta
rm data/ITS2_Database_ready.fasta
rm data/ITS2_Database_trimF.fasta
rm data/ITS2_Database_trimF2.fasta
rm data/ITS2_Database_trimF2_trimR.fasta
