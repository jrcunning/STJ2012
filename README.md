[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.803992.svg)](https://doi.org/10.5281/zenodo.803992)

This repository contains all of the raw data, code, and supplemental information for the manuscript:

# Using high-throughput sequencing of ITS2 to describe *Symbiodinium* metacommunities in St. John, U.S. Virgin Islands
### Authors: Ross Cunning, Ruth D. Gates, Peter J. Edmunds
### Journal: _PeerJ_
### Link: [](http://dx.doi.org/)

-----

### Description:
This project investigates *Symbiodinium* communities in scleractinian and milleporine corals at St. John, US Virgin Islands using high-throughput sequencing of *Symbiodinium* ITS2 nrDNA. The major aims of this work are to:

1. Apply a novel bioinformatic approach (within-sample clustering) to better understand *Symbiodinium* ecology using ITS2 metabarcoding data

2. Characterize *Symbiodinium* communities in a range of coral species using high-throughput sequencing (many for the first time)

The bioinformatic pipeline\* is coded using a Makefile that calls a series of Shell and R scripts that generate the reference database, *Symbiodinium* OTU tables, and taxonomic assignments from the raw sequence data. Downstream statistical analysis and figure generation is then performed in the **R/analysis.R** script.

*\*Code for the bioinformatic pipeline utilized here has been generalized to be readily used on new datasets at [http://github.com/jrcunning/SymITS2](http://github.com/jrcunning/SymITS2)* 

### Repository contents:

#### Makefile: 
- This Makefile codes the entire bioinformatic pipeline utilized in this study with calls to scripts in the **Shell/** and **R/** directories. Given the raw sequence reads (**data/Cunning_3967Raw10232015.zip**), the list of samples contained therein (**data/fastq_list.txt**), and the list of accession numbers to include in the reference database (**data/accn_nos.txt**), execution of this Makefile using `make` will generate OTUs, build the ITS2 reference database, and assign taxonomy using 100% clustering, 97% clustering, and 97%-within-sample clustering.

#### Shell/:
- **fetch_seqs.sh:** This script downloads sequence data from NCBI given a list of named sequences with accession numbers (e.g., >C1_AF333515)
- **merge_qc_reads.sh:** This script unzips the fastq archive, merges and QC's reads, filters out chimeras, and trims primer sequences.
- **otus_100.sh:** This script clusters ITS2 sequences at 100% identity across all samples.
- **otus_97.sh:** This script clusters ITS2 sequences at 97% identity across all samples.
- **otus_97_bysample.sh:** This script separates ITS2 sequences by sample and clusters at 97% identity within each sample independently.
- **prep_ITS2db.sh:** This script prepares the ITS2 reference database from downloaded sequences by aligning sequences within clades, trimming to equal length, and removing duplicated sequences.
- **run_blast_poormatches.sh:** This script identifies OTU representative sequences with poor matches to the ITS2 reference database (<90% similarity to best match) and BLASTs these sequences against the NCBI nt database.

#### R/:
- **analysis.R:** This script codes all of the analysis downstream of OTU clustering and taxonomic assignment, including OTU quality filtering, summary statistics, and analysis of *Symbiodinium* community composition using multivariate statistics and network analysis. This script also generates all of the figures presented in the published manuscript (**figures/Fig\*.png**).
- **functions.R:** This script contains several accessory functions that are utilized in the analysis.R script.
- **run_nw.R:** This script assigns taxonomy to all OTUs by performing Needleman-Wunsch global alignments of representative sequences against the ITS2 reference database.

#### data/:
Primary data:

- **Cunning_3967Raw20232015.zip:** Archive of raw sequencing reads in fastq format.
- **coastline/:** Contains data downloaded from the [NOAA GSHHG database](https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html) used to create map figures.
- **EDNAFULL:** Nucleotide substitution matrix utilized in global alignment for taxonomic assignments.
- **STJ_envdata2.csv:** Contains SST and chlorophyll a data for north and south shores of St. John
- **accn_nos.txt:** The list of named *Symbiodinium* ITS2 sequences with accession numbers (e.g., >C1_AF333515) to be included in the reference database.
- **fastq_list.txt:** The names of the fastq files in which paired reads for each sample can be found (needed for merging paired reads).
- **mapping_file.txt:** Metadata (e.g. species, collection site) for each sample in the dataset.

Intermediate data (these files are all created by executing `make`):

- **ITS2db_trimmed.fasta:** The trimmed ITS2 reference database
- **ITS2db_trimmed_derep.fasta::** The trimmed and dereplicated ITS2 reference database
- **ITS2db_trimmed_notuniques_otus.txt:** A list of the non-unique sequences in the ITS2 reference database
- **otus_100/:** OTUs, representative sequences, and taxonomic assignments generated by 100% clustering
- **otus_97/:** OTUs, representative sequences, and taxonomic assignments generated by 97% clustering across samples
- **otus_97_bysample/:** OTUs, representative sequences, and taxonomic assignments generated by 97%-within-sample-clustering

#### figures/:
- **Fig\*.png:** Each of the figures presented in the published manuscript. The code that produces each of these figures is in **R/analysis.R**.

#### supp/:
- **STJ2012_supp.Rmd:** R Markdown document containing all supplemental information and analysis.
- **STJ2012_supp.html:** Output of supplemental R Markdown document.

### Software dependencies:
The following software is utilized by components of this code: 

- [GNU parallel](https://www.gnu.org/software/parallel/)
- [illumina-utils](https://github.com/merenlab/illumina-utils)
- [oligotyping](https://github.com/merenlab/oligotyping)
- [usearch](http://drive5.com/usearch/)
- [cutadapt](https://github.com/marcelm/cutadapt)
- [QIIME](http://qiime.org)
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [R](https://www.r-project.org) (see session info below)

~~~~
R > sessionInfo()

R version 3.3.1 (2016-06-21)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.12.3 (Sierra)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     methods   base     

other attached packages:
 [1] Biostrings_2.40.2   XVector_0.12.1      IRanges_2.6.1       S4Vectors_0.10.3   
 [5] BiocGenerics_0.18.0 pBrackets_1.0       RColorBrewer_1.1-2  stringr_1.1.0      
 [9] igraph_1.0.1        reshape2_1.4.1      multcompView_0.1-7  vegan_2.4-1        
[13] lattice_0.20-34     permute_0.9-4       rgdal_1.1-10        sp_1.2-3           
[17] phyloseq_1.16.2
~~~~

