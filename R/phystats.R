# Function to calculate summary statistics from phyloseq object

phystats <- function(phy) {
  return(
    data.frame(
      'Total count in OTU table'=sum(otu_table(phy)),
      'Number of OTUs'=ntaxa(phy),
      'Range of OTU counts'=paste0(range(taxa_sums(phy)), collapse=" - "),
      'Number of singleton OTUs'=length(taxa_sums(phy)[taxa_sums(phy) <= 1]),
      'Number of samples'=nsamples(phy),
      'Range of reads per sample'=paste0(range(sample_sums(phy)), collapse=" - "),
      'Arithmetic mean (±SD) reads per sample'=paste(as.integer(mean(sample_sums(phy))), 
                                                     as.integer(sd(sample_sums(phy))), sep=" ± "),
      'Geometric mean (±SD) reads per sample'=paste(as.integer(exp(mean(log(sample_sums(phy))))),
                                                    as.integer(exp(sd(log(sample_sums(phy))))), sep=" ± "),
      check.names=F)
  )
}