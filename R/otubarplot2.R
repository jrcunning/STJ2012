otubarplot2 <- function(samples) {
  # Get phyloseq objects for selected samples from each clustering approach
  samples <<- samples
  phy100 <- subset_samples(phy100.f.p, rownames(sample_data(phy100.f.p)) %in% samples)
  phy100 <- filter_taxa(phy100, function(x) sum(x)!=0, TRUE)
  phy97 <- subset_samples(phy97.f.p, rownames(sample_data(phy97.f.p)) %in% samples)
  phy97 <- filter_taxa(phy97, function(x) sum(x)!=0, TRUE)
  phy97bs <- subset_samples(phy97bs.f.p, rownames(sample_data(phy97bs.f.p)) %in% samples)
  phy97bs <- filter_taxa(phy97bs, function(x) sum(x)!=0, TRUE)
  
  # Get sample data
  samdat <- data.frame(sample_data(phy100))
  # Set sample order to go by dominant 97% OTU symbiont subtype then coral species
  sord <- rownames(sample_data(phy97)[order(factor(tax_table(phy97)[apply(otu_table(phy97), MARGIN=2, FUN=which.max), "Subtype"]), sample_data(phy97)$Species), ])
  samdat <- samdat[sord, ]
  # Get tax table and otu table
  taxdat <- data.frame(tax_table(phy100))
  otutab <- as.matrix(otu_table(phy100)[, rownames(samdat)])

  # Get top n most abundant sequence variants and unique colors for each one
  ranks <- order(rowSums(otutab), decreasing=T)
  n <- 45
  topn <- rownames(otutab[ranks[1:n], ])
  topn <- data.frame(tax_table(phy100)[topn, "Clade"])
  # Get 'qualitative' color palettes from color brewer to assign to each sequence variant
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',][-c(3:5), ]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colors <- sample(col_vector, n)
  topn$colors <- colors
  # Plot bars
  par(mgp=c(0,0.2,0), cex.axis=0.6, mar=c(6,2,3,1))
  bars <- barplot(sort(otutab[, 1], decreasing=T), width=1/(ncol(otutab)*1), xlim=c(0, 1),
                  col=topn[rownames(otutab)[order(otutab[, 1], decreasing=T)], "colors"], las=1, space=0, xaxs="i",
                  xlab="", ylab="", xaxt="n", yaxt="n")
  axis(side=2, tcl=-0.15, at=c(0.2, 0.4, 0.6, 0.8, 1.0), las=1)
  for (i in 2:ncol(otutab)) {
    barplot(sort(otutab[, i], decreasing=T), add=T, yaxt="n", xaxt="n", width=1/(ncol(otutab)*1), space=1*(i-1),
            col=topn[rownames(otutab)[order(otutab[, i], decreasing=T)], "colors"])
  }
  # Get midpoint, start, and end x-coords of each bar for annotating plot
  midpts <- seq(from=bars, by=bars*2, length.out=nsamples(phy100))
  starts <- midpts - bars
  ends <- midpts + bars
  # Annotate plot
  text(-0.03, 1.18, label=expression(paste(bold("A. "), "Relative abundance of unique sequence variants in samples:", sep=" ")), 
       xpd=NA, adj=c(0,0))
  # Add species names
  spnames <- as.character(sample_data(phy97)$Species[match(sord, rownames(sample_data(phy97)))])
  fspnames <- c("D. strigosa", "M. alcicornis", "P. furcata", "O. annularis", "S. siderea", "F. fragum", "S. radians", "P. astreoides", "D. cylindrus", "M. cavernosa")
  spnames <- sapply(spnames, function(x) str_subset(fspnames, x))
  newsp <- c(1, which(diff(as.numeric(factor(spnames)))!=0)+1, nsamples(phy100)+1)
  i <- 1
  while (i < length(newsp)) {
    brackets(starts[newsp[i]], 1, ends[newsp[i+1]-1], 1, xpd=NA)
    text((starts[newsp[i]]+ends[newsp[i+1]-1])/2, 1.1, labels=substitute(italic(x), list(x=spnames[newsp[i]])), xpd=NA, cex=0.8)
    i <- i + 1
  }
  # Add dominant OTU assigned by 97% clustering across samples
  text(-0.03, -0.125, label=expression(paste(bold("B. "), "Dominant OTU assigned by 97% clustering across samples:")), xpd=NA, adj=c(0,0))
  names97 <- unlist(lapply(strsplit(as.character(data.frame(tax_table(phy97)[apply(data.frame(otu_table(phy97), check.names=F)[,sord], MARGIN=2, FUN=which.max), "Subtype"])$Subtype), split="_"), "[[", 1))
  otus97 <- rownames(tax_table(phy97)[apply(data.frame(otu_table(phy97), check.names=F)[,sord], MARGIN=2, FUN=which.max), ])
  newotus97 <- c(1, which(diff(as.numeric(factor(otus97)))!=0)+1, nsamples(phy100)+1)
  # Get color of dominant 97% OTU representative sequence, match to color of representative unique sequence variant
  # Get representative sequence for dominant 97% OTU
  repseqs97 <- sapply(otus97, function(dom97otu) {
    system(sprintf("awk '/%s / {print $2}' data/otus_97/97_otus_rep_set.fasta", dom97otu), intern=T)
  })
  # Find 100% OTU to which representative sequence for dominant 97% OTU belongs
  equivotu100 <- sapply(repseqs97, function(repseq) {
    system(sprintf("awk '/%s\t/ {print $1}' data/otus_100/nosingles_otus.txt", repseq), intern=T)
  })
  # Get matching color for dominant dominant 97% OTU
  matchcols <- topn[equivotu100, "colors"]
  # Plot rectangles with appropriate fill colors for 97% OTUs
  i <- 1
  while (i < length(newotus97)) {
    polygon(x=c(rep(starts[newotus97[i]], 2), rep(ends[newotus97[i+1]-1], 2)), y=c(-0.15, -0.25, -0.25, -0.15), 
            col=unique(matchcols)[i], xpd=NA)
    text((starts[newotus97[i]]+ends[newotus97[i+1]-1])/2, -0.2, labels=names97[newotus97[i]], xpd=NA, cex=0.8)
    i <- i + 1
  }
  # Add dominant OTU assigned by 97% clustering within samples
  text(-0.03, -0.375, label=expression(paste(bold("C. "), "Dominant OTU assigned by 97% clustering within samples:")), xpd=NA, adj=c(0,0))
  names97bs <- unlist(lapply(strsplit(as.character(data.frame(tax_table(phy97bs)[apply(data.frame(otu_table(phy97bs), check.names=F)[,sord], MARGIN=2, FUN=which.max), "Subtype"])$Subtype), split="_"), "[[", 1))
  otus97bs <- rownames(tax_table(phy97bs)[apply(data.frame(otu_table(phy97bs), check.names=F)[,sord], MARGIN=2, FUN=which.max), ])
  newotus97bs <- c(1, which(diff(as.numeric(factor(otus97bs)))!=0)+1, nsamples(phy100)+1)
  # Get color of dominant 97% by-sample OTU representative sequence, match to color of representative unique sequence variant
  # Get representative sequence for dominant 97% by-sample OTU
  repseqs97bs <- sapply(otus97bs, function(dom97bsotu) {
    withinotu <- system(sprintf("awk '/%s / {print $2}' data/otus_97_bysample/all_rep_set_rep_set.fasta", dom97bsotu), intern=T)
    sample <- strsplit(withinotu, split="_")[[1]][1]
    repseqwithin <- system(sprintf("awk '/%s / {print $2}' data/otus_97_bysample/%s_rep_set.fasta", withinotu, sample), intern=T)
    return(repseqwithin)
  })
  # Find 100% OTU to which representative sequence for dominant 97% by-sample OTU belongs
  equivotu100 <- sapply(repseqs97bs, function(repseq) {
    system(sprintf("awk '/%s\t/ {print $1}' data/otus_100/nosingles_otus.txt", repseq), intern=T)
  })
  # Get matching color for dominant dominant 97% by-sample OTU
  matchcols <- topn[equivotu100, "colors"]
  # Plot rectangles with appropriate fill colors for 97% by-sample OTUs
  i <- 1
  while (i < length(newotus97bs)) {
    polygon(x=c(rep(starts[newotus97bs[i]], 2), rep(ends[newotus97bs[i+1]-1], 2)), y=c(-0.4, -0.5, -0.5, -0.4), 
            col=unique(matchcols)[i], xpd=NA)
    text((starts[newotus97bs[i]]+ends[newotus97bs[i+1]-1])/2, -0.45, labels=names97bs[newotus97bs[i]], xpd=NA, cex=0.8)
    i <- i + 1
  }
}
