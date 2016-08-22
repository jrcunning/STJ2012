
phy100 <- subset_samples(phy100.f.p, rownames(sample_data(phy100.f.p)) %in% c("64", "72", "29",
                                                                              "45", "73",
                                                                              "59",
                                                                              "47", "26",
                                                                              "48", "66",
                                                                              "63"))
phy100 <- filter_taxa(phy100, function(x) sum(x)!=0, TRUE)
phy97 <- subset_samples(phy97.f.p, rownames(sample_data(phy97.f.p)) %in% c("64", "72", "29",
                                                                           "45", "73",
                                                                           "59",
                                                                           "47", "26",
                                                                           "48", "66",
                                                                           "63"))
phy97 <- filter_taxa(phy97, function(x) sum(x)!=0, TRUE)
phy97bs <- subset_samples(phy97bs.f.p, rownames(sample_data(phy97bs.f.p)) %in% c("64", "72", "29",
                                                                           "45", "73",
                                                                           "59",
                                                                           "47", "26",
                                                                           "48", "66",
                                                                           "63"))
phy97bs <- filter_taxa(phy97bs, function(x) sum(x)!=0, TRUE)



otubarplot2 <- function(phy100, phy97, phy97bs) {
  par(mgp=c(0,0.5,0), cex.axis=0.75, mar=c(6,2,3,1))
  samdat <- data.frame(sample_data(phy100))
  typerelabund <- as.matrix(otu_table(phy100)[order(data.frame(tax_table(phy100))$Subtype), 
                                           rownames(samdat)])
  typerelabund <- typerelabund[, rownames(sample_data(phy100)[with(sample_data(phy100), 
                                                                order(Shore, Site, InputFileName))])]
  # Get info for plotting OTU names and blast hits on top of barplot
  blasthits <- as.character(data.frame(tax_table(phy100))[order(data.frame(tax_table(phy100))$Subtype), "Subtype"])
  names <- as.character(rownames(data.frame(tax_table(phy100)))[order(data.frame(tax_table(phy100))$Subtype)])
  samples <- colnames(typerelabund)
  shores <- data.frame(sample_data(phy100))[samples, "Shore"]
  sites <- data.frame(sample_data(phy100))[samples, "Site"]
  divides <- apply(typerelabund, 2, cumsum)
  heights <- diff(divides)
  heights[heights < 0.04] <- NA
  # Plot barplot and OTU names and blast hits
  ranks <- order(rowSums(typerelabund), decreasing=T)
  topten <- rownames(typerelabund[ranks[1:45], ])
  topten <- data.frame(tax_table(phy100)[topten, "Clade"])
  nclade <- table(data.frame(topten)$Clade)
  colors <- rev(colorRampPalette(brewer.pal(9,"Blues")[-1])(nclade[1]+1))
  colors <- rainbow(nclade[1]+1)
  color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  colors <- sample(color, 45)
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',][-c(3:5), ]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colors <- sample(col_vector, 45)
  
  
  topten$colors <- colors
  topten
  rownames(topten)
  sord <- rownames(sample_data(phy97)[order(factor(tax_table(phy97)[apply(otu_table(phy97), MARGIN=2, FUN=which.max), "Subtype"]), sample_data(phy97)$Species), ])
  #sord <- rownames(sample_data(phy97))[order(sample_data(phy97)$Species)]
  bars <- barplot(typerelabund[order(typerelabund[, sord[1]], decreasing=TRUE), sord[1]], width=1/(ncol(typerelabund)*1), xlim=c(0, 1),
                  col=topten[rownames(typerelabund)[order(typerelabund[, sord[1]], decreasing=TRUE)], "colors"], las=1, space=0, xaxs="i",
                  xlab="", ylab="", xaxt="n")
  for (i in 2:ncol(typerelabund)) {
    barplot(typerelabund[order(typerelabund[, sord[i]], decreasing=TRUE), sord[i]], add=T, yaxt="n", xaxt="n", width=1/(ncol(typerelabund)*1), space=1*(i-1),
            col=topten[rownames(typerelabund)[order(typerelabund[, sord[i]], decreasing=TRUE)], "colors"])
  }
  midpts <- seq(from=bars, by=bars*2, length.out=nsamples(phy100))
  starts <- midpts - bars
  ends <- midpts + bars
  text(0, 1.20, label="Composition of unique sequence variants:", xpd=NA, adj=c(0,0))
  # Plot species names
  spnames <- as.character(sample_data(phy97)$Species[match(sord, rownames(sample_data(phy97)))])
  newsp <- c(1, which(diff(as.numeric(factor(spnames)))!=0)+1, nsamples(phy100)+1)
  i <- 1
  while (i < length(newsp)) {
    brackets(starts[newsp[i]], 1, ends[newsp[i+1]-1], 1, xpd=NA)
    text((starts[newsp[i]]+ends[newsp[i+1]-1])/2, 1.1, labels=spnames[newsp[i]], xpd=NA, cex=0.8)
    i <- i + 1
  }
  
  text(0, -0.10, label="Dominant OTU assigned by 97% clustering across samples:", xpd=NA, adj=c(0,0))
  names97 <- unlist(lapply(strsplit(as.character(data.frame(tax_table(phy97)[apply(data.frame(otu_table(phy97), check.names=F)[,sord], MARGIN=2, FUN=which.max), "Subtype"])$Subtype), split="_"), "[[", 1))
  otus97 <- rownames(tax_table(phy97)[apply(data.frame(otu_table(phy97), check.names=F)[,sord], MARGIN=2, FUN=which.max), ])
  newotus97 <- c(1, which(diff(as.numeric(factor(otus97)))!=0)+1, nsamples(phy100)+1)
  i <- 1
  while (i < length(newotus97)) {
    rect(starts[newotus97[i]], -0.125, ends[newotus97[i+1]-1], -0.225, xpd=NA)
    text((starts[newotus97[i]]+ends[newotus97[i+1]-1])/2, -0.175, labels=names97[newotus97[i]], xpd=NA, cex=0.8)
    i <- i + 1
  }
  
  text(0, -0.35, label="Dominant OTU assigned by 97% clustering within samples:", xpd=NA, adj=c(0,0))
  names97bs <- unlist(lapply(strsplit(as.character(data.frame(tax_table(phy97bs)[apply(data.frame(otu_table(phy97bs), check.names=F)[,sord], MARGIN=2, FUN=which.max), "Subtype"])$Subtype), split="_"), "[[", 1))
  otus97bs <- rownames(tax_table(phy97bs)[apply(data.frame(otu_table(phy97bs), check.names=F)[,sord], MARGIN=2, FUN=which.max), ])
  newotus97bs <- c(1, which(diff(as.numeric(factor(otus97bs)))!=0)+1, nsamples(phy100)+1)
  i <- 1
  while (i < length(newotus97bs)) {
    rect(starts[newotus97bs[i]], -0.375, ends[newotus97bs[i+1]-1], -0.475, xpd=NA)
    text((starts[newotus97bs[i]]+ends[newotus97bs[i+1]-1])/2, -0.425, labels=names97bs[newotus97bs[i]], adj=0.5, xpd=NA, cex=0.8)
    i <- i + 1
  }
}


set.seed(42432)
otubarplot2(phy100, phy97, phy97bs)
