betad <- function(phy, group) {
  # Calculate distances to species centroids for each sample in principal coordinate space
  samdat <- data.frame(sample_data(phy))
  braydist <- vegdist(t(otu_table(phy)), method="bray")
  sambd <- betadisper(d=braydist, group=samdat[, group], 
                      type="centroid", bias.adjust=TRUE)
  # Calculate each species' mean distance to centroid and standard error
  sambdsumm <- aggregate(data.frame(mean=sambd$distances), by=list(group=sambd$group), FUN=mean)
  sambdsumm$se <- aggregate(sambd$distances, by=list(group=sambd$group), 
                            FUN=function(x) sd(x)/sqrt(length(x)))$x
  # Determine order of decreasing mean dispersion
  sambdsumm.ord <- sambdsumm[order(sambdsumm$mean, decreasing=T), ]
  # Update betadisper object with species in decreasing order of mean dispersion
  sambd <- betadisper(d=braydist, 
                      group=factor(samdat[, group], levels=as.character(sambdsumm.ord$group)), 
                      type="centroid", bias.adjust=TRUE)
  # Use permutations to perform pairwise comparisons of group mean dispersions
  sampt <- permutest(sambd, pairwise=T)
  saml <- multcompLetters(sampt$pairwise$permuted)
  # Return results as a list
  return(list(sambdsumm.ord=sambdsumm.ord, sambdsumm=sambdsumm, saml=saml))
}

composition <- function(phy, col, legend=T) {
  samdat <- data.frame(sample_data(phy))
  samdat$Genus <- factor(samdat$Genus, levels=rev(levels(samdat$Genus)))
  samdat$Species <- factor(samdat$Species, levels=rev(levels(samdat$Species)))
  samdat$Shore <- factor(samdat$Shore, levels=rev(levels(samdat$Shore)))
  samdat <- samdat[with(samdat, order(Genus, Species, Shore)), ]
  typerelabund <- as.matrix(otu_table(phy)[order(data.frame(tax_table(phy))$Subtype), 
                                           rownames(samdat)])
  shorebreaks <- c(as.character(samdat$Shore), "X")==c("X", as.character(samdat$Shore))
  shorebreaks <- which(shorebreaks==F) - 1
  spbreaks <- c(which(duplicated(samdat$Species)==F) - 1, nrow(samdat))
  # Make Barplot
  barplot(typerelabund, horiz=T, space=0, axes=F,axisnames=F, yaxs="i", col=col)
  rect(0, 0, par("usr")[2], par("usr")[4], lwd=1, xpd=T)
  axis(side=1, at=seq(0, 1, 0.1), line=0, tck=-0.025, mgp=c(0,0.25,0), cex.axis=0.7)
  mtext(side=1, "Relative abundance", cex=0.7, line=1)
  # Add legend
  if (legend==T) {
    legend(x=par("usr")[2]/2, y=par("usr")[4], xjust=0.5, yjust=0.25, horiz=T, bty="n", xpd=T, 
           cex=0.7, legend=c("A", "B", "C", "D", "F", "G"), fill=taxcolors, x.intersp=0.5)
    legend(x=par("usr")[2]-0.03, y=par("usr")[4], xjust=0, yjust=0.1, bty="n", xpd=T, cex=0.6, 
           pt.cex=0.6, legend=c("N Shore", "S Shore"), fill=c(0, "gray40"), y.intersp=0.7, 
           x.intersp=0.3)
  }
  # Add grouping bars for Shore
  for (i in 1:length(shorebreaks)) {
    lines(c(0, 1), c(shorebreaks[i], shorebreaks[i]), lty=2, lwd=0.25)
    rect(1.01, shorebreaks[i], 1.04, shorebreaks[i+1], col=rep(c("gray40", 0), 50)[i], 
         lwd=0.25, xpd=T)
  }
  # Add lines to separate species and species names
  for (i in 1:length(spbreaks)) {
    lines(c(0, 1.07), c(spbreaks[i], spbreaks[i]), xpd=T, type="l", lwd=0.4)
    text(1.03, (spbreaks[i] + spbreaks[i+1]) / 2, xpd=T, pos=4, cex=0.6,
         labels=paste(samdat$Genus[which(duplicated(samdat$Species)==F)][i], "\n",
                      samdat$Species[which(duplicated(samdat$Species)==F)][i], sep=""))
  }
}

makenet <- function(physeq, n) {
  require(phyloseq)
  require(reshape2)
  require(igraph)
  # Remove taxa that are not present in subset
  physeq <- prune_taxa(rowSums(otu_table(physeq))!=0, physeq)
  # Convert OTU table to "edges" table (weight = relative abundance)
  otudf <- data.frame(t(otu_table(physeq)), check.names=F)
  edges <- setNames(melt(as.matrix(otudf)), nm=c("source", "target", "weight"))  # Melt data
  edges <- edges[edges$weight>n, ]
  # Create "nodes" table from sample and tax data
  sdf <- data.frame(sample_data(physeq))
  sdf$id <- rownames(sdf)
  tdf <- data.frame(tax_table(physeq))
  tdf$id <- rownames(tdf)
  nodes <- merge(sdf, tdf, all=T)
  # Create network
  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
  return(net)
}

otubarplot <- function(phy, main) {
  samdat <- data.frame(sample_data(phy))
  typerelabund <- as.matrix(otu_table(phy)[order(data.frame(tax_table(phy))$Subtype), 
                                           rownames(samdat)])
  typerelabund <- typerelabund[, rownames(sample_data(phy)[with(sample_data(phy), 
                                                                order(Shore, Site, InputFileName))])]
  # Get info for plotting OTU names and blast hits on top of barplot
  blasthits <- as.character(data.frame(tax_table(phy))[order(data.frame(tax_table(phy))$Subtype), "Subtype"])
  names <- as.character(rownames(data.frame(tax_table(phy)))[order(data.frame(tax_table(phy))$Subtype)])
  samples <- colnames(typerelabund)
  shores <- data.frame(sample_data(phy))[samples, "Shore"]
  sites <- data.frame(sample_data(phy))[samples, "Site"]
  divides <- apply(typerelabund, 2, cumsum)
  heights <- diff(divides)
  heights[heights < 0.04] <- NA
  # Plot barplot and OTU names and blast hits
  bars <- barplot(typerelabund, col=rainbow(ntaxa(phy)), las=1, space=0.1, xaxs="i",
                  xlab="", ylab="", xaxt="n")
  text(bars, -0.06, labels=samples, xpd=T)
  text(bars, -0.10, labels=shores, xpd=T)
  text(bars, -0.14, labels=substr(sites, 1, 2), xpd=T)
  text(bars[1]/4, c(-0.06, -0.10, -0.14), labels=c("sample:", "shore:", "site:"), xpd=T, pos=2)
  mtext(side=2, text = "Relative Abundance", line=2)
  for (i in 1:length(bars)) {
    text(rep(bars[i], length(heights[which(!is.na(heights[,i])),i])), 
         divides[which(!is.na(heights[,i])),i] + heights[which(!is.na(heights[,i])),i] / 2, 
         labels=paste(names[which(!is.na(heights[,i]))+1],blasthits[which(!is.na(heights[,i]))+1],sep="\n"),
         cex=0.75)
  }
  mtext(side=3, main, xpd=T, adj=0, line=0.5, font=2)
}

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

# Function to test differences in community composition within and between groups by PERMANOVA

perms <- function(phy, group, trt) {
  groups <- data.frame(sample_data(phy))[, group]
  # Create data frame to populate with results
  resdf <- data.frame(
    matrix(ncol=8, nrow=10,
           dimnames=list(levels(data.frame(sample_data(phy))[, group]), 
                         c(group, "n", "overall", "within", "between", "R2", "bd", "p")))
  )
  # Compute within- and between-treatment statistics for each group
  for (i in 1:nlevels(groups)) {
    group.name <<- levels(groups)[i]
    resdf[i, group] <- group.name
    # make subsetted phyloseq object for individual group
    phy.group <- phy 
    samdat <- data.frame(sample_data(phy))
    sample_data(phy.group) <- sample_data(samdat[samdat[[group]]==group.name, ])
    samdat <- as(sample_data(phy.group), "data.frame") # sample data for subsetted phyloseq object
    # calculate distances
    braydist <- phyloseq::distance(phy.group, "bray")
    resdf$n[i] <- nsamples(phy.group)
    md <- meandist(braydist, samdat[, trt])
    resdf$within[i] <- summary(md)$W  # Weighted mean dissimilarity within treatments
    resdf$between[i] <- summary(md)$B  # Mean dissimilarity between treatments
    resdf$overall[i] <- summary(md)$D  # Overall dissimilarity
    if (nlevels(as(sample_data(phy.group), "data.frame")[, trt]) > 1) {
      permanova <- adonis(braydist ~ get(trt), data=samdat, permutations=99999)
      resdf$R2[i] <- permanova$aov.tab$"R2"[1]  # PERMANOVA partial R-squared
      resdf$p[i] <- permanova$aov.tab$"Pr(>F)"[1]  # PERMANOVA p-value
      bd <- betadisper(braydist, group=samdat[, trt])
      resdf$bd[i] <- TukeyHSD(bd)$group[4]
    } else {
      resdf$R2[i] <- NA
      resdf$p[i] <- NA
      resdf$bd[i] <- NA
    }
  }
  return(resdf)
}

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

plotnet <- function(net) {
  plot(net,
       edge.curved=0.1, 
       edge.width=15*(E(net)$weight)^0.35,
       vertex.label=parse(text=V(net)$Subtype2),
       vertex.color=ifelse(is.na(V(net)$Clade), "white",
                           taxcolors[factor(V(net)$Clade, levels=c("A","B","C","D","F","G"))]),
       vertex.size=ifelse(is.na(V(net)$Clade), 15, 10*sqrt(degree(net))))
}

sppnet <- function(phy, fun, layout, plot=T) {
  # Get otu table, sample data, and tax table from phyloseq object
  otudf <- data.frame(t(otu_table(phy)), check.names=F)
  samdf <- data.frame(sample_data(phy))
  taxdf <- data.frame(tax_table(phy))
  # Aggregate by species according to function provided (fun)
  sp.ag <- aggregate(otudf, by=list(Species=samdf$Species), FUN=fun)
  rownames(sp.ag) <- sp.ag$Species
  sp.ag <- sp.ag[, -1]
  edges <- setNames(melt(as.matrix(sp.ag)), nm=c("source", "target", "weight"))  # Melt data
  edges <- droplevels(edges[edges$weight>0, ]) # Remove zeros
  # Create "nodes" table from sample and tax data
  sdf <- data.frame(id=levels(samdf$Species))
  tdf <- taxdf[rownames(taxdf) %in% edges$target, ]
  tdf$id <- rownames(tdf)
  nodes <- merge(sdf, tdf, all=T)
  nodes$type <- ifelse(is.na(nodes$Clade), 1, 2)
  # Create network
  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
  # Modify node (vertex) attributes
  V(net)$shape <- c("square", "circle")[factor(V(net)$type)]
  V(net)$size <- ifelse(is.na(V(net)$Clade), 15, 10*(degree(net))^0.75) #5*sqrt(degree(net))
  V(net)$label <- parse(text=str_replace_all(
    ifelse(is.na(V(net)$Clade), names(V(net)), V(net)$Subtype2),
    c("alcicornis"="M.alc.", "annularis"="O.ann",
      "astreoides"="P.ast.", "cavernosa"="M.cav.",
      "cylindrus"="D.cyl.", "fragum"="F.fra.", 
      "furcata"="P.fur.", "radians"="S.rad.", 
      "siderea"="S.sid.", "strigosa"="D.str.")))
  V(net)$color <- ifelse(is.na(V(net)$Clade), "white", taxcolors[factor(V(net)$Clade, levels=c("A","B","C","D","F","G"))])
  E(net)$color <- "gray60"
  E(net)$width <- 15 * E(net)$weight
  # Plot network
  if (plot==T) {
    par(mar=c(0,0,0,0))
    l <- layout(net)
    l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    plot(net, rescale=F, layout=l, edge.curved=0.1, vertex.label.cex=V(net)$size/20,
         vertex.label.color="black")
  }
  return(net)
}

