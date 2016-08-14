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