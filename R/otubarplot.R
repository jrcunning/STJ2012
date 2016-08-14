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