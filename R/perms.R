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