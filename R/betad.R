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