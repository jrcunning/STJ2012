plotnet <- function(net) {
  plot(net,
       edge.curved=0.1, 
       edge.width=15*(E(net)$weight)^0.35,
       vertex.label=unlist(lapply(strsplit(V(net)$Subtype, split="_"), "[", 1)),
       vertex.color=ifelse(is.na(V(net)$Clade), "white",
                           taxcolors[factor(V(net)$Clade, levels=c("A","B","C","D","F","G"))]),
       vertex.size=ifelse(is.na(V(net)$Clade), 15, 10*sqrt(degree(net))))
}
