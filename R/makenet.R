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