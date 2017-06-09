# Set options
knitr::opts_knit$set(root.dir = normalizePath(".."))
knitr::opts_chunk$set(echo = F, warning = F, message = F, cache = T)
options(scipen = 3, digits = 9)

# Load functions
source("R/functions.R", .GlobalEnv)

# Load package libraries
library(phyloseq); library(rgdal)
library(vegan); library(multcompView)
library(reshape2); library(igraph)
library(stringr); library(RColorBrewer)
library(pBrackets)

# Set colors for plotting clades
taxcolors <- matrix(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462"), 
                    dimnames=list(c("CladeA", "CladeB", "CladeC", "CladeD", "CladeF", "CladeG")))

# ==========================
# Import data
# ==========================
# Import taxonomic assignment data
read.nw <- function(file) {
  nw <- read.table(file, stringsAsFactors=FALSE)
  nw <- nw[order(nw$otu),]
  return(nw)
}

tax97 <- read.nw("data/otus_97/nw_tophits.tsv")
tax100 <- read.nw("data/otus_100/nw_tophits.tsv")
tax97bs <- read.nw("data/otus_97_bysample/nw_tophits.tsv")

# Deal with identical taxonomic assignments (because some reference sequences are not unique...)
# Create names for groups of identical sequences based on members of group...
ident <- readLines("data/ITS2db_trimmed_notuniques_otus.txt")
ident <- gsub("denovo[0-9]*\t", "", ident)
ident <- strsplit(ident, split="\t")
ident2 <- lapply(ident, str_match_all, pattern="[A-I]{1}[0-9]{1,3}.*[_/]")
ident2 <- lapply(ident2, function(g) gsub(pattern="\\..*_$", "", x=unlist(g)))
ident2 <- lapply(ident2, function(g) gsub(pattern="_$", "", g))
ident2 <- lapply(ident2, function(g) unlist(strsplit(g, split="/")))
subtypes <- lapply(ident2, function(x) levels(factor(unlist(x))))
subtypes <- lapply(subtypes, function(s) paste(paste(s, collapse="/"), "_multiple", sep=""))
names(ident) <- subtypes
ident <- melt(do.call(rbind, ident))
ident <- unique(ident[order(ident[,1]), c(3,1)])

# Replace any sequence name in taxonomy assignment that is a member of a group of identical sequences with the name of the group
collapse.idents <- function(df) {
  within(df, {
    for (i in 1:nrow(ident)) {
      hit <- gsub(ident[i,1], ident[i,2], hit)
    }
  })
}

tax97 <- collapse.idents(tax97)
tax97bs <- collapse.idents(tax97bs)
tax100 <- collapse.idents(tax100)

tax97 <- with(tax97, tax97[order(otu, -sim), ])
tax97 <- tax97[!duplicated(tax97$otu), ]
rownames(tax97) <- tax97$otu

tax97bs <- with(tax97bs, tax97bs[order(otu, -sim), ])
tax97bs <- tax97bs[!duplicated(tax97bs$otu), ]
rownames(tax97bs) <- tax97bs$otu

tax100 <- with(tax100, tax100[order(otu, -sim), ])
tax100 <- tax100[!duplicated(tax100$otu), ]
rownames(tax100) <- tax100$otu

# Import OTU tables
otu97 <- otu_table(read.table("data/otus_97/97_otus.tsv", header=T, check.names=F, row.names=1,
                              skip=1, sep="\t", comment.char=""), taxa_are_rows=T)
otu97bs <- otu_table(read.table("data/otus_97_bysample/97_otus_bysample.tsv", header=T, check.names=F, row.names=1,
                              skip=1, sep="\t", comment.char=""), taxa_are_rows=T)
otu100 <- otu_table(read.table("data/otus_100/100_otus.tsv", header=T, check.names=F, row.names=1,
                               skip=1, sep="\t", comment.char=""), taxa_are_rows=T)

# Import sample data
sam <- sample_data(read.delim("data/mapping_file.txt", sep="\t", header=T, row.names=1))

# Build phyloseq objects
phy97 <- phyloseq(otu97, tax_table(as.matrix(tax97)), sam)
phy97bs <- phyloseq(otu97bs, tax_table(as.matrix(tax97bs)), sam)
phy100 <- phyloseq(otu100, tax_table(as.matrix(tax100)), sam)

# Filter out OTUs that are not Symbiodinium
## If the top hit in NCBI does not contain the string "Symbiodinium", then this sequence is assumed to not be Symbiodinium.
poormatch97 <- readLines("data/otus_97/poormatch_IDs.txt")
poormatch97 <- data.frame(otu=str_extract(poormatch97, "denovo[^ ]*"),
                          symbio=str_detect(poormatch97, "Symbiodinium"),   # TRUE if Symbiodinium
                          stringsAsFactors = FALSE)
notsymbio97 <- poormatch97[which(poormatch97$symbio==F), 1]

poormatch97bs <- readLines("data/otus_97_bysample/poormatch_IDs.txt")
poormatch97bs <- data.frame(otu=str_extract(poormatch97bs, "denovo[^ ]*"),
                            symbio=str_detect(poormatch97bs, "Symbiodinium"),   # TRUE if Symbiodinium
                            stringsAsFactors = FALSE)
notsymbio97bs <- poormatch97bs[which(poormatch97bs$symbio==F), 1]

poormatch100 <- readLines("data/otus_100/poormatch_IDs.txt")
poormatch100 <- data.frame(otu=str_extract(poormatch100, "denovo[^ ]*"),
                           symbio=str_detect(poormatch100, "Symbiodinium"),   # TRUE if Symbiodinium
                           stringsAsFactors = FALSE)
notsymbio100 <- poormatch100[which(poormatch100$symbio==F), 1]

# Remove OTUs that are not Symbiodinium
phy97.f <- prune_taxa(!taxa_names(phy97) %in% notsymbio97, phy97)
phy97bs.f <- prune_taxa(!taxa_names(phy97bs) %in% notsymbio97bs, phy97bs)
phy100.f <- prune_taxa(!taxa_names(phy100) %in% notsymbio100, phy100)

# Filter OTUs by minimum count
# Set threshold count
n <- 10
# Identify OTUs below threshold count
taxa97 <- taxa_sums(phy97.f)[which(taxa_sums(phy97.f) >= n)]
taxa97bs <- taxa_sums(phy97bs.f)[which(taxa_sums(phy97bs.f) >= n)]
taxa100 <- taxa_sums(phy100.f)[which(taxa_sums(phy100.f) >= n)]
# Remove taxa below threshold count
phy97.f <- prune_taxa(names(taxa97), phy97.f)
phy97bs.f <- prune_taxa(names(taxa97bs), phy97bs.f)
phy100.f <- prune_taxa(names(taxa100), phy100.f)

# Filter samples by minimum count
# Set threshold number of reads
sn <- 200
# Remove samples with fewer reads than threshold
phy97.f <- prune_samples(sample_sums(phy97.f)>=sn, phy97.f)
phy97bs.f <- prune_samples(sample_sums(phy97bs.f)>=sn, phy97bs.f)
phy100.f <- prune_samples(sample_sums(phy100.f)>=sn, phy100.f)

# Filter OTUs by minimum count again in case any dropped below threshold after filtering samples
# Identify OTUs below threshold count
taxa97 <- taxa_sums(phy97.f)[which(taxa_sums(phy97.f) >= n)]
taxa97bs <- taxa_sums(phy97bs.f)[which(taxa_sums(phy97bs.f) >= n)]
taxa100 <- taxa_sums(phy100.f)[which(taxa_sums(phy100.f) >= n)]
# Remove taxa below threshold count
phy97.f <- prune_taxa(names(taxa97), phy97.f)
phy97bs.f <- prune_taxa(names(taxa97bs), phy97bs.f)
phy100.f <- prune_taxa(names(taxa100), phy100.f)

# Label clades and subtypes for filtered phyloseq object tax_tables
get.st <- function(df) {
  within(df, {
    Clade <- substr(hit, 1, 1)
    Subtype <- gsub(hit, pattern="_[A-Z]{2}[0-9]{6}", replacement="")
    Subtype <- gsub(Subtype, pattern="_multiple", replacement="")
    Subtype2 <- ifelse(as.numeric(sim)==100, paste0("'", Subtype, "'"),
                       paste0("'", rep(rle(sort(Subtype))$values, times=rle(sort(Subtype))$lengths), "'^", 
                              unlist(lapply(rle(sort(Subtype))$lengths, seq_len)))[order(order(Subtype))])
    #Subtype <- ifelse(as.numeric(sim)==100, Subtype, paste("*", Subtype, sep=""))
  })
}

tax_table(phy97.f) <- as.matrix(get.st(data.frame(tax_table(phy97.f), stringsAsFactors=FALSE)))
tax_table(phy97bs.f) <- as.matrix(get.st(data.frame(tax_table(phy97bs.f), stringsAsFactors=FALSE)))
tax_table(phy100.f) <- as.matrix(get.st(data.frame(tax_table(phy100.f), stringsAsFactors=FALSE)))

# Transform count data
# Convert to proportion (relative abundance)
phy97.f.p <- transform_sample_counts(phy97.f, function(x) x/sum(x))
phy97bs.f.p <- transform_sample_counts(phy97bs.f, function(x) x/sum(x))
phy100.f.p <- transform_sample_counts(phy100.f, function(x) x/sum(x))
# Apply transformation function
transform <- function(x) sqrt(x/sum(x))  # Set transformation function
phy97.f.t <- transform_sample_counts(phy97.f, transform)  # Transform data
phy97bs.f.t <- transform_sample_counts(phy97bs.f, transform) 
phy100.f.t <- transform_sample_counts(phy100.f, transform)


# ==========================
# Figure 1
# ==========================
stjmap <- readOGR("data/coastline/zipfolder", layer="Coastal_Coastline_line", verbose=F)  # units are in meters
stjmap <- spTransform(stjmap, CRS("+proj=longlat +datum=WGS84"))  # transform to lat/long

sitecoords <- data.frame(
  name=c("Booby Rock", "Tektite", "Kiddel", "Groot Pan", "Yawzi Point", "Anna's Point", "Haulover Bay"),
  lat=c(18.30327, 18.31057, 18.30770, 18.31067, 18.31515, 18.36710, 18.35008),
  lon=c(-64.71057, -64.72203, -64.71398, -64.71813, -64.72513, -64.73337, -64.67933)
)

viers <- c(-64.722697, 18.322277)

png(file="figures/Fig1.png", height=3, width=6.85, units="in", res=300)
layout(mat=matrix(c(1,1,6,2,1,1,5,3,1,1,5,4), nrow=3, byrow=T))

# PLOT MAP
par(mar=c(0,0,0,0), oma=c(0,0,0,0), cex=1, xpd=NA)
plot(stjmap, lwd=1.5, col="gray40",
     xlim=c(-64.77, -64.70), ylim=c(18.303, 18.370))
#text(-64.735, 18.34, expression(italic("       St. John\nU.S. Virgin Islands")), col="gray40")
axis(side=3, at=seq(-64.78, -64.66, 0.01), line=0, tck=0.02, labels=FALSE)
axis(side=3, at=seq(-64.78, -64.66, 0.04), line=-2, lwd=0, cex.axis=0.6)
axis(side=2, at=seq(18.3, 18.4, 0.01), line=0, tck=0.02, labels=FALSE)
axis(side=2, at=seq(18.3, 18.4, 0.02), line=-2, lwd=0, cex.axis=0.6)
with(sitecoords, {
  points(lon, lat, pch=21, bg="black", cex=0.8)
  text(lon, lat+c(-0.001,0,-0.003,0.002,0.003,0,0), name, pos=c(4,2,2,4,4,1,3), cex=0.8, offset=0.3)
})
NS <- c(-64.711967, 18.368141)
SS <- c(-64.726840, 18.315465)
points(NS[1], NS[2], pch="N", cex=0.9, bg="gray80", col="black", lwd=1)
points(SS[1], SS[2], pch="S", cex=0.9, bg="gray40", col="black", lwd=1)
points(NS[1], NS[2], pch=1, cex=2.2, bg="gray80", col="black", lwd=1)
points(SS[1], SS[2], pch=1, cex=2.2, bg="gray40", col="black", lwd=1)

# Plot histogram of sample collections
par(new=T)
samdat <- data.frame(sample_data(phy97))
nsam <- with(samdat, aggregate(Species, by=list(Species=Species, Shore=Shore), FUN=length))
nsam <- dcast(nsam, Species~Shore)
rownames(nsam) <- nsam$Species
rownames(nsam) <- str_replace_all(rownames(nsam), 
                c("alcicornis"="M.alcicornis", "annularis"="O.annularis", "astreoides"="P.astreoides","cavernosa"="M.cavernosa",
                  "cylindrus"="D.cylindrus", "fragum"="F.fragum", "furcata"="P.furcata", "radians"="S.radians",
                  "siderea"="S.siderea", "strigosa"="P.strigosa"))
nsam <- nsam[,2:3]
par(mar=c(7,5,5,6), mgp=c(0.5,0.2,0))
b <- barplot(t(as.matrix(nsam)), las=2, ylab="# samples", xpd=NA, names.arg=rep("",10), cex.axis=0.5, cex.names=0.5, cex.lab=0.5, tcl=0.1)
text(b+1.2, -0.5, labels=rownames(nsam), srt=45, pos=2, cex=0.5, adj=1)
legend("topleft", pch=22, pt.bg=c("gray90", "gray40"), pt.cex=1, cex=0.5, legend=c("N","S"), bty="n", inset=c(0.12,-0.1), xpd=T, y.intersp=0.8, x.intersp=0.75)
text(par("usr")[1], par("usr")[4], "A)", cex=0.75, adj=2.5)

# Plot environmental data
envdata <- read.csv("data/STJ_envdata2.csv")
envdata$Site <- factor(envdata$Site, levels=c("S", "N"))

#SST
par(cex=0.7, mar=c(3,1,1,1), mgp=c(1,0.25,0), xpd=NA, cex.axis=0.7, cex.lab=0.7, tcl=-0.25)
plot(NA, xlim=c(25,30), ylim=c(0.5,2.5), bty="n", yaxt="n", ylab="", xlab="SST (°C)", xaxt="n")
axis(side=1, at=seq(25,30,1))
text(par("usr")[1], c(1,2), c("S","N"), pos=2, offset=0.5)
points(rep(par("usr")[1]-strwidth("cc"),2), c(1,2), pch=1, cex=2.5)
boxplot(envdata$SST ~ envdata$Site, range=0, add=T, axes=F, horizontal=T)
text(par("usr")[1], par("usr")[4], "B)", adj=3)

#chlA
plot(NA, xlim=c(0,0.1), ylim=c(0.5,2.5), bty="n", yaxt="n", ylab="", xlab="ChlA (mg m-3)", xaxt="n")
axis(side=1, at=seq(0,0.1,0.01))
text(par("usr")[1], c(1,2), c("S","N"), pos=2, offset=0.5)
points(rep(par("usr")[1]-strwidth("cc"),2), c(1,2), pch=1, cex=2.5)
boxplot(envdata$chlA ~ envdata$Site, range=0, add=T, axes=F, horizontal=T)
text(par("usr")[1], par("usr")[4], "C)", adj=3)

#Wave exposure
we <- list(N=3.234787362,	S=4.025022708)  # Wave exposure data provided by Iliana Chollett
barplot(c(we$S, we$N), horiz=T, xlim=c(0,5), ylab="", xlab="Wave exposure (ln J m-3)")
text(par("usr")[1], c(0.7,1.9), c("S","N"), pos=2, offset=0.5, xpd=NA)
points(rep(par("usr")[1]-strwidth("cc"),2), c(0.7,1.9), pch=1, cex=2.5, xpd=NA)
text(par("usr")[1], par("usr")[4], "D)", adj=3)

dev.off()

# ==========================
# Table 1
# ==========================
# Compute summary statistics 
stats97.f <- data.frame(`97% OTUs<br>(across samples)`=t(phystats(phy97.f)), check.names=F)
stats97bs.f <- data.frame(`97% OTUs<br>(within samples)`=t(phystats(phy97bs.f)), check.names=F)
stats100.f <- data.frame(`100% OTUs`=t(phystats(phy100.f)), check.names=F)

knitr::kable(cbind(stats97.f, stats97bs.f, stats100.f)[c(2, 3, 6, 8), ],
             caption = "Summary statistics for each clustering approach.")

# ==========================
# Figure 2
# ==========================
png("figures/Fig2.png", width=6, height=4, units="in", res=300)
set.seed(42432) 
par(mfrow=c(1,1))
otubarplot2(samples=c("64","72","29","45","73","59","47","26","48","66","63"))
dev.off()

# ==========================
# Symbiodinium community composition
# ==========================
## Calculate relative abundance of each clade
cladeAbund <- aggregate(data.frame(RelAbund=rowSums(otu_table(phy97bs.f.p))),
                        by=list(Clade=data.frame(tax_table(phy97bs.f.p))$Clade), FUN=sum)
cladeAbund$Prop <- round(prop.table(cladeAbund$RelAbund) * 100, 1)
cladeAbund$Notus <- table(data.frame(tax_table(phy97bs.f.p))$Clade)
clademax <- aggregate(data.frame(max=apply(otu_table(phy97bs.f.p), 1, FUN=max)),
                      by=list(Clade=data.frame(tax_table(phy97bs.f.p))$Clade), FUN=max)

# ==========================
# Figure 3
# ==========================
png(file="figures/Fig3.png", width = 4, height=6, units="in", res=300)
par(mfrow=c(1,1), mar=c(2, 1.5, 2, 5), lwd=0.1, cex=1)
# Plot composition of 97% within-sample OTUs colored by clade
composition(phy97bs.f.p, col=taxcolors[factor(data.frame(tax_table(phy97bs.f.p))[order(data.frame(tax_table(phy97bs.f.p))$Subtype), ]$Clade, levels=c("A","B","C","D","F","G"))], legend=T)
dev.off()

# ==========================
# Table 2
# ==========================
set.seed(43789)  
shorestats97bs <- perms(phy97bs.f.p, group="Species", trt="Shore")
shorestats97bs$Species <- str_replace_all(
  shorestats97bs$Species,
  c("alcicornis"="*Millepora alcicornis*",
    "annularis"="*Orbicella annularis*",
    "astreoides"="*Porites astreoides*",
    "cavernosa"="*Montastraea cavernosa*",
    "cylindrus"="*Dendrogyra cylindrus*",
    "fragum"="*Favia fragum*",
    "furcata"="*Porites furcata*",
    "radians"="*Siderastrea radians*",
    "siderea"="*Siderastrea siderea*",
    "strigosa"="*Pseudodiplora strigosa*")
)
knitr::kable(shorestats97bs, digits=3, row.names=F, 
             caption="Mean overall, within-shore ('within'), and between-shore ('between') Bray-Curtis dissimilarities of the *Symbiodinium* communities in each host species, and PERMANOVA tests (partial R^2^ and p-values) for a difference between shores. Between-shore tests were not possible for *O. annularis* and *P. astreoides* since they were only sampled from one shore.")


# ==========================
# Network analysis of *Symbiodinium* metacommunity
# ==========================
## Create dominant symbionts network
set.seed(12374)
domnet <- sppnet(phy=phy97bs.f.p, plot=F,
                 fun=function(x) length(x[x>0.5])/length(x), 
                 layout=layout.fruchterman.reingold)
# Calculate number of dominant symbionts connected to each species
dom.Nedges <- sapply(levels(data.frame(sample_data(phy97bs.f.p))$Species),
                     FUN=function(x) length(E(domnet)[from(x)]))
dom.Nedges <- sort(dom.Nedges, decreasing=T)

# Create abundant symbionts network
set.seed(12374)
abunet <- sppnet(phy=phy97bs.f.p, plot=F,
                 fun=function(x) length(x[x>0.01])/length(x),
                 layout=layout.fruchterman.reingold)
# Calculate number of abundant symbionts connected to each species
abu.Nedges <- sapply(levels(data.frame(sample_data(phy97bs.f.p))$Species),
                     FUN=function(x) length(E(abunet)[from(x)]))
abu.Nedges <- sort(abu.Nedges, decreasing=T)

# Create background symbionts network
# filter out OTUs that only occurred in one sample
phy97bs.f.p2 <- filter_taxa(phy97bs.f.p, function(x) sum(x>0) > 1, prune=TRUE)
bkgnet <- sppnet(phy=phy97bs.f.p2, plot=F,
                 fun=function(x) any(x>0 & x<0.01),
                 layout=layout.fruchterman.reingold)
# Merge D nodes (assuming they are all D1a/trenchii)
# Get original node map (sequence starting at one to the total number of nodes)
nodemap <- seq(1:length(V(bkgnet)))
# Make a copy of the nodemap that will be modified
nodemapnew <- nodemap
# Find all nodes beginning with D, and give them the same node ID number
nodemapnew[grep("^D", V(bkgnet)$Subtype)] <- min(grep("^D", V(bkgnet)$Subtype))
# Replace the now missing values in the node map so node IDs remain sequential
while (max(setdiff(nodemap, nodemapnew)) < max(nodemap)) {
  nodemapnew[which(nodemapnew==max(setdiff(nodemap, nodemapnew))+1)] <- min(setdiff(nodemap, nodemapnew))
}
# Contract the network using the new nodemap (nodes with same node ID are merged)
bkgnet2 <- contract(bkgnet, nodemapnew, vertex.attr.comb=list(size="sum", "last"))
# Simplify network so that multiple edges are combined
bkgnet2 <- simplify(bkgnet2, remove.multiple=T, edge.attr.comb=list(weight="mean", width="mean", "ignore"))
# remove symbionts that are only connected to one or two coral species
bkgnet3 <- delete.vertices(bkgnet2, degree(bkgnet2)<=2) # & V(net2)$type==2)
V(bkgnet3)$size <- ifelse(is.na(V(bkgnet3)$Clade), 15, 1.5*degree(bkgnet3)^1.5)
E(bkgnet3)$width <- E(bkgnet3)$weight * 2

# Calculate number of edges for background symbionts
nspp <- sapply(V(bkgnet3)[V(bkgnet3)$type==2], function(x) length(E(bkgnet3)[from(x)]))
names(nspp) <- V(bkgnet3)$label[match(names(nspp), V(bkgnet3)$name)]
nspp <- sort(nspp, decreasing=T)

# ==========================
# Figure 4
# ==========================
# Plot abundant symbionts network
png("figures/Fig4.png", width=5, height=5, units="in", res=300)
par(mfrow=c(1,1), mar=c(0,0,0,0), lwd=1)
V(abunet)$size <- ifelse(is.na(V(abunet)$Clade), 10, 13*(degree(abunet))^0.4)
set.seed(12374)
l <- layout.fruchterman.reingold(abunet)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(abunet, rescale=F, layout=l, edge.curved=0.1, 
     vertex.label.cex=ifelse(is.na(V(abunet)$Clade), V(abunet)$size/15, V(abunet)$size/18),
     vertex.label.color="black") 
dev.off()

# ==========================
# Figure 5
# ==========================
png("figures/Fig5.png", width=7.20472, height=3.60236, units="in", res=300)
par(mfrow=c(1,2), mar=c(0.1,1,1.9,1), lwd=1, xpd=NA) 
# Plot dominant symbionts network
set.seed(12374)
l <- layout.fruchterman.reingold(domnet)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
V(domnet)$size <- ifelse(is.na(V(domnet)$Clade), 16, 18*(degree(domnet))^0.5)
plot(domnet, rescale=F, layout=l, edge.curved=0.1, vertex.label.cex=V(domnet)$size/25,
     vertex.label.color="black")
title("A.)", adj=0)

# Plot background symbionts network
set.seed(42384)
l3 <- layout.fruchterman.reingold(bkgnet3)
l3 <- norm_coords(l3, ymin=-1, ymax=1, xmin=-1, xmax=1)
V(bkgnet3)$size <- ifelse(is.na(V(bkgnet3)$Clade), 16, 8*(degree(bkgnet3))^0.8)
plot(bkgnet3, rescale=F, layout=l3, edge.curved=0.1, vertex.label.cex=V(bkgnet3)$size/25,
     vertex.label.color="black")
points(0.9,0.9, pch=0, cex=3.5)
text(0.9,0.9, "Oa", cex=0.7)
title("B.)", adj=0)
dev.off()

# ==========================
# Figure 6
# ==========================
betad97bs <- betad(phy97bs.f.t, group="Species")

png("figures/Fig6.png", width=3.5, height=3.5, units="in", res=300)
par(mar=c(7,3,1,1), mgp=c(1.8,0.1,0), tcl=0.25)
with(betad97bs$sambdsumm.ord, {
  plot(mean, type="n", ylim=c(0, 0.65), ylab="Distance to centroid", xaxt="n", xlab="")
  arrows(1:10, mean - se, 1:10, mean + se, length=0.05, angle=90, code=3)
  points(1:10, mean, cex=1, pch=21, bg="white")
  axis(side=1, at=1:10, cex.axis=0.75, labels=NA)
  text(1.4:10.4, -0.04, srt=45, xpd=T, pos=2, font=3,
       labels=str_replace_all(levels(sam$Species)[order(betad97bs$sambdsumm$mean, decreasing=T)],
                              c("alcicornis"="M. alcicornis",
                                "annularis"="O. annularis",
                                "astreoides"="P. astreoides",
                                "cavernosa"="M. cavernosa",
                                "cylindrus"="D. cylindrus",
                                "fragum"="F. fragum",
                                "furcata"="P. furcata",
                                "radians"="S. radians",
                                "siderea"="S. siderea",
                                "strigosa"="P. strigosa")))
  text(1:10, mean + se + 0.03, labels=betad97bs$saml$Letters[as.character(group)], cex=0.5)
})
dev.off()
