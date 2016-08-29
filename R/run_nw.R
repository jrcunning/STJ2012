# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
# If two arguments not provided, return an error
if (length(args) < 2) {
  stop("must specify query sequences (.fasta) and reference sequences (.fasta)", call.=FALSE)
}

# Import query sequences
library(Biostrings)
otus <- readDNAStringSet(args[1])
names(otus) <- gsub(" .*$", "", names(otus))

# Import reference sequences
ref <- readDNAStringSet(args[2])

# import substitution matrix
EDNAFULL <- as.matrix(read.table("data/EDNAFULL", skip=8, header=T, row.names=1))

# Build function to find which reference sequences produces best alignment score using
# Needleman-Wunsch algorithm with default EMBOSS parameters for DNA sequences.
# -type="global-local" forces whole reference (pattern) sequence to be used, but subject/query may be 
# a concescutive subsequence. This is because queries may include slightly longer region of gene 
# than references because references were trimmed by o-smart-trim. This means that query sequence 
# can have more bases than the reference at beginning or end of the sequence (=end gaps in 
# alignment) without penalties.
get.best <- function(x) {
  scores <- pairwiseAlignment(subject=x, pattern=ref, type="global-local", scoreOnly=TRUE,
                           substitutionMatrix=EDNAFULL,
                           gapOpening=10.0, gapExtension=0.5)
  best <- paste(names(ref)[which(scores==max(scores))], collapse=";")
  return(best)
}

# Set up for parallel computation
library(parallel)
cl <- makeCluster(detectCores() - 1)  # Initiate cluster
clusterEvalQ(cl, library(Biostrings))  # Make Biostrings library available to cluster
clusterExport(cl, c("EDNAFULL", "ref", "otus", "get.best"))  # Export objects to cluster

# Run function in parallel to get best hit for each sequence
best <- parSapply(cl, otus, FUN=get.best)
stopCluster(cl)  # Stop cluster

# Collect results in a data frame
results <- data.frame(otu=names(otus), nw.hits=best, stringsAsFactors=FALSE)

# Get percent similarity of global alignment of each sequence with its best hit from database
# If multiple best hits, select first one for comparison
results$nw.hit1 <- unlist(lapply(strsplit(results$nw.hits, split=";"), "[", 1))
# Align each otu sequence to its best hit
psa <- pairwiseAlignment(subject=otus[results$otu], pattern=ref[results$nw.hit1], type="global-local",
                         substitutionMatrix=EDNAFULL,
                         gapOpening=10.0, gapExtension=0.5)

# Calculate percent similarity to reference with indels counting as single difference regardless of length
results <- within(results, {
  nw.len <- width(pattern(psa))
  nw.mis <- nmismatch(psa)
  nw.ins <- unlist(lapply(as.list(insertion(psa)), length))
  nw.del <- unlist(lapply(as.list(deletion(psa)), length))
  nw.sim <- round((nw.len-nw.mis-nw.ins-nw.del)/nw.len, 4) * 100
})

# Write results to .tsv file
write.table(results, file=paste(dirname(args[1]), "/", "nw_tophits.tsv", sep=""), 
            sep="\t", quote=FALSE)

