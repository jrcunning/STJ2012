# get command line arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("must specify query sequences and reference sequences", call.=FALSE)
}

# Import query sequences
library(Biostrings)
otus <- readDNAStringSet(args[1])
names(otus) <- gsub(" .*$", "", names(otus))

# Import reference sequences
ref <- readDNAStringSet(args[2])

# import substitution matrix
EDNAFULL <- as.matrix(read.table("data/EDNAFULL", skip=8, header=T, row.names=1))

# Build function to find which reference sequences produces best alignment score
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
results <- data.frame(otu=names(otus), hits=best, stringsAsFactors = FALSE)
# If multiple hit, select first one
results$hit <- unlist(lapply(strsplit(results$hits, split=";;"), "[", 1))
results$hit <- gsub(";;", ";", paste(results$hit, ";", sep="")) # temp fix to end with ;....get rid of sizeout!!!

# Get percent identity of global alignment of each sequence with its best hit from database

# align each otu sequence to its best needle hit, using same parameters as command line needle
# "global-local" forces whole reference/pattern sequence to be used, but subject/query can be a concescutive subsequence -- this is because queries may include slightly longer region of gene than references because references were trimmed by o-smart-trim. this means that query sequence can have more bases than the reference at beginning or end of the sequence (=end gaps in alignment) without penalties.
psa <- pairwiseAlignment(subject=otus[results$otu], pattern=ref[results$hit], type="global-local",
                         substitutionMatrix=EDNAFULL,
                         gapOpening=10.0, gapExtension=0.5)

# calculate percent similarity to reference with indels counting as single difference regardless of length
results <- within(results, {
  hit.length <- width(pattern(psa))
  mismatches <- nmismatch(psa)
  insertions <- unlist(lapply(as.list(insertion(psa)), length))
  deletions <- unlist(lapply(as.list(deletion(psa)), length))
  similarity <- round((hit.length-mismatches-insertions-deletions)/hit.length, 4)
})

# Write results to .tsv file
write.table(results, file=paste(dirname(args[1]), "/", "nw_tophits.tsv", sep=""), 
            sep="\t", quote=FALSE)

