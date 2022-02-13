#### Kakapo 12S metabarcoding data processing script ####

## Load required packages
library(insect)
library(dada2)

## Set working directory to folder containing demuxed fastq files in a folder named 'demux', eDNA-matching and taxonomy rds files
setwd(dirname(file.choose()))

## DADA2 filtering
path <- "demux"
list.files(path)
fnFs <- sort(list.files(path, pattern="", full.names = TRUE))
sample.names <- basename(sub("\\.fastq", "", fnFs))
# check quality profiles using plotQualityProfile(fnFs[1:4])
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs,  truncLen=0,
                     maxN=0, maxEE = 2, truncQ = 2, rm.phix=TRUE,
                     compress = TRUE, multithread = T)
if(!all(file.exists(filtFs))){
  message("some samples had zero reads")
  message(sample.names[!file.exists(filtFs)])
  sample.names <- sample.names[file.exists(filtFs)]
  fnFs <- fnFs[file.exists(filtFs)]
  filtFs <- filtFs[file.exists(filtFs)]
  out <- filterAndTrim(fnFs, filtFs,  truncLen=0,
                       maxN=0, maxEE = 2, truncQ = 2, rm.phix=TRUE,
                       compress = TRUE, multithread = TRUE)
}
# learn error function
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF, file = "errF.rds")
derepFs <- derepFastq(filtFs, verbose=F)
names(derepFs) <- sample.names
# run dada function
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = FALSE)
# make sequence table
seqtab <- makeSequenceTable(dadaFs)
saveRDS(seqtab, file = "seqtab.rds")
# remove bimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# output statistics table
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track, 50) # gives you how many were dropped at each stage
write.csv(track, file = "trackprocess.csv", row.names = TRUE)
# save final sequence table
saveRDS(seqtab.nochim, file = "seqtab_nochim.rds")


## Create long-form sequence records table and save
myrecords <- vector(mode = "list", length = nrow(seqtab.nochim))
for(i in seq_len(nrow(seqtab.nochim))){
  tmp <- seqtab.nochim[i, ]
  if(is.null(names(tmp))) names(tmp) <- colnames(seqtab.nochim)
  tmp <- tmp[tmp > 0]
  if(length(tmp) > 0L){
    tmp <- sort(tmp, decreasing = TRUE)
    out <- names(tmp)
    UID <- as.integer(sub(".+_([[:digit:]]{6})$", "\\1", rownames(seqtab.nochim)[i]))
    newrecords <- data.frame(UID = rep(UID, length(tmp)))
    newrecords$PrimerSet <- "RV"
    newrecords$Sequence <- out
    newrecords$Count <- unname(tmp)
    myrecords[[i]] <- newrecords
  }
}
newrecords <- do.call("rbind", myrecords)
saveRDS(newrecords, file = "records.rds")


## Create table of unique sequences and adjoin taxon information and counts within each sample
records <- readRDS("records.rds")
eDNA <- readRDS("eDNA.rds")
taxa <- readRDS("taxa.rds")
newtb <- data.frame(Sequence = unique(records$Sequence))
newtb$Target <- "RV"
taxids <- eDNA$TaxID[match(newtb$Sequence, eDNA$Sequence)]
taxids[is.na(taxids)] <- 1
newtb$ScientificName <- taxa$name[match(taxids, taxa$taxID)]
newtb$Rank <- taxa$rank[match(taxids, taxa$taxID)]
newtb$TaxID <- taxids
newtb$CommonName <- taxa$commonname[match(taxids, taxa$taxID)]
countmat <- matrix(0L, nrow = nrow(newtb), ncol = length(uids))
colnames(countmat) <- uids
for(i in seq_along(uids)){
  tmp <- records[records$UID == uids[i], ]
  countmat[, i] <- tmp$Count[match(newtb$Sequence, tmp$Sequence)]
  countmat[, i][is.na(countmat[, i])] <- 0L
}
newtb <- cbind(newtb, countmat)
# order table by taxon rank and sequence count
seqcounts <- apply(newtb[grep("^[125][[:digit:]]{5}$", colnames(newtb))], 1,  sum)
newtb <- newtb[seqcounts > 0, ]
seqcounts <- seqcounts[seqcounts > 0]
seqranks <- rep(0L, nrow(newtb))
seqranks[!grepl("unranked", newtb$Rank)] <- 1L
seqranks[grepl("family", newtb$Rank)] <- 2L
seqranks[grepl("genus", newtb$Rank)] <- 3L
seqranks[grepl("species", newtb$Rank)] <- 4L
newtb <- newtb[order(seqranks, seqcounts, decreasing = TRUE),]


## Output sample metadata for upload to SRA
path <- "demux2"
fnames <- list.files(path)
out <- data.frame(SampleName = sub("WL...._.._(......)\\.fastq", "\\1", fnames),
                  LibraryID = sub("\\.fasta", "", fnames),
                  Filename = fnames)
write.csv(out, file = "../SRAmetadata.csv")

########################################
