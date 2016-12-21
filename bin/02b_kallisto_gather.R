## Kallisto Gather was created by Dennis Wylie and edited by Rayna
## original kallisto gater script found here:
## location of script "/Volumes/HofmannLab/rmharris/singlecellseq/scripts.dennis"
library("dplyr")

setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444/04_kallistoquant/")
## this will create lists of all the samples
kallistoDirs = dir(".")
kallistoDirs = kallistoDirs[!grepl("\\.(R|py|pl|sh|xlsx?|txt|tsv|csv|org|md|obo|png|jpg|pdf)$",
        kallistoDirs, ignore.case=TRUE)]

kallistoFiles = paste0(kallistoDirs, "/abundance.tsv")
names(kallistoFiles) = kallistoDirs
if(file.exists(kallistoFiles))
  kallistoData = lapply(
    kallistoFiles,
    read.table,
    sep = "\t",
    row.names = 1,
    header = TRUE
)

## this for loop uses the reduce function to make two data frame with counts or tpm from all the samples
ids = Reduce(f=union, x=lapply(kallistoData, rownames))
if (all(sapply(kallistoData, function(x) {all(rownames(x)==ids)}))) {
    count = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$est_counts}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    tpm = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$tpm}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

## make a dataframe with the parts of the gene id as columns
geneids <- count[c(1)] 
geneids$ENSMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 1)
geneids$ENSMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 2)
geneids$OTTMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 3)
geneids$OTTMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 4)
geneids$transcript <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 5)
geneids$gene <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 6)
geneids$length <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 7)
geneids$structure1 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 8)
geneids$structure2 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 9)
geneids$structure3 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 10)
geneids$transcript_lenght <- as.factor(paste(geneids$transcript, geneids$length, sep="_"))

## prep data for wgcna
countswgcna <- count
row.names(countswgcna) <- geneids$transcript_lenght
countswgcna[1] <- NULL
countswgcna <- round(countswgcna)
summary(countswgcna)

## prep data for wgcna
tpmswgcna <- tpm
row.names(tpmswgcna) <- geneids$transcript_lenght
tpmswgcna[1] <- NULL
tpmswgcna <- round(tpmswgcna)
summary(tpmswgcna)

## save output file for easy loading later
#setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444")
#write.csv(count, "count.csv", row.names = TRUE)
#write.csv(countswgcna, "countswgcna.csv", row.names = TRUE)
#write.csv(tpm, "tpm.csv", row.names = TRUE)
#write.csv(tpmswgcna, "tpmswgcna.csv", row.names = TRUE)

# save Rdata for loadinhg
#save.image("~/Github/BehavEphyRNAseq/bin/02b_kallisto_gather.Rdata")
