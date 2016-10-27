## Kallisto Gather was created by Dennis Wylie and edited by Rayna
## original kallisto gater script found here:
## location of script "/Volumes/HofmannLab/rmharris/singlecellseq/scripts.dennis"

setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444/06_kallistoquantcandidategenes/")
## this will create lists of all the samples
kallistoDirs = dir("~/Github/BehavEphyRNAseq/TACC-copy/JA16444/06_kallistoquantcandidategenes/")
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

# view transcripts per million dataframe
head(tpm)
tail(tpm)

## summarys: col and row sums
names(tpm)
str(tpm)
tpmsummary <- tpm %>% 
  mutate(sum = rowSums(.[2:29])) %>% 
  select(id, sum) %>% 
  filter(sum > 0)
tpmcolsums <- colSums(tpm[,-1])
tpmcolsums <- as.data.frame(tpmcolsums)

## summarys: col and row sums
names(counts)
str(counts)
countssummary <- counts %>% 
  mutate(sum = rowSums(.[2:29])) %>% 
  select(id, sum) %>% 
  filter(sum > 0)
countscolsums <- colSums(counts[,-1])
countscolsums <- as.data.frame(countscolsums)


## separating the compentents of the gene name
count$ENSMUST <- sapply(strsplit(as.character(count$id),'\\|'), "[", 1)
count$ENSMUSG <- sapply(strsplit(as.character(count$id),'\\|'), "[", 2)
count$OTTMUSG <- sapply(strsplit(as.character(count$id),'\\|'), "[", 3)
count$OTTMUST <- sapply(strsplit(as.character(count$id),'\\|'), "[", 4)
count$transcript <- sapply(strsplit(as.character(count$id),'\\|'), "[", 5)
count$gene <- sapply(strsplit(as.character(count$id),'\\|'), "[", 6)
count$length <- sapply(strsplit(as.character(count$id),'\\|'), "[", 7)
count$type <- sapply(strsplit(as.character(count$id),'\\|'), "[", 8)

## separating the compentents of the gene name
tpm$ENSMUST <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 1)
tpm$ENSMUSG <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 2)
tpm$OTTMUSG <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 3)
tpm$OTTMUST <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 4)
tpm$transcript <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 5)
tpm$gene <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 6)
tpm$length <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 7)
tpm$type <- sapply(strsplit(as.character(tpm$id),'\\|'), "[", 8)

head(tpm)

## summary by transcript
tpmsummary <- tpm %>% 
  mutate(sum = rowSums(.[2:29])) %>% 
  select(transcript, type, sum) %>% 
  filter(sum > 0)


