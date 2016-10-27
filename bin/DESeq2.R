## DESeq2 for differential gene expression analysis

#source("https://bioconductor.org/biocLite.R") ## needed to install deseq2
#biocLite("DESeq2") ## needed to install deseq2
library(DESeq2)
library(dplyr) ## for filtering and selecting rows
library(magrittr) ## to use the weird pipe

## set paths to sample info  and read in sample data
setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444")
samples <- read.csv("JA16444samples.csv", sep=",", header = TRUE, stringsAsFactors=FALSE, na.string = "NA")
str(samples)

#set reference levels (ie one that other will be compared to), default is alphabetical
#samples$Conflict <- relevel(samples$Conflict, ref="NoConflict") 
#samples$APA <- relevel(samples$APA, ref="Yoked") 

samples <- dplyr::rename(samples, sample = RNAseqID) ## unique identifier must be called sample
samples$APA[is.na(samples$APA)] <- "noAPA" ## make NA more meaningful
samples$Behavior[is.na(samples$Behavior)] <- "noAPA" ## make NA more meaningful
samples$E.phy[is.na(samples$E.phy)] <- "noAPA" ## make NA more meaningful
samples$Conflict[is.na(samples$Conflict)] <- "noAPA" ## make NA more meaningful

## add more grouping columns
samples$APAconflict <- as.factor(paste(samples$APA, samples$Conflict, sep="_")) 
samples$APApunch <- as.factor(paste(samples$APA, samples$Punch, sep="_")) 
samples$ConflictPunch <- as.factor(paste(samples$Conflict, samples$Punch, sep="_")) 
samples$APAconflictPunch <- as.factor(paste(samples$APA, samples$Conflict, samples$Punch, sep="_")) 
str(samples)
names(samples)

## make mutilple columsn factors, with 2 commands
cols = c(1:13)
samples[,cols] %<>% lapply(function(x) as.factor(as.character(x)))
str(samples)

## to calculate counts, use the kallisto_gather.R script (in the same directory) to calculate tpm and counts
## then, reset working directory and set either counts or
setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444")
counts <- count

## make "trainscript column" the rownames
rownames(counts) <- counts[,73]

## remove non ount data
counts <- counts[ -c(1,69:76) ]

## round the decimal values
counts <- round(counts)

##  remove two samples from both the samples and the ocunts file
samples <- samples %>% filter(!grepl("147D-CA1-1|145B-CA3-1", sample)) 
counts <- counts[ -c(36,56) ]


## construct DESeq database
d <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~APAconflictPunch )
head(d)

## calculate differential expression!
d <- DESeq(d)

# examine results results and also get some basic summary information
res <- results(d)
res
summary(res)

##output
#out of 142 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 16, 11% 
#low counts [2]   : 0, 0% 
#(mean count < 0)

# get the names for the specific results from the factors
resultsNames(d)
#[1] "Intercept"                              "APAconflictPunchnoAPA_noAPA_CA1"       
#[3] "APAconflictPunchnoAPA_noAPA_CA3"        "APAconflictPunchnoAPA_noAPA_DG"        
#[5] "APAconflictPunchTrained_Conflict_CA1"   "APAconflictPunchTrained_Conflict_CA3"  
#[7] "APAconflictPunchTrained_Conflict_DG"    "APAconflictPunchTrained_NoConflict_CA1"
#[9] "APAconflictPunchTrained_NoConflict_CA3" "APAconflictPunchTrained_NoConflict_DG" 
#[11] "APAconflictPunchYoked_Conflict_CA1"     "APAconflictPunchYoked_Conflict_CA3"    
#[13] "APAconflictPunchYoked_Conflict_DG"      "APAconflictPunchYoked_NoConflict_CA1"  
#[15] "APAconflictPunchYoked_NoConflict_CA3"   "APAconflictPunchYoked_NoConflict_DG"   
## examine reults for one factor
res <- results(d, name="APAconflictPunchYoked_NoConflict_CA3") 

# note you can change the p-value here to your desired cutoff
table(res$padj<0.1)

# make single file with LFC and pvals by pulling those variables and then binding them together (if you have other effects of interested you could employ a similar approach)
LFC = res$log2FoldChange
pval = res$pvalue
padj = res$padj
result = data.frame(cbind(LFC,pval,padj))
row.names(result) = row.names(res)
head(result)

# write your results into a single output file (for later, for other analyses, etc.)
#write.csv(as.data.frame(result), file="DESeq2Summary.csv") 


# let's also have a look at the data in a few different ways ...

# the MA plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts. points with an adjusted p-value < 0.1 are colored red. 
plotMA(res, main="DESeq2", ylim=c(-2,2))


# we may also want to look at expression differences in a specific gene (i.e. one gene at a time). we can call genes by name, but here we will just plot the gene with the greatest expression change (i.e. the minimum p-value)
plotCounts(d, gene=which.min(res$padj), intgroup="pop")

# PCA
rld <- rlog(d, blind=FALSE)
plotPCA(rld, intgroup=c("pop"))
