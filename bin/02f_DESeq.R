setwd("~/Github/BehavEphyRNAseq/bin")
load("~/Github/BehavEphyRNAseq/bin/02b_kallisto_gather.Rdata")
library("DESeq2")
library("ggplot2")
library("dplyr")
library(magrittr)


# 1.3.3 Count matrix input ----
countData <- countbygene 
colData <- Traits %>%
  arrange(RNAseqID)

row.names(colData) <- colData$RNAseqID
colData$RNAseqID <- NULL

# factors must be factors
cols = c(1:4,6)
colData[,cols] %<>% lapply(function(x) as.factor(as.character(x)))

str(colData)
head(colData)
head(countData)
names(countData)
rownames(colData)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ APA + Punch + APA:Punch)
dds

#class: DESeqDataSet 
#dim: 22485 45 
#metadata(1): version
#assays(1): counts
#rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
#rowData names(0):
#  colnames(45): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
#colData names(7): Mouse Conflict ... APAconflict ID


## 1.3.6 Pre-filtering

dds <- dds[ rowSums(counts(dds)) > 1, ]

## 1.3.7 Note on factor levels

dds$APA <- factor(dds$APA, levels=c("Yoked","Trained"))
levels(dds$APA )

## 1.4  Differential expression analysi

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]

summary(res)
## APA CA1 and DG
#out of 17316 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 100, 0.58% 
#LFC < 0 (down)   : 34, 0.2% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 7426, 43% 
#(mean count < 6)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## APA ca1 ca2 and dg
#out of 18074 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 65, 0.36% 
#LFC < 0 (down)   : 21, 0.12% 
#outliers [1]     : 20, 0.11% 
#low counts [2]   : 5254, 29% 

# APA only
#out of 17688 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 32, 0.18% 
#LFC < 0 (down)   : 3, 0.017% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 11373, 64% 
#(mean count < 19)



sum(res$padj < 0.1, na.rm = TRUE)
# 134 apa ca1 dg
# 86 apa ca1 dg ca3
# 35 apa only

res05 <- results(dds, alpha=0.05)
summary(res05)

#out of 17316 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 73, 0.42% 
#LFC < 0 (down)   : 22, 0.13% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 10783, 62% 
#(mean count < 18)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## ca1 dg ca3
#out of 18074 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 50, 0.28% 
#LFC < 0 (down)   : 10, 0.055% 
#outliers [1]     : 20, 0.11% 
#low counts [2]   : 8021, 44% 

sum(res05$padj < 0.05, na.rm=TRUE)
# 95 ca1 dg
# 60 ca1 dg ca3

## 1.5 exploring and reporting results
plotMA(res, main="DESeq2", ylim=c(-1.5,1.5))

resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)

plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))


plotCounts(dds, gene=which.min(res$padj), intgroup="APA")
plotCounts(dds, gene=which.min(res$padj), intgroup="Punch")


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="APA",
                returnData=TRUE)

ggplot(d, aes(x=APA, y=count, color=APA)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

## 1.5 more info
mcols(res)$description


ddsMF <- dds


design(ddsMF) <- formula(~APA + Punch)

ddsMF <- DESeq(ddsMF)
resultsddsMF <- results(ddsMF)
head(resultsddsMF)
#resMFAPA <- results(ddsMF, contrast=c("APA","Trained", "Yoked"))
head(resMFAPA)

rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
head(assay(rld), 3)

# 2.2.1 Heatmap of the count matrix

install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:100]
head(select)

nt <- normTransform(dds) # defaults to log2(x+1) 
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("APA","Punch")])

pheatmap(log2.norm.counts, cluster_rows=TRUE,show_rownames=FALSE,
         show_colnames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
#pheatmap(assay(rld)[select], cluster_rows=FALSE)
#pheatmap(vsd.norm.counts, cluster_rows=FALSE)




library("RColorBrewer")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Punch, rld$APA, rld$Conflict, rld$Mouse, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pheatmap(cor(assay(rld)))

plotPCA(rld, intgroup=c("APA", "Punch"))

data <- plotPCA(rld, intgroup=c("APA", "Punch"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Punch, shape=APA)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


p78 <- fviz_pca_biplot(pc, axes = c(7, 8), geom = c("point"),
                       label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                       habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                       col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                       col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                       #col.circle = "grey70", 
                       select.var = list(contrib = 5), 
                       select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                       jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) 



## stats and individual genes
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

dds$group <- factor(paste0(dds$Punch, dds$APA))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
results(dds, contrast=c("group", "groupCA1Trained", "groupCA1Yoked"))

## 3.10 Access to all calculated values
deseq2results <- mcols(dds,use.names=TRUE)[1:4,]
write.csv(deseq2results, "deseq2results.csv", row.names = F)

# here using substr() only for display purposes
substr(names(mcols(dds)),1,10)

head(assays(dds)[["mu"]])



head(mcols(dds)$dispersion)

sizeFactors(dds)
head(coef(dds))
