setwd("~/Github/BehavEphyRNAseq/bin")

## if loadinging, proceed to line 102
load("~/Github/BehavEphyRNAseq/bin/deseq.Rdata")

library("DESeq2")
library("ggplot2")
library("dplyr")
library("magrittr")


# 1.3.3 Count matrix input ----
countData <- countbygene 
colData <- Traits %>%
  arrange(RNAseqID)
head(countData)
head(colData)

countData[countData < 2] <- 0

## remove outliers
# 146B-DG, 146D-DG, 148B-CA1
dropcols <- c("146B-DG-2","146D-DG-3", "148B-CA3-4", "146C-CA3-4")
countData <- countData %>% select(-one_of(dropcols))
colData <- colData %>% filter(!grepl("146B-DG-2|146D-DG-3|148B-CA3|146C-CA3-4", RNAseqID))

row.names(colData) <- colData$RNAseqID
colData$RNAseqID <- NULL

# factors must be factors
cols = c(1:4,7)
colData[,cols] %<>% lapply(function(x) as.factor(as.character(x)))
colData$Slice <- as.factor(colData$Slice)
str(colData)


dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Punch + TrainGroup + Punch*TrainGroup)
dds

#class: DESeqDataSet 
#dim: 22485 45 
#metadata(1): version
#assays(1): counts
#rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
#rowData names(0):
#  colnames(42): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
#colData names(7): Mouse Conflict ... APAconflict ID


## 1.3.6 Pre-filtering

dds <- dds[ rowSums(counts(dds)) > 1, ]

## 1.3.7 Note on factor levels
str(dds$TrainGroup)
levels(dds$TrainGroup)
str(dds$Punch)
dds$TrainGroup <- factor(dds$TrainGroup, levels=c("Yoked","Trained"))
dds$Punch <- factor(dds$Punch, levels=c("DG","CA3", "CA1"))


## 1.4  Differential expression analysi

dds <- DESeq(dds)

# general deseq
res <- results(dds, independentFiltering = F)
resOrdered <- res[order(res$padj),]
summary(res)
sum(res$padj < 0.1, na.rm = TRUE) #14
res05 <- results(dds, alpha=0.05)
summary(res05) #23
sum(res05$padj < 0.05, na.rm=TRUE)

## 1.5 exploring and reporting results

plotMA(res, main="DESeq2")
resMLE <- results(dds)
head(resMLE, 4)

plotCounts(dds, gene=which.min(res$padj), intgroup="TrainGroup")
plotCounts(dds, gene=which.min(res$padj), intgroup="Punch")


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="TrainGroup",
                returnData=TRUE)

ggplot(d, aes(x=TrainGroup, y=count, color=TrainGroup)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

## 1.5 more info
mcols(res)$description

## for variance stablized gene expression and log transformed data
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
head(assay(rld), 3)

#####save data space!! ####
#save.image("~/Github/BehavEphyRNAseq/bin/deseq.Rdata")

## DEG by contrasts
resAPATY <- results(dds, contrast = c("TrainGroup", "Yoked", "Trained"), independentFiltering = F)
sum(resAPATY$padj < 0.1, na.rm = TRUE) #118

resPunchCA1DG <- results(dds, contrast = c("Punch", "CA1", "DG"), independentFiltering = F)
sum(resPunchCA1DG$padj < 0.1, na.rm = TRUE) #3219
resPunchCA1CA3 <- results(dds, contrast = c("Punch", "CA1", "CA3"), independentFiltering = F)
sum(resPunchCA1CA3$padj < 0.1, na.rm = TRUE) #2266
resPunchCA3DG <- results(dds, contrast = c("Punch", "CA3", "DG"), independentFiltering = F)
sum(resPunchCA3DG$padj < 0.1, na.rm = TRUE) #4029

##  now bind the table of 
valsAPATY <- cbind(resAPATY$pvalue, resAPATY$padj) 
valsPunchCA1DG <- cbind(resPunchCA1DG$pvalue, resPunchCA1DG$padj) 
valsPunchCA1CA3 <- cbind(resPunchCA1CA3$pvalue, resPunchCA1CA3$padj) 
valsPunchCA3DG <- cbind(resPunchCA3DG$pvalue, resPunchCA3DG$padj) 

colnames(valsAPATY)=c("pval.APATY", "padj.APATY")
colnames(valsPunchCA1DG)=c("pval.CA1DG", "padj.CA1DG")
colnames(valsPunchCA1CA3)=c("pval.CA1CA3", "padj.CA1CA3")
colnames(valsPunchCA3DG)=c("pval.CA3DG", "padj.CA3DG")



rldd <- assay(rld)
rldpvals <- cbind(rldd, valsAPATY, valsPunchCA1DG, valsPunchCA1CA3, valsPunchCA3DG)
head(rldpvals)
dim(rldpvals)
#14347    56
table(complete.cases(rldpvals))
#   TRUE 
#  17422 

### venn diagram
#install.packages("VennDiagram")
library(VennDiagram)

rldpvals <- as.data.frame(rldpvals)

APATY <- row.names(rldpvals[rldpvals$padj.APATY<0.1 & !is.na(rldpvals$padj.APATY),])
PunchCA1DG <- row.names(rldpvals[rldpvals$padj.CA1DG<0.1 & !is.na(rldpvals$padj.CA1DG),])
PunchCA1CA3 <- row.names(rldpvals[rldpvals$padj.CA1CA3<0.1 & !is.na(rldpvals$padj.CA1CA3),])
PunchCA3DG <- row.names(rldpvals[rldpvals$padj.CA3DG<0.1 & !is.na(rldpvals$padj.CA3DG),])

candidates <- list("CA1 v. DG" = PunchCA1DG, "CA1 v. CA3" = PunchCA1CA3,"CA3 v. DG"= PunchCA3DG, "Trained v. Yoked" = APATY)

dev.off()
prettyvenn <- venn.diagram(
  x = candidates, filename=NULL, lwd=4,
  col = "transparent",
  fill = (values=c("#00441b", "#238b45", "#66c2a4", "#ccece6")),
  alpha = 0.5,
  cex = 1, fontfamily = "sans", #fontface = "bold",
  cat.default.pos = "text",
  #cat.col = c("darkred", "darkgreen", "blue4", "orange"),
  #cat.dist = c(0.08, 0.08, 0.08, 0.08), cat.pos = 1,
  cat.cex = 1, cat.fontfamily = "sans"
);
grid.draw(prettyvenn)


#save_plot("rnaseqvennV2.pdf", prettyvenn, base_aspect_ratio = 1.4)
#save_plot("rnaseqvennV2.png", prettyvenn, base_aspect_ratio = 1.4)


# 2.2.1 Heatmap of the count matrix
#install.packages("pheatmap")
library("pheatmap")
nt <- normTransform(dds) # defaults to log2(x+1) 

## making pheatmaps annoations
df <- as.data.frame(colData(dds)[,c("Punch","TrainGroup")])
str(df)

ann_colors = list(
  TrainGroup =  c(Yoked = (values=c("#f1a340")), Trained = (values=c("#9970ab"))),
  Punch =  c(DG = (values=c("#006837")),  CA3 = (values=c("#41ab5d")), 
             CA1 = (values=c("#d9f0a3"))))

#### top variance genes des by brain region
DEGes <- as.data.frame(rldpvals) # convert matrix to dataframe
DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
DEGes$padjmin <- with(DEGes, pmin(padj.APATY, padj.CA1DG, padj.CA1CA3, padj.CA3DG)) # put the min pvalue in a new column
DEGes <- DEGes %>% filter(padjmin < 0.001)
rownames(DEGes) <- DEGes$rownames
drop.cols <- c("padj.APATY", "padj.CA1DG" , "padj.CA1CA3" , "padj.CA3DG" , "pval.APATY", "pval.CA1DG" , "pval.CA1CA3" , "pval.CA3DG" , "rownames", "padjmin")
DEGes <- DEGes %>% select(-one_of(drop.cols))
DEGes <- as.matrix(DEGes)
DEGes <- DEGes - rowMeans(DEGes)


pheatmap(DEGes, show_colnames=TRUE, show_rownames = F,
         annotation_col=df, annotation_colors = ann_colors,
         fontsize = 12, fontsize_row = 10, 
         #cellwidth=10, cellheight=10,
         #width = 10,
         border_color = "grey60"
)

DEGes <- as.data.frame(rldpvals) # convert matrix to dataframe
DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
DEGes$padjmin <- with(DEGes, pmin(padj.APATY)) # put the min pvalue in a new column
DEGes <- DEGes %>% filter(padjmin < 0.001)
rownames(DEGes) <- DEGes$rownames
drop.cols <- c("padj.APATY", "padj.CA1DG" , "padj.CA1CA3" , "padj.CA3DG" , "pval.APATY", "pval.CA1DG" , "pval.CA1CA3" , "pval.CA3DG" , "rownames", "padjmin")
DEGes <- DEGes %>% select(-one_of(drop.cols))
DEGes <- as.matrix(DEGes)
DEGes <- DEGes - rowMeans(DEGes)


pheatmap(DEGes, show_colnames=TRUE, show_rownames = F,
         annotation_col=df, annotation_colors = ann_colors,
         fontsize = 12, fontsize_row = 10, 
         #cellwidth=10, cellheight=10,
         #width = 10,
         border_color = "grey60"
)


## regular pca
plotPCA(rld, intgroup=c("TrainGroup", "Punch"), returnData=TRUE)
pcadata <- plotPCA(rld, intgroup=c("TrainGroup", "Punch"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))

ggplot(pcadata, aes(PC1, PC2, color=Punch, shape=TrainGroup)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  stat_ellipse(level = 0.95, (aes(color=Punch)),size=1.5)   + 
  scale_color_manual(values=c("#006837", "#41ab5d", "#d9f0a3")) +
  geom_text(aes(label=name),vjust=2)


## stats
library(edgeR)
counts <- countData
dim( counts )
colSums( counts ) / 1e06  # in millions of reads
table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

rowsum <- as.data.frame(colSums( counts ) / 1e06 )
names(rowsum)[1] <- "millioncounts"
rowsum$sample <- row.names(rowsum)

ggplot(rowsum, aes(x=millioncounts)) + 
  geom_histogram(binwidth = 1, colour = "black", fill = "darkgrey") +
  theme_classic() +
  scale_x_continuous(name = "Millions of Gene Counts per Sample",
                     breaks = seq(0, 8, 1),
                     limits=c(0, 8)) +
  scale_y_continuous(name = "Number of Samples")
 