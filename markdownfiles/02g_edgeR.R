#source("http://www.bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)

counts <- countData
dim( counts )
colSums( counts )

lowcounts <- as.data.frame(table( rowSums( counts ) )[ 1:30 ]) # Number of genes with low counts

millioncounts <- as.data.frame(colSums( counts ) / 1e06)  # in millions of reads
colnames(millioncounts) <- c("totalmilcounts")
colSums(millioncounts)
millioncounts$fracoftotal <- round((millioncounts$totalmilcounts / colSums(millioncounts)), digits = 2)
millioncounts$fracoftotalminus1 <- round((millioncounts$totalmilcounts / ((colSums(millioncounts)) - 11)), digits = 2)


group <- c(rep("C", 4) , rep("T", 3))
colData$Method

cds <- DGEList( counts , group = colData$Method)
names( cds )
head(cds$counts)
cds$samples
sum( cds$all.zeros )
cds

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
dim( cds )

cds <- calcNormFactors( cds )
cds$samples

cds$samples$lib.size*cds$samples$norm.factors
# broken
# plotMDS.dge( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )

cds <- estimateCommonDisp( cds )
cds$common.dispersion

cds <- estimateTagwiseDisp( cds , prior.n = 10 )
summary( cds$tagwise.dispersion )

meanVarPlot <- plotMeanVar( cds ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
                            pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ) #these arguments are to make it look prettier
