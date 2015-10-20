## This script imports and pre-processes of microarray data from multiple hippocampus datasets. 
## The goal is to identify learning-induced gene modules that are preserved in single CA1 neurons, CA1 tissue samples, and whole hippocampi.
## The data can also be ranked by fold change (FC) or mag(FC) and graphed.
## The second function performs a similar analysis with expression and log2(expression) values for each time point.

# source("https://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("simpleaffy")
# biocLite("biomaRt")

## Import necessary methods
library('affy')
library("simpleaffy")
library("biomaRt")

## Read and pre-process the chips given in the targets file in the directory
setwd("~/Dropbox/BehavEphyRNAseq/MetaAnalysisData")
raw.haberman08 <- read.affy(path="Haberman08/", covdesc = 'covdesc_DG')
eset.haberman08 <- call.exprs(raw.haberman08,algorithm="gcrma")
  
## Compute Fold Change (FC) between control and learning groups
results.haberman08 <- pairwise.comparison(eset.haberman08, "Behavior")
fc.results.haberman08 <- fc(results.haberman08)

## Query BioMart database for the ENSEMBL gene IDs, excluding spike-in controls (prefixed "AFFX")
listMarts() #check available databases
#listDatasets("ensembl") #check species-specific ensembl databases
ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
# assigning filters, attributes, and values for BioMart query
my_filters <- "affy_rae230a" #calling the Rat array specified in GEO
my_attributes <- c("affy_rae230a", "ensembl_gene_id")
my_values <- names(fc.results.haberman08[grep("AFFX", names(fc.results.haberman08), invert = TRUE)])
gene.ids.query <- getBM(attributes = my_attributes, filters = my_filters, values = my_values, mart = ensembl)
  
## Make a vector with ENSEMBL gene IDs as values and the corresponding Affy gene IDs as keys
ensembl.ids.haberman08 <- as.vector(gene.ids.query$ensembl_gene_id)
names(ensembl.ids.haberman08) <- as.vector(gene.ids.query$affy_rae230a)


## Convert the probe IDs to ENSEMBL using the vector and omit any genes that don't convert
names(fc.results.haberman08) <- ensembl.ids.haberman08[names(fc.results.haberman08)]
fc.results.haberman08 <- fc.results.haberman08[na.omit(names(fc.results.haberman08))]
  
#looking at our current data
plot(fc.results.haberman08, main="DG")
sorted <- sort(fc.results.haberman08, decreasing = TRUE)
head(sorted)

## Average expression over probes with identical gene IDs
# not doing this yet!
# need to Add file average_redundant_probes.r to source
source("C:/Users/Bridget Kajs/Documents/Hofmann Lab/Hippocampal Gene Expression/Data for Meta-Analysis")
fc.results.haberman08 <- average_redundant_probes_vec(fc.results.haberman08)

## Read data from the RNA-seq study, stored in a csv
eset.yang <- read.csv("zfish-2h-7d-RNAseq-p-yang-test.csv", stringsAsFactors = FALSE)
eset.yang <- eset.yang[names(eset.yang) %in% c("ensembl_gene_id", "log2_fold_change")]
fc.results.yang <- as.numeric(eset.yang$log2_fold_change)
names(fc.results.yang) <- as.vector(eset.yang$ensembl_gene_id)
fc.results.yang <- na.omit(fc.results.yang)

## Restrict each data set to the genes they have in common
common.genes <- colnames(fc.results.roux)[colnames(fc.results.roux) %in% colnames(fc.results.yang)]
fc.results.roux <- fc.results.roux[common.genes]
fc.results.yang <- fc.results.yang[common.genes]
#
by_mag_FC <- userQuery("Rank genes by magnitude of FC instead of FC?", allowed = c("y", "n"), default = "n", case.sensitive = FALSE)
  
if (by_mag_FC == 'y') {
    
    ## This will order genes by the magnitude of FC (high to low), useful for manipulating and comparing gene lists
  Microarray_Rank <- abs(common.roux)
  Microarray_Rank <- Microarray_Rank[order(Microarray_Rank, decreasing = TRUE)]
  RNA_seq_Rank <- abs(common.yang)
  RNA_seq_Rank <- RNA_seq_Rank[order(RNA_seq_Rank, decreasing = TRUE)]
    
  ## Rank genes in order and sort alphabetically
  Microarray_Rank[] <- seq(1, length(Microarray_Rank))
  Microarray_Rank <- Microarray_Rank[order(names(Microarray_Rank))]
  RNA_seq_Rank[] <- seq(1, length(RNA_seq_Rank))
  RNA_seq_Rank <- RNA_seq_Rank[order(names(RNA_seq_Rank))]
    
  ## Plot RNA-seq Abs(FC) vs MA Abs(FC) between the two time-points
  plot(Microarray_Rank, RNA_seq_Rank, main = "RNA-seq vs. Microarray Rank Abs(FC)", xlab = "Microarray Gene Rank", 
       ylab = "RNA-seq Gene Rank")
    
} else {
    ## This will order genes by FC (low to high), useful for comparing MA and RNA-seq expression to ranked expression
  Microarray_Rank <- common.roux
  Microarray_Rank <- Microarray_Rank[order(Microarray_Rank)]
  RNA_seq_Rank <- common.yang
  RNA_seq_Rank <- RNA_seq_Rank[order(RNA_seq_Rank)]
  
  ## Rank genes in order and sort alphabetically
  Microarray_Rank[] <- seq(1, length(Microarray_Rank))
  Microarray_Rank <- Microarray_Rank[order(names(Microarray_Rank))]
  RNA_seq_Rank[] <- seq(1, length(RNA_seq_Rank))
  RNA_seq_Rank <- RNA_seq_Rank[order(names(RNA_seq_Rank))]
  
  ## Plot RNA-seq FC vs MA FC between the two time-points
  plot(Microarray_Rank, RNA_seq_Rank, main = "RNA-seq vs. Microarray FC Rank", xlab = "Microarray Gene Rank", 
       ylab = "RNA-seq Gene Rank")
}
  

## This script takes MA_RNA_startup() and adapts it to plot expression vs expression rather than FC or mag(FC)

MA_RNA_expr_startup <- function() {
  
  ## Import necessary methods
  library("simpleaffy", lib.loc="~/R/R-3.2.1/library")
  library("biomaRt", lib.loc="~/R/R-3.2.1/library")
  source('H:/Undergraduate Students/Pranav Bhamidipati/Rscripts/average_redundant_probes.R')
#
  ## Read and pre-process the chips given in the targets file in the directory
  setwd("H:/All projects/Phylotypic Stage/data/")
  raw.roux <- read.affy(path="zfish-zyg-adult-micro-r-roux/")
  eset.roux <- call.exprs(raw.roux,algorithm="gcrma")
  data.roux <- as(eset.roux, "data.frame")
#
  ## Average the expression of replicate chips for each time point, assuming 2 chips per time point and rows sorted by time point
  ## (This is kind of a crappy assumption that only works for this data set. Maybe I'll write a better function at some point.)
  n_tp <- dim(data.roux)[1] %/% 2
  for (i in seq(n_tp)) {
    x <- colMeans(rbind(as.numeric(data.roux[2*i - 1,]), as.numeric(data.roux[2*i,])))
    data.roux <- rbind(data.roux, x)
    row.names(data.roux)[nrow(data.roux)] <- i
  }
  seq_rows <- seq(nrow(data.roux))
  data.roux <- subset(data.roux, seq_rows > n_tp * 2)
#
  ## Query BioMart database for the ENSEMBL gene IDs, excluding spike-in controls (prefixed "AFFX")
  my_filters <- "affy_zebrafish"
  my_attributes <- c("affy_zebrafish", "ensembl_gene_id")
  my_values <- colnames(data.roux[grep("AFFX", colnames(data.roux), invert = TRUE)])
  ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
  gene.ids.query <- getBM(attributes = my_attributes, filters = my_filters, values = my_values, mart = ensembl)
  
  ## Make a vector with ENSEMBL gene IDs as values and the corresponding Affy gene IDs as keys
  ensembl.ids.roux <- as.vector(gene.ids.query$ensembl_gene_id)
  names(ensembl.ids.roux) <- as.vector(gene.ids.query$affy_zebrafish)
#
  ## Average expression over probes for the same gene (internally omits probes that correspond to multiple or zero gene IDs). 
  results.roux <- average_redundant_probes_df(data.roux, ensembl.ids.roux)
#
  ## Read data from the RNA-seq study, stored in a csv. Convert expression for the two time points to a data frame results.yang
  eset.yang <- read.csv("zfish-2h-7d-RNAseq-p-yang-test.csv", stringsAsFactors = FALSE)
  results.yang <- list()
  for (i in 1:nrow(eset.yang)) {
    results.yang[[eset.yang$ensembl_gene_id[i]]] <- c(eset.yang$X15.somite[i], eset.yang$X48hpf[i])
  }
  results.yang <- as.data.frame(results.yang)
  rownames(results.roux) <- c('16 hpf', '48 hpf')
#
  ## Restrict each data set to the genes they have in common
  common.genes <- colnames(results.roux)[colnames(results.roux) %in% colnames(results.yang)]
  results.roux <- results.roux[common.genes]
  results.yang <- results.yang[common.genes]