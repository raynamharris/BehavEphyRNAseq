# ## Generates the text for the covdesc file to import data from Xue, et al.
# stages <- c('TS01', 'TS09', 'TS11','TS13','TS16','TS19','TS21','TS22','TS23','TS25','TS27')
# chips <- list.files(path = 'H:/All projects/Phylotypic Stage/data/mus-TS01-TS27-micro-r-xue/')
# names(chips) <- rep(stages, each = 3)
# writeLines(paste(mapply(paste, chips, names(chips), MoreArgs = list(sep = '\t')), collapse = '\n'))

## Boot up data from the WGCNA tutorial
setwd('H:/Undergraduate Students/Pranav Bhamidipati/WGCNA Tutorial/');
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH;
rownames(datExpr0) = names(femData)[-c(1:8)];

## Import necessary methods
library("simpleaffy", lib.loc="~/R/R-3.2.1/library")
library("biomaRt", lib.loc="~/R/R-3.2.1/library")
source('H:/Undergraduate Students/Pranav Bhamidipati/Rscripts/average_redundant_probes.R')

## Import data
setwd('H:/All projects/Phylotypic Stage/data/mus-TS01-TS27-micro-r-xue/')
# raw.xue <- read.affy(covdesc = 'covdesc') # final form
raw.xue <- read.affy(covdesc = 'covdesc_test') # for testing
eset.xue <- call.exprs(raw.xue, algorithm = "rma", method = 'quantiles')
data.xue<- as(eset.xue, "data.frame")
colnames(data.xue) <- substr(colnames(data.xue), 2, nchar(colnames(data.xue)))
rownames(data.xue) <- c('rep1','rep2','rep3')

## Query BioMart database for the ENSEMBL gene IDs, excluding spike-in controls (prefixed "AFFX")
my_filters <- "affy_mouse430_2"
my_attributes <- c("affy_mouse430_2", "ensembl_gene_id")
my_values <- colnames(data.xue[grep("AFFX", colnames(data.xue), invert = TRUE)])
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene.ids.query <- getBM(attributes = my_attributes, filters = my_filters, values = my_values, mart = ensembl)

## Make a vector with ENSEMBL gene IDs as values and the corresponding Affy probe IDs as keys
gene.ids <- as.vector(gene.ids.query$ensembl_gene_id)
names(gene.ids) <- as.vector(gene.ids.query$affy_mouse430_2)

## Average expression over all probes for the same gene (internally omits probes that correspond to multiple gene IDs). 
results.xue <- average_redundant_probes_df(data.xue, gene.ids)