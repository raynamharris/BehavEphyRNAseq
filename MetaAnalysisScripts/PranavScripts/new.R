## This is my attempt to turn some of the code blocks in 9-9-2015-compressed.R into a library of functions

___ <- function (path, targets, )


## Import necessary methods
library("simpleaffy", lib.loc="~/R/R-3.2.1/library")
library("biomaRt", lib.loc="~/R/R-3.2.1/library")
source('H:/Undergraduate Students/Pranav Bhamidipati/Rscripts/average_redundant_probes.R')
#
## Read and pre-process the chips given in the targets file in the directory
setwd("H:/All projects/Phylotypic Stage/data/")
raw.roux <- read.affy(path="zfish-zyg-adult-micro-r-roux/")
eset.roux <- call.exprs(raw.roux,algorithm="gcrma")

## Compute FC between 2 time points
results.roux <- pairwise.comparison(eset.roux, "time.point")
fc.results.roux <- fc(results.roux)
