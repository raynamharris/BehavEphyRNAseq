## This script is a compressed workflow. This combines the importing and pre-processing of MA and RNA-seq data from zebrafish
## (Roux = MA; Yang = RNA-seq) by graphing their log2-fold-changes between two common time points parametrically against each other!
## The data can also be ranked by FC or mag(FC) and graphed.
## The second function performs a similar analysis with expression and log2(expression) values for each time point.

MA_RNA_FC_startup <- function() {

  ## Import necessary methods
  library('affy', lib.loc="~/R/R-3.2.1/library")
  library("simpleaffy", lib.loc="~/R/R-3.2.1/library")
  library("biomaRt", lib.loc="~/R/R-3.2.1/library")
  source('H:/Undergraduate Students/Pranav Bhamidipati/Rscripts/average_redundant_probes.R')
#
  ## Read and pre-process the chips given in the targets file in the directory
  setwd("H:/All projects/Phylotypic Stage/data/")
  raw.roux <- read.affy(path="zfish-zyg-adult-micro-r-roux/", covdesc = 'covdesc')
  eset.roux <- call.exprs(raw.roux,algorithm="gcrma")
  
  ## Compute FC between 2 time points
  results.roux <- pairwise.comparison(eset.roux, "time.point")
  fc.results.roux <- fc(results.roux)
#
  ## Query BioMart database for the ENSEMBL gene IDs, excluding spike-in controls (prefixed "AFFX")
  my_filters <- "affy_zebrafish"
  my_attributes <- c("affy_zebrafish", "ensembl_gene_id")
  my_values <- names(fc.results.roux[grep("AFFX", names(fc.results.roux), invert = TRUE)])
  ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
  gene.ids.query <- getBM(attributes = my_attributes, filters = my_filters, values = my_values, mart = ensembl)
  
  ## Make a vector with ENSEMBL gene IDs as values and the corresponding Affy gene IDs as keys
  ensembl.ids.roux <- as.vector(gene.ids.query$ensembl_gene_id)
  names(ensembl.ids.roux) <- as.vector(gene.ids.query$affy_zebrafish)
#
  ## Convert the probe IDs to ENSEMBL using the vector and omit any genes that don't convert
  names(fc.results.roux) <- ensembl.ids.roux[names(fc.results.roux)]
  fc.results.roux <- fc.results.roux[na.omit(names(fc.results.roux))]
#
  ## Average expression over probes with identical gene IDs
  fc.results.roux <- average_redundant_probes_vec(fc.results.roux)
#
  ## Read data from the RNA-seq study, stored in a csv
  eset.yang <- read.csv("zfish-2h-7d-RNAseq-p-yang-test.csv", stringsAsFactors = FALSE)
  eset.yang <- eset.yang[names(eset.yang) %in% c("ensembl_gene_id", "log2_fold_change")]
  fc.results.yang <- as.numeric(eset.yang$log2_fold_change)
  names(fc.results.yang) <- as.vector(eset.yang$ensembl_gene_id)
  fc.results.yang <- na.omit(fc.results.yang)
#
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

  
  

## Below are a bunch of different ways of plotting the data.
  
#   ## Plot RNA-seq vs. Microarray expression and log2(expression) for each time point
#   ## NOTE: x_ denotes Microarray, y_ denotes RNA-seq, _1 denotes time-point 1 (16hpf), _2 denotes time-point 2 (48hpf)
#   x1 <- as.double(results.roux[1,])
#   names(x1) <- colnames(results.roux)
#   x1 <- x1[order(names(x1))]
#   y1 <- as.double(results.yang[1,])
#   names(y1) <- colnames(results.yang)
#   y1 <- y1[order(names(y1))]
#   
#   plot(x1, y1, main = "RNA-seq vs. Microarray Expression, 16hpf", xlab = "Microarray Gene Expression (concentration)", 
#        ylab = "RNA-seq Gene Expression (RPKM)")
#   
#   x2 <- as.double(results.roux[2,])
#   names(x2) <- colnames(results.roux)
#   x2 <- x2[order(names(x2))]
#   y2 <- as.double(results.yang[2,])
#   names(y2) <- colnames(results.yang)
#   y2 <- y2[order(names(y2))]
#   
#   plot(x2, y2, main = "RNA-seq vs. Microarray Expression, 48hpf", xlab = "Microarray Gene Expression (concentration)", 
#        ylab = "RNA-seq Gene Expression (RPKM)")
#   
#   x1log2 <- log2(x1)
#   names(x1log2) <- names(x1)
#   y1log2 <- log2(y1)
#   names(y1log2) <- names(y1)
#   
#   plot(x1log2, y1log2, main = "RNA-seq vs. Microarray Log2(Expression), 16hpf", xlab = "Microarray Gene Expression Log2(concentration)", 
#        ylab = "RNA-seq Gene Expression Log2(RPKM)")
#   
#   x2log2 <- log2(x2)
#   names(x2log2) <- names(x2)
#   y2log2 <- log2(y2)
#   names(y2log2) <- names(y2)
#   
#   plot(x2log2, y2log2, main = "RNA-seq vs. Microarray Log2(Expression), 48hpf", xlab = "Microarray Gene Expression Log2(concentration)", 
#        ylab = "RNA-seq Gene Expression Log2(RPKM)")
  
#   ## Order by decreasing expression
#   x1 <- as.double(results.roux[1,])
#   y1 <- as.double(results.yang[1,])
#   names(x1) <- colnames(results.roux)
#   names(y1) <- colnames(results.yang)
#   x1 <- x1[order(x1, decreasing = T)]
#   y1 <- y1[order(y1, decreasing = T)]
#   
#   ## Replace expression with ranks and plot (NOTE: highest expression = rank of 1 = lowest rank)
#   x1[] <- seq(1, length(x1))
#   y1[] <- seq(1, length(y1))
#   x1 <- x1[order(names(x1))]
#   y1 <- y1[order(names(y1))]
#   
#   plot(x1, y1, main = "RNA-seq vs. Microarray Ranked Expression, 16hpf", xlab = "Microarray Gene Expression Rank", 
#        ylab = "RNA-seq Gene Expression Rank")
  
#   ## Order by decreasing expression
#   x2 <- as.double(results.roux[2,])
#   y2 <- as.double(results.yang[2,])
#   names(x2) <- colnames(results.roux)
#   names(y2) <- colnames(results.yang)
#   x2 <- x2[order(x2, decreasing = T)]
#   y2 <- y2[order(y2, decreasing = T)]
#   
#   ## Replace expression with ranks and plot (NOTE: highest expression = rank of 1 = lowest rank)
#   x2[] <- seq(1, length(x2))
#   y2[] <- seq(1, length(y2))
#   x2 <- x2[order(names(x2))]
#   y2 <- y2[order(names(y2))]
#   
#   plot(x2, y2, main = "RNA-seq vs. Microarray Ranked Expression, 48hpf", xlab = "Microarray Gene Expression Rank", 
#        ylab = "RNA-seq Gene Expression Rank")
  
#   ## Plot log2(expression) vs. rank(expression) for each technique at each time point
#   
#   ## Order by decreasing log2(expression). For RNA-seq, some RPKM values are zero, so the log2 transform produces values of -Inf. While this isn't inherently 
#   ## problematic, it impedes future computation with that data (ex. relative expression calculation). To avoid this, I replaced every zero value with the lowest
#   ## non-zero value in the vector, reasoning that the lowest detectable level of expression is comparable to no detectable expression.
#   x1 <- as.double(results.roux[1,])
#   y1 <- as.double(results.yang[1,])
#   x2 <- as.double(results.roux[2,])
#   y2 <- as.double(results.yang[2,])
#   x1log2 <- log2(x1)
#   y1[y1 == 0] <- min(y1[y1 != 0])
#   y1log2 <- log2(y1)
#   x2log2 <- log2(x2)
#   y2[y2 == 0] <- min(y2[y2 != 0])
#   y2log2 <- log2(y2)
#   names(x1log2) <- colnames(results.roux)
#   names(y1log2) <- colnames(results.yang)
#   names(x2log2) <- colnames(results.roux)
#   names(y2log2) <- colnames(results.yang)
#   x1log2 <- x1log2[order(x1log2, decreasing = T)]
#   y1log2 <- y1log2[order(y1log2, decreasing = T)]
#   x2log2 <- x2log2[order(x2log2, decreasing = T)]
#   y2log2 <- y2log2[order(y2log2, decreasing = T)]
#   
#   ## Create another vector with ranked expression (NOTE: highest expression = rank of 1 = lowest rank)
#   x1rank <- x1log2
#   y1rank <- y1log2
#   x2rank <- x2log2
#   y2rank <- y2log2
#   x1rank[] <- seq(1, length(x1log2))
#   y1rank[] <- seq(1, length(y1log2))
#   x2rank[] <- seq(1, length(x2log2))
#   y2rank[] <- seq(1, length(y2log2))
#   
#   ## Order all the vectors alphabetically by gene ID
#   x1log2 <- x1log2[order(names(x1log2))]
#   y1log2 <- y1log2[order(names(y1log2))]
#   x2log2 <- x2log2[order(names(x2log2))]
#   y2log2 <- y2log2[order(names(y2log2))]
#   x1rank <- x1rank[order(names(x1rank))]
#   y1rank <- y1rank[order(names(y1rank))]
#   x2rank <- x2rank[order(names(x2rank))]
#   y2rank <- y2rank[order(names(y2rank))]
# 
#   
#   ## For each time point, plot log2(expr) vs. rank for MA (blue) and RNA-seq (red)
#   plot(y1rank, y1log2, col = "red", 
#        main = "Microarray and RNA-seq: Log2(Expression) vs. Ranked Expression, 16hpf", 
#        xlab = "Expression Rank (high to low)", 
#        ylab = "Log2(Expression)")
#   points(x1rank, x1log2, col = "blue")
# 
#   plot(y2rank, y2log2, col = "red", 
#        main = "Microarray and RNA-seq: Log2(Expression) vs. Ranked Expression, 48hpf", 
#        xlab = "Expression Rank (high to low)", 
#        ylab = "Log2(Expression)")
#   points(x2rank, x2log2, col = "blue")
  
#   ## Compute "relative expression" from log2(Expression) for each time point (see Zhao, et al, Jan. 2014, PLoS One). Relative 
#   ## expression seems to be a measure of where the expression stands in the range of data (ignoring center and spread). In the 
#   ## paper they average log2(Expr) over all time points, but I think it makes more sense to view this at each time point.
#   x1relex <- (x1log2 - min(x1log2)) / (max(x1log2) - min(x1log2))
#   y1relex <- (y1log2 - min(y1log2)) / (max(y1log2) - min(y1log2))
#   x2relex <- (x2log2 - min(x2log2)) / (max(x2log2) - min(x2log2))
#   y2relex <- (y2log2 - min(y2log2)) / (max(y2log2) - min(y2log2))
#
#   ## Plot log2(Expression) vs. relative expression
#   plot(y1relex, y1log2, col = "red", 
#        main = "Microarray and RNA-seq: Log2(Expression) vs. Ranked Expression, 16hpf", 
#        xlab = "Relative Expression", 
#        ylab = "Log2(Expression)")
#   points(x1relex, x1log2, col = "blue")
#   
#   plot(y2relex, y2log2, col = "red", 
#        main = "Microarray and RNA-seq: Log2(Expression) vs. Ranked Expression, 48hpf", 
#        xlab = "Relative Expression", 
#        ylab = "Log2(Expression)")
#   points(x2relex, x2log2, col = "blue")
#   ## The previous lines generated plots that looked linear, very different from the published plots in Zhao, et al, which were
#   ## non-linear). Still very unsure what they mean by "relative expression."

}