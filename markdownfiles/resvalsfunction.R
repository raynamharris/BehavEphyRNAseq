resvals <- function(contrastvector, mypadj){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = T)
  sum <- sum(res$padj < mypadj, na.rm = TRUE)
  print(sum)
  vals <- cbind(res$pvalue, res$padj)
  pvalcolname <- as.character(paste("pval",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  padjcolname <- as.character(paste("padj",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  colnames(vals) <- c(pvalcolname, padjcolname)
  return(vals)
}



