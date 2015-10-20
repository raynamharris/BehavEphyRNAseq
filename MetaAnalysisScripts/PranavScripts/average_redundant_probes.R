## This function takes a microarray expression set vector for which the probe names have already been converted to gene ids and returns
## an eset for which the expression values of all probes for the same gene have been averaged.

average_redundant_probes_vec <- function (eset) {

  for (n in names(eset)) {
    red_probes <- eset[names(eset) %in% n]
    
    if (length(red_probes) > 1) {
      exp_avg <- mean(red_probes)
      eset <- eset[! names(eset) %in% n]
      eset <- append(eset, exp_avg)
      names(eset)[length(eset)] <- n
    }
  }
  return(eset)
}

## Analagous to the above for data frame expression sets, except the df argument takes a data frame named by probe id, not ENSEMBL id. 
## The reason is that data frames do not allow redundancy in column naming, so the mean expression calculation cannot happen before
## id conversion. Argument 'df' is a data frame of expression in the format observation by probe id, row by col. Argument 'id_vector'
## is a vector associating probe IDs to gene IDs, key to value.

average_redundant_probes_df <- function (data, id_vector, gene_by_obs = FALSE) {
  
  ## Import library for data frame manipulation
  library(data.table, lib.loc="~/R/R-3.2.1/library")
  
  ## For every probe, if it maps to more than one gene, add it to the mult_mapped_ids vector
  mult_mapped_ids <- c()
  for (probe_id in unique(names(id_vector))) {
    if (sum(names(id_vector) %in% probe_id) > 1) {
      mult_mapped_ids <- c(mult_mapped_ids, id_vector[names(id_vector) %in% probe_id])
    }
  }

  ## Eliminate multiple-mapped probes from the data set and conversion vector
  data <- data[!(colnames(data) %in% names(mult_mapped_ids))]
  id_vector <- id_vector[!(id_vector %in% mult_mapped_ids)] 

  ## Initialize a list to hold the rows of the output data frame.
  new_df_list <- list()

  ## Eliminate probes not in the conversion vector
  data <- data[unique(names(id_vector))]

  ## Iterate over each chip in the probe-level data and assign a subset of the data containing just the chip's row to df.
  for (repl in rownames(data)) {
    
    df <- data[repl,]
    
    ## Initialize a list. For every gene (key), add a vector containing all corresponding probes to the list (value)
    gene_list <- list()
    for (gene in unique(id_vector)) {
      gene_list[[gene]] <- names(id_vector[id_vector %in% gene])
    }

    ## For each gene, replace the probe vector with the average expression over the probes for every time point (vector)
    for (gene in names(gene_list)) {
      gene_list[[gene]] <- rowMeans(df[gene_list[[gene]]])
    }

    ## Create a new expression df from this gene list
    new_row <- as.data.frame(gene_list)

    ## Store this chip's gene-level expression data in the list of rows for the output data frame
    new_df_list[[repl]] <- new_row
  }
  
  ## rbind the rows into one output data frame, new_df, and return the result
  new_df <- rbindlist(new_df_list)
  rownames(new_df) <- names(new_df_list)
  return(new_df)
  
}