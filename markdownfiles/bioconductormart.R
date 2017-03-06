# install the biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# load biomaRt
library(biomaRt)

# look at top 10 databases
head(listMarts(host = "www.ensembl.org"), 10)


head(listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")), 10)

head(listAttributes(useDataset(dataset = "hsapiens_gene_ensembl", 
                               mart    = useMart("ENSEMBL_MART_ENSEMBL",host = "www.ensembl.org"))), 10)

head(listFilters(useDataset(dataset = "hsapiens_gene_ensembl", 
                            mart    = useMart("ENSEMBL_MART_ENSEMBL",
                                              host = "www.ensembl.org"))), 10)

# 1) select a mart and data set
mart <- useDataset("hsapiens_gene_ensembl", 
                   mart = useMart("ENSEMBL_MART_ENSEMBL",
                                  host = "www.ensembl.org"))

# 2) run a biomart query using the getBM() function
# and specify the attributes and filter arguments
geneSet <- "GUCA2A"

resultTable <- getBM(attributes = c("start_position","end_position","description"),
                     filters = "hgnc_symbol", values = geneSet, mart = mart)
resultTable 




# load the biomartr package
biocLite("biomartr")
library(biomartr)

# list all available databases
getMarts()
