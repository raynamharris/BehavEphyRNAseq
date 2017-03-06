Not all of my samples need to be processed together. This script will
make the appropriate subset files. First, read in the new data

    library("dplyr") ## for filtering and selecting rows

    DissociationCountData <- read.csv("../data/rnaseq/countbygene.csv", header=T, check.names = F, row.names = 1)
    DissociationColData <- read.csv("../data/rnaseq/kallistosamples.csv", header=T)

    DissociationColData <- DissociationColData %>%
      filter(Group %in% c("homecage", "control")) %>%
      filter(Genotype != "FMR1") %>%
      filter(!grepl("C|D|16-", Mouse))
    savecols <- as.character(DissociationColData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    DissociationCountData <- DissociationCountData %>% select(one_of(savecols)) # keep good samples

    # save files to dissociation project directory
    write.csv(DissociationColData, "~/Github/DissociationTest/data/DissociationColData.csv", row.names=F)
    write.csv(DissociationCountData, "~/Github/DissociationTest/data/DissociationCountData.csv", row.names=T)

Session Info
------------

    sessionInfo()

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: OS X 10.10.5 (Yosemite)
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_0.5.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.7     digest_0.6.11   rprojroot_1.2   assertthat_0.1 
    ##  [5] R6_2.2.0        DBI_0.5-1       backports_1.0.5 magrittr_1.5   
    ##  [9] evaluate_0.10   stringi_1.1.2   lazyeval_0.2.0  rmarkdown_1.3  
    ## [13] tools_3.3.1     stringr_1.1.0   yaml_2.1.14     htmltools_0.3.5
    ## [17] knitr_1.15.1    tibble_1.2
