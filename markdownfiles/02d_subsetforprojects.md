Not all of my samples need to be processed together. This script will
make the appropriate subset files. First, read in the new data

    library("dplyr") ## for filtering and selecting rows
    library("plyr") ## for renaming factors

    DissociationCountData <- read.csv("../data/rnaseq/countbygene.csv", header=T, check.names = F, row.names = 1)
    DissociationColData <- read.csv("../data/rnaseq/kallistosamples.csv", header=T)

    DissociationColData <- DissociationColData %>%
      filter(Group %in% c("homecage", "control")) %>%
      filter(Genotype != "FMR1") %>%
      filter(Mouse == "15-100") %>%
      filter(!grepl("C|D|16-", Mouse))  
    savecols <- as.character(DissociationColData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    DissociationCountData <- DissociationCountData %>% select(one_of(savecols)) # keep good samples

    # For GEO database
    DissociationMetaData <- DissociationColData %>%
      select(RNAseqID, Mouse, Region, Method) 
    DissociationMetaData$title <- as.factor(paste(DissociationMetaData$Mouse, DissociationMetaData$Region, DissociationMetaData$Method,sep=" "))
    DissociationMetaData$Organism <- "Mus musculus"
    DissociationMetaData$molecule <- "RNA"
    DissociationMetaData$file <- "/abuntance.txt"
    DissociationMetaData$alignment <- as.factor(paste(DissociationMetaData$RNAseqID,DissociationMetaData$file, sep=""))
    DissociationMetaData$strain <- "C57BL/6"

    #write.csv(DissociationMetaData, "~/Github/DissociationTest/data/DissociationMetaData.csv", row.names=T)

    ## rename, relevel, slim, and prep for DESEq2 
    rownames(DissociationColData) <- DissociationColData$RNAseqID
    DissociationColData <- DissociationColData %>%
      filter(Mouse %in% c("15-100")) %>% droplevels()
    savecols <- as.character(DissociationColData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    DissociationCountData <- DissociationCountData %>% select(one_of(savecols)) # keep good samples
    DissociationColData <- rename(DissociationColData, c("Method"="Treatment"))
    colnames(DissociationColData)

    ##  [1] "RNAseqID"  "Mouse"     "year"      "Genotype"  "jobnumber"
    ##  [6] "Region"    "Group"     "Conflict"  "APA"       "Treatment"
    ## [11] "dodgy"     "daytime"   "Slice"     "Date"

    DissociationColData <- DissociationColData %>%
      select(RNAseqID,Mouse, Genotype,jobnumber, Region, Treatment, Date)
    write.csv(DissociationColData, "../data/DissociationColData.csv",row.names=F)

    # save files to dissociation project directory
    #write.csv(DissociationColData, "~/Github/DissociationTest/data/DissociationColData.csv", row.names=F)
    #write.csv(DissociationCountData, "~/Github/DissociationTest/data/DissociationCountData.csv", row.names=T)

    fmr1CountData <- read.csv("../data/rnaseq/countbygene.csv", header=T, check.names = F, row.names = 1)
    fmr1ColData <- read.csv("../data/rnaseq/kallistosamples.csv", header=T)

    fmr1ColData <- fmr1ColData %>%
      filter(jobnumber == "JA17009")
    savecols <- as.character(fmr1ColData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    fmr1CountData <- fmr1CountData %>% select(one_of(savecols)) # keep good samples

    # save files to dissociation project directory
    write.csv(fmr1ColData, "../data/rnaseq/fmr1ColData.csv", row.names=F)
    write.csv(fmr1CountData, "../data/rnaseq/fmr1CountData.csv", row.names=T)

    BehaviorSlimCountData <- read.csv("../data/rnaseq/countbygene.csv", header=T, check.names = F, row.names = 1)
    BehaviorSlimColData <- read.csv("../data/rnaseq/kallistosamples.csv", header=T)

    BehaviorSlimColData <- BehaviorSlimColData %>%
      filter(Group %in% c("consistent", "control")) %>%
      filter(Genotype != "FMR1") %>%
      filter(!grepl("16-", Mouse)) %>%
      filter(Conflict != "Conflict")
    savecols <- as.character(BehaviorSlimColData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    BehaviorSlimCountData <- BehaviorSlimCountData %>% select(one_of(savecols)) # keep good samples

    # save files to dissociation project directory
    write.csv(BehaviorSlimColData, "~/Github/DissociationTest/data/BehaviorSlimColData.csv", row.names=F)
    write.csv(BehaviorSlimCountData, "~/Github/DissociationTest/data/BehaviorSlimCountData.csv", row.names=T)

    WT2015countData <- read.csv("../data/rnaseq/countbygene.csv", header=T, check.names = F, row.names = 1)
    WT2015colData <- read.csv("../data/rnaseq/kallistosamples.csv", header=T)

    WT2015colData <- WT2015colData %>%
      filter(!grepl("16-", Mouse)) %>%
      filter(!grepl("100", Mouse)) %>%
      filter(!grepl("homecage", Group))
      
    savecols <- as.character(WT2015colData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    WT2015countData <- WT2015countData %>% select(one_of(savecols)) # keep good samples

    # save files to dissociation project directory
    #write.csv(WT2015colData, "~/Github/IntegrativeProjectWT2015/data/WT2015ColData.csv", row.names=F)
    #write.csv(WT2015countData, "~/Github/IntegrativeProjectWT2015/data/WT2015CountData.csv", row.names=T)

    ## WithHomecage
    IntegrativeWT2015countData <- read.csv("../data/rnaseq/countbygene.csv", header=T, check.names = F, row.names = 1)
    IntegrativeWT2015colData <- read.csv("../data/rnaseq/kallistosamples.csv", header=T)
    IntegrativeWT2015colData <- IntegrativeWT2015colData %>%
      filter(!grepl("16-", Mouse)) %>%
      filter(!grepl("100", Mouse)) 
      
    savecols <- as.character(IntegrativeWT2015colData$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector
    IntegrativeWT2015countData <- IntegrativeWT2015countData %>% select(one_of(savecols)) # keep good samples

    IntegrativeWT2015colData <- IntegrativeWT2015colData %>%
      select(RNAseqID, Mouse, year, Genotype, Region, jobnumber, Group, APA, Conflict)

    IntegrativeWT2015colData$APA_Conflict <- as.factor(paste(IntegrativeWT2015colData$APA, IntegrativeWT2015colData$Conflict, sep="_"))
    levels(IntegrativeWT2015colData$APA_Conflict)

    ## [1] "NA_NA"              "Trained_Conflict"   "Trained_NoConflict"
    ## [4] "Yoked_Conflict"     "Yoked_NoConflict"

    IntegrativeWT2015colData$Project <- ifelse(grepl("NA_NA", IntegrativeWT2015colData$APA_Conflict), "homecage", 
                                        ifelse(grepl("Trained_Conflict", IntegrativeWT2015colData$APA_Conflict), "conflict",
                                        ifelse(grepl("Trained_NoConflict", IntegrativeWT2015colData$APA_Conflict), "trained",
                                        ifelse(grepl("Yoked_NoConflict", IntegrativeWT2015colData$APA_Conflict), "yoked",
                                        ifelse(grepl("Yoked_Conflict", IntegrativeWT2015colData$APA_Conflict), "stressed", "NA")))))
                                                                                                      
    # save files to IntegrativeWT2015 project directory
    write.csv(IntegrativeWT2015colData, "~/Github/IntegrativeProjectWT2015/data/IntegrativeWT2015ColData.csv", row.names=F)
    write.csv(IntegrativeWT2015countData, "~/Github/IntegrativeProjectWT2015/data/IntegrativeWT2015CountData.csv", row.names=T)

    # save files to Dissociation project directory
    write.csv(IntegrativeWT2015colData, "~/Github/DissociationTest/data/StressCognitionTestColData.csv", row.names=F)
    write.csv(IntegrativeWT2015countData, "~/Github/IntegrativeProjectWT2015/data/StressCognitionTestCountData.csv", row.names=T)

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
    ## [1] plyr_1.8.4  dplyr_0.5.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.9     digest_0.6.12   rprojroot_1.2   assertthat_0.1 
    ##  [5] R6_2.2.0        DBI_0.6         backports_1.0.5 magrittr_1.5   
    ##  [9] evaluate_0.10   stringi_1.1.2   lazyeval_0.2.0  rmarkdown_1.3  
    ## [13] tools_3.3.1     stringr_1.2.0   yaml_2.1.14     htmltools_0.3.5
    ## [17] knitr_1.15.1    tibble_1.2
