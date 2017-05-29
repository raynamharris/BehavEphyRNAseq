    library("dplyr") ## for filtering and selecting rows
    library("plyr")  ## for renmaing factors
    library("reshape2") ##  for melting dataframe
    knitr::opts_chunk$set(fig.path = '../results/02b_KallistoGather/') # set output file for figures 

RNAseq samples
--------------

In the summers of 2015 and 2016, I processed a bunch of hippocampal
tissue samples from 59 mice. Most mice were trained in an active place
avoidance task or used as yoked controls; however, a few animals were
taken straight from the home cage. This data has been cleaned using a
different script `02a_tidysamples.R`. After gathering the kallisto data,
I will do one more cleaning step at the end.

    tidysamples <- read.csv("../data/rnaseq/tidysamples.csv", header=TRUE)
    rownames(tidysamples) <- tidysamples$RNAseqID   # make RNAseqID the row name
    tidysamples$year <- as.factor(tidysamples$year) # was being read as intger
    tidysamples$Slice <- as.factor(tidysamples$Slice) # was being read as intger
    tidysamples <- tidysamples %>% arrange(RNAseqID) # sort by RNAseqID

Kallisto Gather
---------------

The kallisto output gives you read counts for sample in an abundance
file for every single sample. This portion of the code goes through and
finds each samples' abundance.tsv file, extracts the data, and combines
it all into a dataframe. The `counts` file is unnormalized, but the
`tpm` is the data after being normalized by transcripts per million.  
(P.S. Unfortunately, I have no idea how to do this next part without
changing directories. This script was developed with assistance from
Anna Batthenhouse and Dennis Whylie.)

    setwd("../data/rnaseq/04_kallistoquant/")
    ## this will create lists of all the samples
    kallistoDirs = dir(".")
    kallistoDirs = kallistoDirs[!grepl("\\.(R|py|pl|sh|xlsx?|txt|tsv|csv|org|md|obo|png|jpg|pdf)$",
    kallistoDirs, ignore.case=TRUE)]

    kallistoFiles = paste0(kallistoDirs, "/abundance.tsv")
    names(kallistoFiles) = kallistoDirs
    if(file.exists(kallistoFiles))
      kallistoData = lapply(
      kallistoFiles,
      read.table,
      sep = "\t",
      row.names = 1,
      header = TRUE
    )

    ## this for loop uses the reduce function to make two data frame with counts or tpm from all the samples
    ids = Reduce(f=union, x=lapply(kallistoData, rownames))
    if (all(sapply(kallistoData, function(x) {all(rownames(x)==ids)}))) {
      count = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$est_counts}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
      tpm = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$tpm}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    }

Tranform transcipts-level measurements to gene-level measurements
-----------------------------------------------------------------

Kallisto is cool because it does give you transcript level stuff, but
right now I think I have more power if I just look at gene level stuff.
I'll go back and look at transcripts if I want to. I use some string
splitting to take the very long transcript identifying and create a
`geneids` file that has all the database identifiers for each
transcript.

    ## make a dataframe with the parts of the gene id as columns
    geneids <- count[c(1)] 
    geneids$ENSMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 1)
    geneids$ENSMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 2)
    geneids$OTTMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 3)
    geneids$OTTMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 4)
    geneids$transcript <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 5)
    geneids$gene <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 6)
    geneids$length <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 7)
    geneids$structure1 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 8)
    geneids$structure2 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 9)
    geneids$structure3 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 10)
    geneids$transcript_lenght <- as.factor(paste(geneids$transcript, geneids$length, sep="_"))

    # tpm to tpmbygene
    tpmbygene <-  full_join(geneids, tpm) # merge tpm and genids

    ## Joining, by = "id"

    tpmbygene <- tpmbygene[-c(1:6,8:12)]   ## keep gene name and tpm for samples)
    tpmbygene <- melt(tpmbygene, id=c("gene")) ## lenghten 
    tpmbygene <- dcast(tpmbygene, gene ~ variable, value.var= "value", fun.aggregate=mean) #then widen by sum
    row.names(tpmbygene) <- tpmbygene$gene ## make gene the row name
    tpmbygene[1] <- NULL ## make gene the row name
    tpmbygene <- round(tpmbygene) #round all value to nearest 1s place


    # count to countbygene
    countbygene <- full_join(geneids, count) # merge count and genids

    ## Joining, by = "id"

    countbygene <- countbygene[-c(1:6,8:12)]   ## rkeep gene name and counts for samples)
    countbygene <- melt(countbygene, id=c("gene")) ## lenghten 
    countbygene  <- dcast(countbygene, gene ~ variable, value.var= "value", fun.aggregate=sum) #then widen by sum
    row.names(countbygene) <- countbygene$gene ## make gene the row name
    countbygene[1] <- NULL ## make gene the row name
    countbygene <- round(countbygene) #round all value to nearest 1s place

Final samples removal
---------------------

At this point, we need to remove a few more samples that either didn't
have an abundance tsv file or were deemed outliers.

    ## remove samples that never produced quality reads
    savesamples <- colnames(countbygene)
    kallistosamples <- tidysamples %>%
      filter(RNAseqID %in% savesamples) %>% droplevels()

    ## remove outlier samples
    kallistosamples <- kallistosamples %>%
      filter(RNAseqID != "146C-CA3-4", RNAseqID != "16-116D") # two outliters
    savecols <- as.character(kallistosamples$RNAseqID) #selects all good samples
    savecols <- as.vector(savecols) # make it a vector

    countbygene <- countbygene %>% select(one_of(savecols)) # keep good samples
    count <- count %>% select(one_of(savecols)) # keep good samples
    tpm <- tpm %>% select(one_of(savecols)) # keep good samples
    tpmbygene <- tpmbygene %>% select(one_of(savecols)) # keep good samples

    summary(kallistosamples)

    ##       RNAseqID      Mouse      year    Genotype    jobnumber  Region  
    ##  100-CA1-1: 1   15-100 :14   2015:68   FMR1: 8   JA16268: 4   CA1:41  
    ##  100-CA1-2: 1   15-143C: 3   2016:16   WT  :76   JA16444:64   CA3:19  
    ##  100-CA1-3: 1   15-144A: 3                       JA17009:16   DG :24  
    ##  100-CA3-1: 1   15-144C: 3                                            
    ##  100-CA3-4: 1   15-145A: 3                                            
    ##  100-DG-2 : 1   15-146A: 3                                            
    ##  (Other)  :78   (Other):55                                            
    ##         Group          Conflict       APA             Method  
    ##  conflict  :14   Conflict  :34   Trained:27   control    : 7  
    ##  consistent:13   NoConflict:30   Yoked  :37   dissociated: 7  
    ##  control   :37   NA's      :20   NA's   :20   homogenized:70  
    ##  homecage  :20                                                
    ##                                                               
    ##                                                               
    ##                                                               
    ##      dodgy          daytime   Slice       Date   
    ##  allgood:73   afternoon : 2   1:20   9/28/15:14  
    ##  ephys  : 2   beforenoon:14   2:31   7/18/15: 8  
    ##  slice  : 9   earlyAM   : 6   3:21   7/23/15: 6  
    ##               evening   : 5   4:12   7/30/15: 6  
    ##               nighttime : 3          7/10/15: 5  
    ##               norecord  :54          7/15/15: 5  
    ##                                      (Other):40

Write files
-----------

Now, let's save the output files for future use

    write.csv(countbygene, "../data/rnaseq/countbygene.csv", row.names=T)
    write.csv(tpmbygene, "../data/rnaseq/tpmbygene.csv", row.names=T)
    write.csv(geneids, "../data/rnaseq/geneids.csv", row.names=FALSE)
    write.csv(kallistosamples, "../data/rnaseq/kallistosamples.csv", row.names=FALSE)

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
    ## [1] reshape2_1.4.2 plyr_1.8.4     dplyr_0.5.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.9     digest_0.6.12   rprojroot_1.2   assertthat_0.1 
    ##  [5] R6_2.2.0        DBI_0.6         backports_1.0.5 magrittr_1.5   
    ##  [9] evaluate_0.10   stringi_1.1.2   lazyeval_0.2.0  rmarkdown_1.3  
    ## [13] tools_3.3.1     stringr_1.2.0   yaml_2.1.14     htmltools_0.3.5
    ## [17] knitr_1.15.1    tibble_1.2
