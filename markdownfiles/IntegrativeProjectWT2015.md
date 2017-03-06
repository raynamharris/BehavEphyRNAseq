#### Prepare the enviornment by loading the necessary packages

    library("tidyr")
    library("dplyr")
    library("plyr")
    library("reshape2")

#### wrangle data

    # all animals that have been process in some shape or form
    animals <- read.csv("../data/sample_info/animals.csv", header=TRUE, sep="," )
    # information on  the tissue samples collected 
    punches <- read.csv("../data/sample_info/punches.csv", header=TRUE, sep="," )
    # merging the two by the Mouse name
    full <- join(animals, punches, by = "Mouse", type = "full", match = "all")
    # filter to just the 2015 WT samples that have been sequenced
    WT2015samples <- full %>%
      filter(jobnumber %in% c("JA16268", "JA16444")) %>%
      distinct(RNAseqID, Tube, Mouse, Genotype, Conflict, APA, Group,
               Behavior, E.phy, Punch, Slice, Date, jobnumber) %>%
      droplevels()
