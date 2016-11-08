
library(ggdendro) ## for dendrograms!!
library(magrittr) ## to use the weird pipe
library(gplots) ##for making awesome plots
library(cowplot) ## for some easy to use themes
library(tidyr)
library(dplyr) ## for filtering and selecting rows
library(reshape2) #@ for melting dataframe
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!


setwd("~/Github/BehavEphyRNAseq/data/sample_info//")
Traits <- read.csv("JA16444samples.csv", sep=",", header = TRUE, stringsAsFactors=FALSE, na.string = "NA")
rownames(Traits) <- Traits$RNAseqID    # set $genoAPAsessionInd as rownames
names(Traits)

#keeping informative volumns
Traits <- Traits[c(1,3,5:6,10:11)] 
tail(Traits)
str(Traits)

## remove 100 and 100 animals because they aren't APA trained
## remove mice 147D_CA1_1 and 145B_CA3_1 because they produced almost 0 zeros
## remove 147 and 148 because they are home cage animals
Traits <- Traits %>% dplyr::filter(!grepl("100-|101-|147D-CA1-1|145B-CA3-1|147-|148-", RNAseqID)) 


## adding combinatorial traits
Traits$APAconflict <- as.factor(paste(Traits$APA, Traits$Conflict, sep="_")) 
summary(Traits)

#make a column that thas id without the dash
Traits$ID <- gsub("[[:punct:]]", "", Traits$Mouse)


## Subset by punch or training
TraitsCA1 <- Traits %>% filter(Punch == "CA1") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict, ID)
TraitsCA3 <- Traits %>% filter(Punch == "CA3") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict, ID)
TraitsDG <- Traits %>% filter(Punch == "DG") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict, ID)
TraitsCA1DG <- Traits %>% filter(Punch %in% c("DG", "CA1")) %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict, ID, Punch)
TraitsTrained <- Traits %>% filter(APA == "Trained") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict, ID)

## pick which data to subset
Traits <- TraitsCA1DG

## making it a numeric
Traits$Mouse <- as.integer(factor(Traits$Mouse))
Traits$Conflict <- as.integer(factor(Traits$Conflict))
Traits$APA <- as.integer(factor(Traits$APA))
Traits$Slice <- as.integer(factor(Traits$Slice))
Traits$APAconflict <- as.integer(factor(Traits$APAconflict))
Traits$ID <- as.integer(factor(Traits$ID))
Traits$Punch <- as.integer(factor(Traits$Punch))

Traits$Mouse <- as.numeric(factor(Traits$Mouse))
Traits$Conflict <- as.numeric(factor(Traits$Conflict))
Traits$APA <- as.numeric(factor(Traits$APA))
Traits$Slice <- as.numeric(factor(Traits$Slice))
Traits$APAconflict <- as.numeric(factor(Traits$APAconflict))
Traits$ID <- as.numeric(factor(Traits$ID))
Traits$Punch <- as.integer(factor(Traits$Punch))

head(Traits)
str(Traits)
summary(Traits)

## make gene the row name then round all value to nearest 1s place
row.names(Traits) <- Traits$RNAseqID
Traits[1] <- NULL





