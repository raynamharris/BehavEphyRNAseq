# Part 1: Reading and analyzing beahvior and physiology data

## load libraries -----
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!
library(reshape2) #@ for melting dataframe
library(ggdendro) ## for dendrograms!!
library(magrittr) ## to use the weird pipe

## wrangle the raw wt and fmr1 dataframes ----

## read the data 
setwd("~/Github/BehavEphyRNAseq/data/sample_info/")
behav <- read.csv("APA_2013-2016.csv", header=TRUE, stringsAsFactors = FALSE, na.strings = c("", "ND", "N/A"))
str(behav)

## rename columns -----
names(behav)[18] <- "NumEntrances"
names(behav)[23] <- "NumShock"
names(behav)[26] <- "Speed_cm_s"
names(behav) # check all good
head(behav)

## remove extra header rows ----
behav <- behav %>% dplyr::filter(Year %in% c("2012", "2013", "2014", "2016", "2017"))

## make strings factors ----
cols = c(1:4,6:14)
cols2 = c(15:58)
behav[,cols] %<>% lapply(function(x) as.factor(as.character(x)))
behav[,cols2] %<>% lapply(function(x) as.numeric(as.character(x)))
str(behav)
names(behav)

## rename some values -----
behav$Genotype <- revalue(behav$Genotype, c("FMR1" = "FMR1-KO")) 
behav$Genotype <- revalue(behav$Genotype, c("FMR1?WT?" = "FMR1-KO")) 
behav$TrainGroup <- revalue(behav$TrainGroup, c("control" = "untrained")) 
behav$TrainSequence[is.na(behav$TrainSequence)] <- "untrained" ## make NA more meaningful
behav$TrainSequence <- as.factor(behav$TrainSequence)

## create columns for genotype*APA, genotype*APA*session, and genotype*APA*session*IND
behav$genoAPA <- as.factor(paste(behav$Genotype,behav$TrainSequence, sep="_"))
behav$genoAPAsession  <- as.factor(paste(behav$genoAPA,behav$TrainSession, sep="_"))
behav$genoAPAsessionInd <- as.factor(paste(behav$genoAPAsession,behav$ID, sep="_"))
head(behav)


## reorders dataframe ----
behav <- behav[c(2,59:61,3:13,1,14:58)]  
names(behav)

### melt the wtfmr1 df to make long for graphics
behav_long <- melt(behav, id=c("ID","genoAPA", "genoAPAsession", "genoAPAsessionInd", 
                               "Genotype", "TrainProtocol", "TrainSequence","TrainGroup",
                               "Day", "TrainSession", "ShockOnOff", "PairedPartner",
                               "Experimenter", "Housing", "TestLocation", "Year", "filename"))
str(behav_long)
head(behav_long)

### Beahvior ggplots!!! -----
## create the color palette
FentonPalette <- c('black','grey50','red','darkorange')
WTPalette <- c('black','red')
FMR1Palette <- c('grey50','darkorange')

## basic format to behavior of groups by session stat smooth

## this plot is no longer legit because of the retention/retest problem
behav %>%
  ggplot(aes(as.numeric(x=TrainSession), y=pTimeTarget, color=genoAPA)) + 
  stat_smooth() + theme_bw() 

behav$Experimenter
behav$TrainSession
behav %>% filter(Experimenter == "Maddy") %>%
  droplevels() %>%
  ggplot(aes(x=TrainSession, y=pTimeTarget, color=genoAPA)) + 
  geom_boxplot() + theme_bw() 
