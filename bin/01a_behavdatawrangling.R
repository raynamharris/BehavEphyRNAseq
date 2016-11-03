# Part 1: Reading and analyzing beahvior and physiology data

## load libraries -----
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!
library(reshape2) #@ for melting dataframe
library(magrittr)

## read and wrangle the data ----

## read the data 
setwd("~/Github/BehavEphyRNAseq/data/behavior/")
behav <- read.csv("Data2013_2016_forAnalysis.csv", header=TRUE, 
                  stringsAsFactors = FALSE, na.strings = c("", "ND", "N/A"))
str(behav)

## drop columns that are not behaviorally relevant
## see apameasurementnotes.doc from maddy
behav <- behav %>%
  select(-p.miss, -TotalTime.s., -filename, -TotalPath.m., -Time1stShock, 
         -Path1stShock, -Speed..cm.s., -SdSpeed, -Linearity,- MaxPathAvoid)
names(behav)

## rename columns 
names(behav)[names(behav)=="TotalPath.Arena."] <- "TotalPathArena.m"
names(behav)[names(behav)=="sd.Speed.Arena."] <- "SdevSpeedArena."
names(behav)[names(behav)=="TotalPath.Arena."] <- "TotalPathArena.m"
names(behav)[names(behav)=="X.Shock"] <- "NumShock"
names(behav)[names(behav)=="X.Entrances"] <- "NumEntrances"
names(behav)[names(behav)=="Speed..Arena."] <- "SpeedArena.cm.s"
names(behav)[names(behav)=="Entr.Dist.1.m."] <- "Dist1stEntr.m."
names(behav)[names(behav)=="Speed..Arena."] <- "SpeedArena.cm.s"
names(behav) # check all good

## remove bad/extra rows
behav <- behav %>% dplyr::filter(Year %in% c("2012", "2013", "2014", "2015", "2016", "2017")) ## remove extra headers
behav <- behav %>% dplyr::filter(!(TrainSession %in% c("C4", "T7"))) ## remove extra sessions

## make characters either factors or numbers, as appropriate
colsfactor = c(1:13)
colsnumeric = c(14:52)
behav[,colsfactor] %<>% lapply(function(x) as.factor(as.character(x)))
behav[,colsnumeric] %<>% lapply(function(x) as.numeric(as.character(x)))
str(behav)
names(behav)
summary(behav)

## some aniamls have a 0 for path to 1st entrance because they were dropped into the schock zone. THose are replaced with NA. 
## all values of -1 are removed because they really mean, couldnt be calculated.

behav$Path1stEntr <-ifelse((behav$Path1stEntr != 0), behav$Path1stEntr, NA)
behav$Speed1stEntr.cm.s. <-ifelse((behav$Speed1stEntr.cm.s. != -1), behav$Speed1stEntr.cm.s., NA)
behav$Speed2ndEntr <-ifelse((behav$Speed2ndEntr != -1), behav$Speed2ndEntr, NA)
behav$Path2ndEntr <-ifelse((behav$Path2ndEntr != -1), behav$Path2ndEntr, NA)
summary(behav)

## recalculate speed 1 and speed 2
behav <- behav %>% mutate(Speed1 = Path1stEntr / Time1stEntr)
behav <- behav %>% mutate(Speed2 = Path2ndEntr / Time2ndEntr)
select(behav, Speed1, Speed2, Speed1stEntr.cm.s., Speed2ndEntr) ## shows definately not the same thign!!!

## create a entra shock diff
behav <- behav %>% mutate(EntrShockDiff = NumShock - NumEntrances)

## rename and revalue some values 
behav$Genotype <- revalue(behav$Genotype, c("FMR1" = "FMR1KO")) 
behav$Genotype <- revalue(behav$Genotype, c("FMR1?" = "FMR1KO")) 
behav$Genotype <- factor(behav$Genotype, levels = c("WT", "FMR1KO"))
levels(behav$Genotype)

behav$TrainGroup <- revalue(behav$TrainGroup, c("control" = "untrained")) 
behav$TrainGroup <- revalue(behav$TrainGroup, c("control?" = "untrained")) 
levels(behav$TrainGroup)

#reaplace NAs in trainSequence with train
behav$TrainSequence[is.na(behav$TrainSequence)] <- "train"

## create combinatorial factor columns 
behav$genoYear <- as.factor(paste(behav$Genotype,behav$Year,sep="_"))
behav$genoYear <- factor(behav$genoYear, 
                    levels = c("WT_2013", "WT_2015", "WT_2016", "FMR1KO_2013","FMR1KO_2016"))
levels(behav$genoYear)

behav$APA <- as.factor(paste(behav$TrainSequence,behav$TrainGroup,sep="_"))
levels(behav$APA)
behav$APA <- revalue(behav$APA, c("train_trained" = "trained")) 
behav$APA <- revalue(behav$APA, c("train_untrained" = "untrained")) 
behav$APA <- revalue(behav$APA, c("train-conflict_trained" = "trained_conflict")) 
behav$APA <- revalue(behav$APA, c("train-conflict_yoked" = "yoked_conflict")) 
behav$APA <- revalue(behav$APA, c("train-train_trained" = "trained_trained")) 
behav$APA <- revalue(behav$APA, c("train-train_yoked" = "yoked_trained")) 
behav$APA <- factor(behav$APA, 
                        levels = c("untrained", "yoked_trained", 
                                   "yoked_conflict", "trained", 
                                   "trained_trained", "trained_conflict"))
levels(behav$APA)

behav$genoAPA <- as.factor(paste(behav$Genotype,behav$APA, sep="_"))
behav$genoAPA <- factor(behav$genoAPA, 
                        levels = c("WT_untrained", "WT_trained", 
                                   "WT_yoked_trained", "WT_trained_trained", 
                                   "WT_yoked_conflict", "WT_trained_conflict", 
                                   "FMR1KO_untrained", "FMR1KO_trained",
                                   "FMR1KO_yoked_trained", "FMR1KO_trained_trained",
                                   "FMR1KO_yoked_conflict", "FMR1KO_trained_conflict"))

behav$TrainSessionCombo <- behav$TrainSession
levels(behav$TrainSessionCombo)
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C1" = "T4_C1")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T4" = "T4_C1")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C2" = "T5_C2")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T5" = "T5_C2")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C3" = "T6_C3")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T6" = "T6_C3")) 
behav$TrainSessionCombo <- factor(behav$TrainSessionCombo, ## set levels
                        levels = c("Hab", "T1","T2","T3","Retest",
                                   "T4_C1","T5_C2","T6_C3","Retention"))
levels(behav$TrainSessionCombo)

behav$genoAPATrainSessionCombo <- as.factor(paste(behav$genoAPA, behav$TrainSessionCombo, sep="_"))
behav$genoAPATrainSessionComboInd <- as.factor(paste(behav$genoAPAsessionCombo, behav$ID, sep="_"))

behav$pair1 <- as.factor(paste(behav$ID,behav$TrainSessionComboDay, sep="_"))
behav$pair2 <- as.factor(paste(behav$PairedPartner,behav$TrainSessionComboDay, sep="_"))

## reorders dataframe 
names(behav)
behav <- behav[c(2,1,3:13,56:63,14:55)]  
names(behav)


## subset the data -----
# by experimenter
#maddy <- behav %>% filter(Experimenter == "Maddy") 
#jma <- behav %>% filter(Experimenter != "Maddy") 

# by genotype
#wt <- behav %>% filter(Genotype == "WT")
#frmr1 <- behav %>% filter(Genotype != "WT")  

# by year
#y2015 <- behav %>%  filter(Year == "2015")

# by experiement by genotype experimenter training
#maddyWT <- behav %>%  filter(Experimenter == "Maddy", Genotype == "WT")
#maddyWTtrained <- behav %>%  filter(Experimenter == "Maddy", Genotype == "WT",  TrainGroup == "trained") 


## Create novel dataframes
## make a df to look at number of shock actually received by the yoked animals
names(behav)
yoked <- behav %>% 
  filter(Experimenter == "Maddy") %>%
  filter(grepl("yoked", APA) ) %>%
  select(pair1, pair2, ID, Genotype, APA, Year, TrainGroup, TrainSequence,TrainSessionCombo, NumShock, NumEntrances, TimeTarget, Time1stEntr, Speed1, Speed2, EntrShockDiff) %>% 
  droplevels()
trainedpair <- behav %>% 
  filter(Experimenter == "Maddy") %>%
  select(pair1, pair2, ID, Genotype, APA, Year, TrainGroup, TrainSequence, TrainSessionCombo, NumShock, NumEntrances, TimeTarget, Time1stEntr, Speed1, Speed2, EntrShockDiff) %>% 
  droplevels()
## rename columns 
names(yoked)[1] <- "yoked"
names(yoked)[2] <- "trained"
names(trainedpair)[1] <- "trained"
names(trainedpair)[2] <- "yoked"
yokedtrainedpair <- left_join(yoked, trainedpair, by = "yoked") ## join an caluclate some values
rm(yoked)  # only need to keep yokedtrainedpair
rm(trainedpair) # only need to keep yokedtrainedpair

#write.csv(yokedtrainedpair, "yokedtrainedpair.csv", row.names = FALSE)


### Making the data long ----
names(behav)
behavbysession <- melt(behav, id = c(1:21))

behavbysession <- filter(behavbysession, !grepl("TotalTime.s|p.miss", variable )) %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>%  droplevels() 
## create the bysession column and wide
behavbysession$bysession <- as.factor(paste(behavbysession$TrainSessionCombo, behavbysession$variable, sep="_"))
behavbysession <- dcast(behavbysession, ID + APA + Genotype + TrainSequence + TrainGroup ~ bysession, value.var= "value", fun.aggregate = mean)
summary(behavbysession) 
head(behavbysession)


#write.csv(behav, "behav.csv", row.names = FALSE)

