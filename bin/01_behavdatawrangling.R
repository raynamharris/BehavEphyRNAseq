# Part 1: Reading and analyzing beahvior and physiology data

## load libraries -----
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!
library(reshape2) #@ for melting dataframe


## read and wrangle the data ----

## read the data 
setwd("~/Github/BehavEphyRNAseq/data/sample_info/")
behav <- read.csv("APA_2013-2016.csv", header=TRUE, stringsAsFactors = FALSE, na.strings = c("", "ND", "N/A"))
str(behav)

behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C1" = "T4_C1")) 


## rename columns 
names(behav)[18] <- "NumEntrances" # previously  X.Entrances 
names(behav)[23] <- "NumShock"    # previously X.Shock 
names(behav)[26] <- "Speed.cm.s" # Speed..cm.s.
names(behav) # check all good
head(behav)

## remove bad/extra rows
behav <- behav %>% dplyr::filter(Year %in% c("2012", "2013", "2014", "2015", "2016", "2017")) ## remove extra headers
behav <- behav %>% dplyr::filter(TrainSession != "C4?") ## remove ?C in train train animal
behav <- behav %>% dplyr::filter(!is.na(TotalTime.s.)) ## remove file with no data
behav <- behav %>% dplyr::filter(!grepl("16-357A|16-357B|16-357D", ID)) ## 3 bad yoked animals

## make characters either factors or numbers, as appropriate
cols = c(1:4,6:14)
cols2 = c(15:58)
behav[,cols] %<>% lapply(function(x) as.factor(as.character(x)))
behav[,cols2] %<>% lapply(function(x) as.numeric(as.character(x)))
str(behav)
names(behav)

## rename and revalue some values 
behav$Genotype <- revalue(behav$Genotype, c("FMR1" = "FMR1-KO")) 
behav$Genotype <- revalue(behav$Genotype, c("FMR1?WT?" = "FMR1-KO")) 
behav$Genotype <- factor(behav$Genotype, levels = c("WT", "FMR1-KO"))

behav$TrainGroup <- revalue(behav$TrainGroup, c("control" = "untrained")) 
behav$TrainSequence[is.na(behav$TrainSequence)] <- "untrained" ## make NA more meaningful
behav$TrainSequence <- as.factor(behav$TrainSequence)

## create combinatorial factor columns 
behav$genoYear <- as.factor(paste(behav$Genotype,behav$Year,sep="_"))
behav$genoYear <- factor(behav$genoYear, 
                    levels = c("WT_2013", "FMR1-KO_2014", "WT_2014", "WT_2015", 
                               "FMR1-KO_2016", "WT_2016"))

behav$genoYear <- factor(behav$genoYear, 
                         levels = c("WT_2013", "WT_2014", "WT_2015", "WT_2016",
                                    "FMR1-KO_2014", "FMR1-KO_2016"))


behav$APA <- as.factor(paste(behav$TrainSequence,behav$TrainGroup,sep="_"))
behav$APA <- revalue(behav$APA, c("train_trained" = "trained")) 
behav$APA <- revalue(behav$APA, c("untrained_untrained" = "untrained")) 
behav$APA <- revalue(behav$APA, c("train-conflict_trained" = "trained_conflict")) 
behav$APA <- revalue(behav$APA, c("train-conflict_yoked" = "yoked_conflict")) 
behav$APA <- revalue(behav$APA, c("train-train_trained" = "trained_trained")) 
behav$APA <- revalue(behav$APA, c("train-train_yoked" = "yoked_trained")) 
behav$APA <- factor(behav$APA, 
                        levels = c("untrained", "yoked_trained", 
                                   "yoked_conflict", "trained", 
                                   "trained_trained", "trained_conflict"))

behav$genoAPA <- as.factor(paste(behav$Genotype,behav$APA, sep="_"))
behav$genoAPA <- factor(behav$genoAPA, 
                        levels = c("WT_untrained", "WT_trained", 
                                   "WT_yoked_trained", "WT_trained_trained", 
                                   "WT_yoked_conflict", "WT_trained_conflict", 
                                   "FMR1-KO_untrained", "FMR1-KO_trained",
                                   "FMR1-KO_yoked_trained", "FMR1-KO_trained_trained",
                                   "FMR1-KO_yoked_conflict", "FMR1-KO_trained_conflict"))

behav$genoAPAyear <- as.factor(paste(behav$genoAPA,behav$Year, sep="_"))
behav$genoAPAsession  <- as.factor(paste(behav$genoAPA,behav$TrainSession, sep="_"))
behav$genoAPAsessionDay  <- as.factor(paste(behav$genoAPAsession,behav$Day, sep="_"))
behav$genoAPAsessionDayInd <- as.factor(paste(behav$genoAPAsessionDay,behav$ID, sep="_"))

behav$TrainSessionCombo <- behav$TrainSession
levels(behav$TrainSessionCombo)
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C1" = "T4_C1")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T4" = "T4_C1")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C2" = "T5_C2")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T5" = "T5_C2")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C3" = "T6_C3")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T6" = "T6_C3")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T7" = "T+_C+")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T8" = "T+_C+"))
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C4" = "T+_C+"))
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C5" = "T+_C+")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C6" = "T+_C+")) 
behav$TrainSessionCombo <- factor(behav$TrainSessionCombo, ## set levels
                        levels = c("Hab", "T1","T2","T3","Retest",
                                   "T4_C1","T5_C2","T6_C3","T+_C+","Retention"))
behav <- behav %>% dplyr::filter(TrainSessionCombo != "T+_C+")  ## remove T+_C+"

behav$TrainSessionComboDay <- as.factor(paste(behav$TrainSessionCombo, behav$Day, sep="_"))
behav$genoAPAsessionCombo <- as.factor(paste(behav$genoAPA, behav$TrainSessionCombo, sep="_"))
behav$genoAPAsessionComboInd <- as.factor(paste(behav$genoAPAsessionCombo, behav$ID, sep="_"))

behav$pair1 <- as.factor(paste(behav$ID,behav$TrainSessionComboDay, sep="_"))
behav$pair2 <- as.factor(paste(behav$PairedPartner,behav$TrainSessionComboDay, sep="_"))

## reorders dataframe 
names(behav)
behav <- behav[c(2,59:71,3:9,1,10:58)]  
names(behav)


## subset the data -----
# by experimenter
maddy <- behav %>% filter(Experimenter == "Maddy") 
jma <- behav %>% filter(Experimenter != "Maddy") 

# by genotype
wt <- behav %>% filter(Genotype == "WT")
frmr1 <- behav %>% filter(Genotype != "WT")  

# by year
y2015 <- behav %>%  filter(Year == "2015")

# by experiement by genotype experimenter training
maddyWT <- behav %>%  filter(Experimenter == "Maddy", Genotype == "WT")
maddyWTtrained <- behav %>%  filter(Experimenter == "Maddy", Genotype == "WT",  TrainGroup == "trained") 

## Create novel dataframes

## make a df to look at number of shock actually received by the yoked animals
names(behav)
yoked <- behav %>% 
  filter(Experimenter == "Maddy") %>%
  filter(grepl("yoked", APA) ) %>%
  select(pair1, pair2, ID, Genotype, APA, Year, TrainGroup, TrainSequence,TrainSessionCombo, NumShock, NumEntrances, TimeTarget, Time1stEntr, Time1stShock) %>% 
  droplevels()
trainedpair <- behav %>% 
  filter(Experimenter == "Maddy") %>%
  select(pair1, pair2, ID, Genotype, APA, Year, TrainGroup, TrainSequence, TrainSessionCombo, NumShock, NumEntrances, TimeTarget, Time1stEntr, Time1stShock) %>% 
  droplevels()
## rename columns 
names(yoked)[1] <- "yoked"
names(yoked)[2] <- "trained"
names(trainedpair)[1] <- "trained"
names(trainedpair)[2] <- "yoked"
## join an caluclate some values
yokedtrainedpair <- left_join(yoked, trainedpair, by = "yoked")
yokedtrainedpair <- yokedtrainedpair %>%  
  mutate(NumShockXYequal = ifelse(NumShock.x == NumShock.y, "equal", "not equal")) %>%  
  arrange(NumShock.x) %>% 
  select(Genotype.x,Year.x,TrainSequence.x,TrainSessionCombo.x,
         ID.x, TrainGroup.x, NumShock.x, NumEntrances.x, TimeTarget.x, Time1stEntr.x, Time1stShock.x,
         ID.y, TrainGroup.y, NumShock.y, NumEntrances.y, TimeTarget.y, Time1stEntr.y, Time1stShock.y)

names(yokedtrainedpair)

#write.csv(yokedtrainedpair, "yokedtrainedpair.csv", row.names = FALSE)


### For heat map melt to make long  ## goes with a heatmap in next script ----
behavbysession <- melt(behav, id=c("ID","APA","genoAPA","genoAPAsession","genoAPAsessionDay", "genoYear", "genoAPAsessionCombo",
                               "genoAPAsessionDayInd","TrainSessionCombo" ,"Genotype", "TrainSessionComboDay", "genoAPAyear",
                               "genoAPAsessionComboInd","TrainProtocol","TrainSequence","TrainGroup","Day","TrainSession",
                               "ShockOnOff","Year","PairedPartner","Experimenter",
                               "Housing","TestLocation", "filename", "pair1", "pair2"))
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

