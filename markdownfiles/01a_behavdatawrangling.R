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
  dplyr::select(-p.miss, -TotalTime.s., -filename, -TotalPath.m., -Time1stShock, 
         -Path1stShock, -Speed..cm.s., -SdSpeed, -Linearity,- MaxPathAvoid)
names(behav)

## rename columns 
names(behav)[names(behav)=="TotalPath.Arena."] <- "TotalPathArena.m"
names(behav)[names(behav)=="sd.Speed.Arena."] <- "SdevSpeedArena"
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


## recalculate speed 1 and speed 2
behav <- behav %>% mutate(Speed1 = Path1stEntr / Time1stEntr)
behav <- behav %>% mutate(Speed2 = Path2ndEntr / Time2ndEntr)

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
behav$APA <- revalue(behav$APA, c("train-conflict_trained" = "Conflict")) 
behav$APA <- revalue(behav$APA, c("train-conflict_yoked" = "Yoked")) 
behav$APA <- revalue(behav$APA, c("train-train_trained" = "Same")) 
behav$APA <- revalue(behav$APA, c("train-train_yoked" = "Yoked")) 
behav$APA <- factor(behav$APA, 
                        levels = c("untrained", "trained", 
                                   "Yoked", "Same", "Conflict"))
levels(behav$APA)

behav$genoAPA <- as.factor(paste(behav$Genotype,behav$APA, sep="_"))
levels(behav$genoAPA)
behav$genoAPA <- factor(behav$genoAPA, 
                        levels = c("WT_untrained", "WT_trained", 
                                   "WT_Yoked", "WT_Same", "WT_Conflict",
                                   "FMR1KO_untrained", "FMR1KO_trained",
                                   "FMR1KO_Yoked", "FMR1KO_Same", "FMR1KO_Conflict"))
levels(behav$genoAPA)


behav$APA2 <- as.factor(paste(behav$TrainSequence,behav$TrainGroup,sep="_"))
behav$APA2 <- revalue(behav$APA2, c("train_trained" = "trained")) 
behav$APA2 <- revalue(behav$APA2, c("train_untrained" = "untrained")) 
behav$APA2 <- revalue(behav$APA2, c("train-conflict_trained" = "Conflict")) 
behav$APA2 <- revalue(behav$APA2, c("train-conflict_yoked" = "YokedConflict")) 
behav$APA2 <- revalue(behav$APA2, c("train-train_trained" = "Same")) 
behav$APA2 <- revalue(behav$APA2, c("train-train_yoked" = "YokedSame")) 
behav$APA2 <- factor(behav$APA2, 
                    levels = c("untrained", "trained", 
                               "YokedSame",  "YokedConflict", "Same", "Conflict"))
levels(behav$APA2)


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
behav$genoAPATrainSessionComboInd <- as.factor(paste(behav$genoAPATrainSessionCombo, behav$ID, sep="_"))

behav$pair1 <- as.factor(paste(behav$ID,behav$TrainSessionCombo, sep="_"))
behav$pair2 <- as.factor(paste(behav$PairedPartner,behav$TrainSessionCombo, sep="_"))

## reorders dataframe 
names(behav)
behav <- behav[c(2,1,3:13,56:64,16:55)]  
names(behav)


## Making the data long then wide----
behavbysession <- melt(behav, id = c(1:22))
behavbysession$bysession <- as.factor(paste(behavbysession$TrainSessionCombo, behavbysession$variable, sep="_"))
behavbysession <- dcast(behavbysession, ID + APA + APA2 + Genotype + TrainProtocol + TrainSequence + TrainGroup + 
                          PairedPartner + Experimenter + Housing + TestLocation + genoYear +  
                          genoAPA + pair1 + pair2 ~ bysession, value.var= "value", fun.aggregate = mean)
summary(behavbysession) 
tail(behavbysession)



## subset the data -----
# by experimenter
#maddy <- behav %>% filter(Experimenter == "Maddy")  
#jma <- behav %>% filter(Experimenter != "Maddy") 


# by year
y2015 <- behav %>%  filter(Year == "2015")
y2015 <- y2015 %>%
  dplyr::select(-genoYear, -genoAPA, -genoAPATrainSessionCombo, 
         -genoAPATrainSessionComboInd, -EntrShockDiff)
names(y2015)

y2015$TrainSessionComboNum <- ifelse(grepl("Hab", y2015$TrainSessionCombo), "1", 
                                        ifelse(grepl("T1", y2015$TrainSessionCombo), "2",
                                               ifelse(grepl("T2", y2015$TrainSessionCombo), "3",
                                                      ifelse(grepl("T3", y2015$TrainSessionCombo), "4", 
                                                             ifelse(grepl("Retest", y2015$TrainSessionCombo), "5",
                                                                    ifelse(grepl("T4_C1", y2015$TrainSessionCombo), "6",
                                                                           ifelse(grepl("T5_C2", y2015$TrainSessionCombo), "7",
                                                                                  ifelse(grepl("T6_C3", y2015$TrainSessionCombo), "8",
                                                                                         ifelse(grepl("Retention", y2015$TrainSessionCombo), "9",NA)))))))))
y2015$TrainSessionComboNum <- as.numeric(as.character(y2015$TrainSessionComboNum))
summary(y2015$TrainSessionComboNum) 
head(y2015)
names(y2015)
y2015 <- y2015[c(1:18,58,19:57)]  
names(y2015)

#write.csv(y2015, '/Users/raynamharris/Github/IntegrativeProjectWT2015/data/01_behaviordata.csv', row.names = F)


# by experiement by genotype experimenter training
#maddyWT <- behav %>%  filter(Experimenter == "Maddy", Genotype == "WT")
#maddyWTtrained <- behav %>%  filter(Experimenter == "Maddy", Genotype == "WT",  TrainGroup == "trained") 


## Actual shocks to yoked animals
## make a df to look at number of shock actually received by the yoked animals
yoked <- behav %>% 
  filter(Experimenter == "Maddy") %>%
  filter(grepl("Yoked", APA)) %>%
  select(pair1, pair2, ID, Genotype, APA, Year, TrainGroup, TrainSequence,TrainSessionCombo, NumShock, NumEntrances, TimeTarget, Time1stEntr, Speed1, Speed2, EntrShockDiff) %>% 
  droplevels()
trainedpair <- behav %>% 
  filter(Experimenter == "Maddy") %>%
  select(pair1, pair2, ID, Genotype, APA, Year, TrainGroup, TrainSequence, TrainSessionCombo, NumShock, NumEntrances, TimeTarget, Time1stEntr, Speed1, Speed2, EntrShockDiff) %>% 
  droplevels()
# rename columns 
names(yoked)[1] <- "yoked"
names(yoked)[2] <- "trained"
names(trainedpair)[1] <- "trained"
names(trainedpair)[2] <- "yoked"
yokedtrainedpair <- left_join(yoked, trainedpair, by = "yoked") ## join an caluclate some values
rm(yoked)  # only need to keep yokedtrainedpair
rm(trainedpair) # only need to keep yokedtrainedpair

## now subset to look at values for the stress test papers
unavoidableshock <- yokedtrainedpair %>%
  select(ID.x, APA.x, TrainSequence.x, TrainSessionCombo.x, NumShock.y) %>%
  filter(TrainSessionCombo.x != "Hab", TrainSessionCombo.x != "Retest", TrainSessionCombo.x != "Retention") %>%
  filter(TrainSequence.x == "train-conflict") %>%
  filter(grepl("143|145|146|148B", ID.x))
str(unavoidableshock)
#write.csv(unavoidableshock, "../../DissociationTest/data/unavoidableshock.csv", row.names = FALSE)



## now subset to look at values for the dissociation test papers
avoidableshock <- yokedtrainedpair %>%
  select(ID.x, APA.x, ID.y, APA.y, TrainSequence.y, TrainSessionCombo.y, TimeTarget.x, TimeTarget.y, NumShock.y) %>%
  filter( TrainSessionCombo.y != "Retest", TrainSessionCombo.y != "Retention") %>%
  filter(TrainSequence.y == "train-train") %>%
  filter(grepl("143|144|146|147", ID.y))
str(avoidableshock)
unavoidable <- avoidableshock %>% 
  select(ID.x, APA.x, TrainSequence.y, TrainSessionCombo.y, TimeTarget.x, NumShock.y)
avoidable <- avoidableshock %>% 
  select(ID.y, APA.y, TrainSequence.y, TrainSessionCombo.y, TimeTarget.y, NumShock.y)
unavoidable <- rename(unavoidable, c("ID.x"="ID"))
avoidable <- rename(avoidable, c("ID.y"="ID"))
unavoidable <- rename(unavoidable, c("APA.x"="Treatment"))
avoidable <- rename(avoidable, c("APA.y"="Treatment"))
unavoidable <- rename(unavoidable, c("TimeTarget.x"="TimeTarget"))
avoidable <- rename(avoidable, c("TimeTarget.y"="TimeTarget"))
avoidableshock <- rbind(unavoidable, avoidable)
avoidableshock$Treatment <- plyr::revalue(avoidableshock$Treatment, c("Yoked"="unavoidable", "Same"="avoidable"))
avoidableshock$pTimeTarget <- avoidableshock$TimeTarget / 600

#write.csv(avoidableshock, "../../DissociationTest/data/avoidableshock.csv", row.names = FALSE)



# by genotype
# wt <- behav %>% filter(Genotype == "WT")
# frmr1 <- behav %>% filter(Genotype != "WT")  
fmr1 <- behav %>%  filter(Year == "2016")
fmr1 <- fmr1 %>%
  dplyr::select(-genoYear, -genoAPA, -genoAPATrainSessionCombo, 
                -genoAPATrainSessionComboInd, -EntrShockDiff)
names(fmr1)

fmr1$TrainSessionComboNum <- ifelse(grepl("Hab", fmr1$TrainSessionCombo), "1", 
                                     ifelse(grepl("T1", fmr1$TrainSessionCombo), "2",
                                            ifelse(grepl("T2", fmr1$TrainSessionCombo), "3",
                                                   ifelse(grepl("T3", fmr1$TrainSessionCombo), "4", 
                                                          ifelse(grepl("Retest", fmr1$TrainSessionCombo), "5",
                                                                 ifelse(grepl("T4_C1", fmr1$TrainSessionCombo), "6",
                                                                        ifelse(grepl("T5_C2", fmr1$TrainSessionCombo), "7",
                                                                               ifelse(grepl("T6_C3", fmr1$TrainSessionCombo), "8",
                                                                                      ifelse(grepl("Retention", fmr1$TrainSessionCombo), "9",NA)))))))))
fmr1$TrainSessionComboNum <- as.numeric(as.character(fmr1$TrainSessionComboNum))
fmr1 <- fmr1[c(1:18,58,19:57)]  
fmr1 <- fmr1 %>% droplevels()

fmr1$APA <- revalue(fmr1$APA, c("Conflict" = "conflict")) 
fmr1$APA <- revalue(fmr1$APA, c("Yoked" = "control")) 
fmr1$APA <- revalue(fmr1$APA, c("Same" = "consistent")) 

fmr1$APA2 <- revalue(fmr1$APA2, c("Conflict" = "conflict")) 
fmr1$APA2 <- revalue(fmr1$APA2, c("YokedConflict" = "controlconflict")) 
fmr1$APA2 <- revalue(fmr1$APA2, c("Same" = "consistent")) 
fmr1$APA2 <- revalue(fmr1$APA2, c("YokedSame" = "controlconsistent"))

write.csv(fmr1, '../data/behavior/fmr1.csv', row.names = F)

