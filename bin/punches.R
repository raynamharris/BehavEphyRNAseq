#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("reshape2")

#script saved in ~/Github/BehavEphyRNAseq/RNAseqSamples
setwd("~/Github/BehavEphyRNAseq/data/sample_info")

#read raw data
punches <- read.csv("punches_060915.csv", header=TRUE)
animals <- read.csv("animals_maddy.csv", header=TRUE)

#combine mine and maddy's notes
full <- join(animals, punches, by = "Mouse", type = "full", match = "all")

#select the best samples
best <- full %>%
  filter(Behavior %in% c("Good","Cage"), best.location > 1, Punch %in% c("CA1","CA3", "DG")) %>%
  select(Paradigm, Mouse, Behavior, E.phy, Punch, photos, Tube, Slice, location.1,best.location, 
         notes,	jobnumber,	RNAseqID) %>% 
  arrange(Mouse, Slice, Punch) 

nophoto_good_yes <- full %>%
  filter(Behavior %in% c("Good","Cage"), photos != "slice photos", 
         Punch %in% c("CA1","CA3", "DG"), E.phy == "Yes", jobnumber != "JA16268") %>%
  select(Paradigm, Mouse, Behavior, E.phy, Punch, photos, Tube, Slice,	
         notes, jobnumber,	RNAseqID) %>% 
  arrange(Mouse, Slice, Punch) 

forRNAseq <- full %>%
  filter(Behavior %in% c("Good","Cage", "None"), best.location >= 1, Punch %in% c("CA1","CA3", "DG")) %>%
  select(Tube, RNAseqID, Paradigm, Mouse, Behavior, E.phy, Punch, Slice, best.location, 
         jobnumber, Random) %>% 
  arrange(Random) 
write.csv(forRNAseq, "forRNAseq.csv", row.names=FALSE)







#combine photos with maddies
CA1toDGphotos <- read.csv("CA1toDGphotos.csv", header=TRUE, sep=",")

#group and count number of animals
animals_good <- animals %>% 
  group_by(Paradigm, Behavior, E.phy) %>%
  filter(Behavior == 'Good', E.phy == 'Yes')%>%
  tally()
head(animals_good)
#write.csv(animals_good, "animals_good.csv", row.names=FALSE)

#show punch location for "good" animals
location <- CA1toDGphotos %>%
  select(Mouse, Paradigm, Behavior, E.phy, L1, L2, L3, L4) %>%
  filter(Behavior == 'Good', E.phy == 'Yes', L1 > 1)
#write.csv(location, "location.csv", row.names=FALSE)

