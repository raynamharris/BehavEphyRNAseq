#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("plyr")
library("reshape2")

#script saved in ~/Github/BehavEphyRNAseq/RNAseqSamples
setwd("~/Github/BehavEphyRNAseq/data/sample_info")

#read raw data
punches <- read.csv("punches.csv", header=TRUE)
animals <- read.csv("animals.csv", header=TRUE)

#combine mine and maddy's notes
full <- join(animals, punches, by = "Mouse", type = "full", match = "all")

#select the best samples & then the 2nd best -> the create the final output "for RNAseq
best <- full %>%
  filter(Behavior %in% c("Good","Cage"), best.location > 1, Punch %in% c("CA1","CA3", "DG")) %>%
  select(Paradigm, Mouse, Behavior, E.phy, Punch, photos, Tube, Slice, location.1,best.location, 
         notes,	jobnumber,	RNAseqID) %>% 
  arrange(Mouse, Slice, Punch) 

secondbest <- full %>%
  filter(Behavior %in% c("Good","Cage"), photos != "slice photos", 
         Punch %in% c("CA1","CA3", "DG"), E.phy == "Yes", 
         jobnumber != "JA16268") %>%
  select(Paradigm, Mouse, Behavior, E.phy, Punch, photos, Tube, Slice,	
         notes, jobnumber,	RNAseqID) %>% 
  arrange(Mouse, Slice, Punch) 

forRNAseq <- full %>%
  filter(Behavior %in% c("Good","Cage", "None"), best.location >= 1, 
         Punch %in% c("CA1","CA3", "DG")) %>%
  select(Tube, RNAseqID, Paradigm, Mouse, Behavior, E.phy, Punch, 
         Slice, best.location, jobnumber, Random) %>% 
  arrange(Random) 
write.csv(forRNAseq, "forRNAseq.csv", row.names=FALSE)

MaxwellRSCtest <- full %>%
  filter(Behavior %in% c("Bad","No"), E.phy %in% c("NG","No"), 
         RNAisolationdate != "?") %>%
  select(Tube, Paradigm, Mouse, Behavior, E.phy, Punch, Slice, RNAisolationdate) %>% 
  arrange(RNAisolationdate, Mouse, Punch, Slice) 
  