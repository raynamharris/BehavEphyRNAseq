#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("plyr")
library("reshape2")

#script saved in ~/Github/BehavEphyRNAseq/RNAseqSamples ----
setwd("~/Github/BehavEphyRNAseq/data/sample_info")

#read raw data -----
punches <- read.csv("punches.csv", header=TRUE)
animals <- read.csv("animals.csv", header=TRUE)

#combine mine and maddy's notes -----
full <- join(animals, punches, by = "Mouse", type = "full", match = "all")
str(full)

## used to create the RNAseq printed program ----
forRNAseq <- full %>%
  filter(Behavior %in% c("Good","Cage", "None"), best.location >= 1, 
         Punch %in% c("CA1","CA3", "DG")) %>%
  select(Tube, RNAseqID, Group, Mouse, Behavior, E.phy, Punch, 
         Slice, best.location, jobnumber, Random) %>% 
  arrange(Random) 
# write.csv(forRNAseq, "forRNAseq.csv", row.names=FALSE)

## used to find some samples that I don't need that can be used for practice ---
MaxwellRSCtest <- full %>%
  filter(Behavior %in% c("Bad","No"), E.phy %in% c("NG","No"), 
         RNAisolationdate != "?") %>%
  select(Tube, Group, Mouse, Behavior, E.phy, Punch, Slice, RNAisolationdate) %>% 
  arrange(RNAisolationdate, Mouse, Punch, Slice) 

## used to tell Maddy which animals have RNAseq data ---
RNAseqAnimals <- full %>%
  filter(jobnumber %in% c("JA16268","JA16443", "JA16444")) %>%
  filter(Mouse.short != 100) %>% filter(Mouse.short != 101) %>%
  distinct(Mouse, Date, Conflict, APA, Behavior, E.phy)
str(RNAseqAnimals)
# write.csv(RNAseqAnimals, "RNAseqAnimals.csv", row.names=FALSE)

JA16444samples <- full %>%
  filter(jobnumber == "JA16444") %>%
  distinct(RNAseqID, Tube, Mouse, Genotype, Conflict, APA, Group, Behavior, E.phy, Punch, Slice, Date, jobnumber)
str(JA16444samples)
tail(JA16444samples)
# write.csv(JA16444samples, "JA16444samples.csv", row.names=FALSE)

