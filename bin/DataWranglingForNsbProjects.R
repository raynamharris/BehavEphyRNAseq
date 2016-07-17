## Data Wrangling for Student Projects

## All punches and tissues that Rayna has collected during the  2015 and 2016 NS&B course
## are described in the same csv file in BehavPhysRNAseq repo
## I'll use this script to create a file for the students and save it in the qPCR-mouse repo

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

## 1. Read and clean data
## 2. Wrangle data for students
## 3. select samples for exercises
## 4. select huntington's mice
## 5. select fmr1 mice

## 1. read data
setwd("~/Github/BehavEphyRNAseq/data/sample_info")
punches <- read.csv("punches.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
animals <- read.csv("animals.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
punches$mouse <- as.factor(punches$Mouse)
animals$mouse <- as.factor(animals$Mouse)
str(punches)
str(animals)

## join animals and punches using the mouse id
alldata <- punches %>% 
  left_join(animals, by="Mouse")
str(alldata)

## clean data 
alldata$Mouse <- as.factor(alldata$Mouse)
alldata$Date <- as.Date(alldata$Date, "%m/%d/%y")
alldata$Slice <- as.factor(alldata$Slice)
alldata$Punch <- as.factor(alldata$Punch)
alldata$year <- as.factor(alldata$year)
alldata$RNAisolationdate <- as.Date(alldata$RNAisolationdate, "%m/%d/%y")
alldata$Punchsize <- as.factor(alldata$Punchsize)
alldata$Tube <- as.factor(alldata$Tube)
alldata$collector <- as.factor(alldata$Tube)
alldata$Group <- as.factor(alldata$Group)
alldata$APA <- as.factor(alldata$APA)
str(alldata)


## 3. Wrangle for students
## inludes Maddies Test samples
NSBpunches <- alldata %>%
  filter(Purpose == "students", Punch != "CA1?", 
         Mouse != "???", Mouse != "15-310", Mouse != "15-331", 
         Slice != "NA") %>%
  select(Tube, Mouse, Slice, Punch, Date, 
         ready.time, start.time, Slice.collector,
         punch.location, slice.location, storagebox, notes.x,
         Genotype, Conflict, APA, Group, collector, Purpose) %>%
  arrange(Mouse)
str(NSBpunches)


## 4. Group Project: FMR1 - conflict - slices
FMR1_conflict_slices <- alldata %>%
  filter(Punch == "slice", Genotype %in% c("WT", "FMR1"), 
         Conflict == "Conflict", Slice != "5", 
         Mouse != "16-119C", Mouse != "16-119D", Mouse != "16-118C", Mouse != "16-118D") %>%
  select(Mouse, Group, Slice,Punch, Tube, storagebox) %>%
  arrange(Group, Mouse) 
View(FMR1_conflict_slices) 

FMR1_conflict_animals <- dcast(FMR1_conflict_slices, Mouse + Group ~ Slice, value.var="Slice")
View(FMR1_conflict_animals)


## 5. just the huntington comparisons
HttMice <- cleanpunches %>%
  select(Date, Mouse, Genotype, Conflict, APA, Group, Slice, Punch, collector, Tube) %>%
  filter(Punch != "slice", collector == "Cycle3") %>%
  distinct(Mouse, Genotype, Group) %>%
  arrange(Mouse)

HttGroups <- cleanpunches %>%
  filter(Punch != "slice", collector == "Cycle3") %>%
  distinct(Mouse, Genotype, Group) %>%
  arrange(Group, Genotype) %>%
  summarise(count(Group))
str(HttGroups)

HttTissues <- NSBpunches %>%
  filter(Date < "2016-07-28", Punch != "slice", collector == "Cycle3") %>%
  select(Mouse, Genotype, Group, Punch) %>%
  arrange(Group, Genotype, Punch) 
HttTissues <- (count(HttTissues, vars=c("Group","Punch")))
HttTissues <- spread(HttTissues, Punch, freq)


NSBpunchesSlice <- NSBpunches %>%
  filter(Punch == "slice") %>%
  distinct(Mouse, Punch, Date, Group, collector, Slice, Tube, storagebox) %>%
  arrange(Mouse)
View(NSBpunchesSlice)

NSBpunchesCA1 <- NSBpunches %>%
  filter(Punch == "CA1") %>%
  distinct(Mouse, Punch, Date, Group, collector) 

NSBpunchesCA2 <- NSBpunches %>%
  filter(Punch == "CA2") %>%
  distinct(Mouse, Punch, Date, Group, collector) 

NSBpunchesCA3 <- NSBpunches %>%
  filter(Punch == "CA3") %>%
  distinct(Mouse, Punch, Date, Group, collector) 

NSBpunchesCA4 <- NSBpunches %>%
  filter(Punch == "CA4") %>%
  distinct(Mouse, Punch, Date, Group, collector) 

NSBpunchesDG <- NSBpunches %>%
  filter(Punch == "DG") %>%
  distinct(Mouse, Punch, Date, Group, collector) 

