## Data Wrangling for Student Projects

## All punches and tissues that Rayna has collected during the  2015 and 2016 NS&B course
## are described in the same csv file in BehavPhysRNAseq repo
## I'll use this script to create a file for the students and save it in the qPCR-mouse repo
library(dplyr)
library(plyr)
library(tidyr)

## 1. Read and clean data
## 2. Wrangle data for students
## 3. select samples for exercises
## 4. select huntington's mice
## 5. select fmr1 mice

## 1. read data
setwd("~/Github/BehavEphyRNAseq/data/sample_info")
punches <- read.csv("punches.csv", header=TRUE, sep="," )
animals <- read.csv("animals.csv", header=TRUE, sep="," )
str(punches)
str(animals)
str(punches$Mouse)
str(animals$Mouse)

## clean data 
cleanpunches <- punches
cleanpunches$Mouse <- as.factor(cleanpunches$Mouse)
cleanpunches$Date <- as.Date(cleanpunches$Date, "%m/%d/%y")
cleanpunches$year <- as.factor(cleanpunches$year)
cleanpunches$year <- as.factor(cleanpunches$year)
cleanpunches$RNAisolationdate <- as.Date(cleanpunches$RNAisolationdate, "%m/%d/%y")
cleanpunches$Punchsize <- as.factor(cleanpunches$Punchsize)
cleanpunches$Tube <- as.character(cleanpunches$Tube)
cleanpunches$ready.time <-  strftime(cleanpunches$ready.time,"%H:%M")
cleanpunches$year <- as.character(cleanpunches$year)
str(cleanpunches)

## join animals and punches using the mouse id
cleanpunches <- cleanpunches %>% 
  left_join(animals, by="Mouse")
  
## 3. Wrangle for students
## inludes Maddies Test samples
NSBpunches <- cleanpunches %>%
  filter(Purpose == "students", Punch != "CA1?", Mouse != "???", Mouse != "15-310", Mouse != "15-331") %>%
  select(Tube, Mouse, Slice, Punch, Date, 
         ready.time, start.time, Slice.collector,
         punch.location, slice.location, storagebox, notes.x,
         Genotype, Conflict, APA, Group, collector, Purpose) %>%
  arrange(Mouse)
str(NSBpunches)
View(NSBpunches)

## Rayna Tests
RaynaTest <- NSBpunches %>%
  filter(Punch != "slice", collector == "MK", Punch %in% c("18", "20", "19")) %>%
  select(Tube, Mouse, Genotype, Group, Punch, Slice, Tube, storagebox) %>%
  arrange(Punch, Mouse)
View(RaynaTest)
write.csv(RaynaTest, "RaynaTest.csv")

## 4. Exercise Animals
PracticeMice <- NSBpunches %>%
  filter(Punch != "slice", collector == "MK") %>%
  select(Mouse, Genotype, Group, Punch, Tube) %>%
  arrange(Mouse)
View(PracticeMice)



## 5. just the huntington comparisons
HttMice <- NSBpunches %>%
  filter(Date < "2016-07-28", Punch != "slice", collector == "Cycle3") %>%
  distinct(Mouse, Genotype, Group) %>%
  arrange(Mouse)

HttGroups <- NSBpunches %>%
  filter(Date < "2016-07-28", Punch != "slice", collector == "Cycle3") %>%
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

