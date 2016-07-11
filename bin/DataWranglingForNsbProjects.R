## Data Wrangling for Student Projects

## All punches and tissues that Rayna has collected during the  2015 and 2016 NS&B course
## are described in the same csv file in BehavPhysRNAseq repo
## I'll use this script to create a file for the students and save it in the qPCR-mouse repo

library(dplyr)

## 1. Read data
## 2. Wrangle/clean for students


## 1. read data
setwd("~/Github/BehavEphyRNAseq/data/sample_info")
punches <- read.csv("punches.csv", header=TRUE, sep="," )
animals <- read.csv("animals.csv", header=TRUE, sep="," )
str(punches)
str(animals)
str(punches$Mouse)
str(animals$Mouse)

## 2. clean data & join
cleanpunches <- punches
cleanpunches$Date <- as.Date(cleanpunches$Date, "%m/%d/%y")
cleanpunches$year <- as.factor(cleanpunches$year)
cleanpunches$year <- as.factor(cleanpunches$year)
cleanpunches$RNAisolationdate <- as.Date(cleanpunches$RNAisolationdate, "%m/%d/%y")
cleanpunches$Punchsize <- as.factor(cleanpunches$Punchsize)
cleanpunches$Tube <- as.character(cleanpunches$Tube)
cleanpunches$ready.time <-  strftime(cleanpunches$ready.time,"%H:%M")
cleanpunches$year <-  as.character(cleanpunches$year)
str(cleanpunches)

cleanpunches <- cleanpunches %>% 
  left_join(animals, by="Mouse")
  
## 3. Wrangle for students

NSBpunches <- cleanpunches %>%
  filter(Purpose == "students", Punch != "CA1?") %>%
  select(Tube, Mouse, Slice, Punch, photos, Date, 
         ready.time, start.time, Slice.collector,
         punch.location, slice.location, storagebox, notes.x,
         Strain, Conflict, APA, Paradigm, behavior, Purpose)
str(NSBpunches)



NSBpunchesMouse <- NSBpunches %>%
  filter(Date < "2016-07-28", Punch != "slice") %>%
  distinct(Mouse, Date, Paradigm, behavior, Purpose) 
NSBpunchesSlice <- NSBpunches %>%
  filter(Punch == "slice") %>%
  distinct(Mouse, Punch, Date, Paradigm, behavior, storagebox) 
NSBpunchesCA1 <- NSBpunches %>%
  filter(Punch == "CA1") %>%
  distinct(Mouse, Punch, Date, Paradigm, behavior, storagebox) 
NSBpunchesCA2 <- NSBpunches %>%
  filter(Punch == "CA2") %>%
  distinct(Mouse, Punch, Date, Paradigm, behavior, storagebox) 
NSBpunchesCA3 <- NSBpunches %>%
  filter(Punch == "CA3") %>%
  distinct(Mouse, Punch, Date, Paradigm, behavior, storagebox) 
NSBpunchesCA4 <- NSBpunches %>%
  filter(Punch == "CA4") %>%
  distinct(Mouse, Punch, Date, Paradigm, behavior, storagebox) 
NSBpunchesDG <- NSBpunches %>%
  filter(Punch == "DG") %>%
  distinct(Mouse, Punch, Date, Paradigm, behavior, storagebox) 

