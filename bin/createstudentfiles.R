## Sample Wrangling for Student Projects

## The script was written to important and organize all tissues for student qPCR projects.
## All NSBpunches and slices were collected during the 2016 NS&B course.
## The NSBanimals.csv file contains the information about each animal used.
## Students will add inforamtion about the slice location to animals.csv
## The punches.csv files contains information about all the tissues collected.

## Parts of this script
## 1. Install packages
## 2. Read and clean data
## 3. Wrangle data for students
## 4. select samples for exercises
## 5. select huntington's mice
## 6. select fmr1 mice

## 1. install packages needed
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

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
NSBanimals <- alldata %>%
  filter(Purpose == "students", Punch %in% c("CA1", "CA2", "CA3", "CA4", "DG", "slice"), Mouse != "???", Mouse != "15-310", Mouse != "15-331") %>%
  distinct(Mouse, Genotype, Conflict, APA, Group) %>%
  arrange(Mouse) 
str(NSBanimals)
View(NSBanimals)

NSBpunches <- alldata %>%
  filter(Purpose == "students", Punch %in% c("CA1", "CA2", "CA3", "CA4", "DG", "slice"), Mouse != "???", Mouse != "15-310", Mouse != "15-331") %>%
  select(Tube, Mouse, Slice, Punch, Date, 
         ready.time, start.time, Slice.collector,
         storagebox, notes.x) %>%
  arrange(Mouse)
str(NSBpunches)
View(NSBpunches)

'''
Never do again!!!
write.csv(NSBanimals, "NSBanimals.csv", row.names = FALSE)
write.csv(NSBpunches, "NSBpunches.csv", row.names = FALSE)
'''
