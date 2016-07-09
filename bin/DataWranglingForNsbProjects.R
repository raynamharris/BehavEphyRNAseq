## Data Wrangling for Student Projects

## All punches and tissues that Rayna has collected during the  2015 and 2016 NS&B course
## are described in the same csv file in BehavPhysRNAseq repo
## I'll use this script to create a file for the students and save it in the qPCR-mouse repo

# 1. Read data
# 2. Wrangle/clean for students

# 1. read data
setwd("~/Github/BehavEphyRNAseq/data/sample_info")
punches <- read.csv("punches.csv", header=TRUE, sep="," )
str(punches)

## 2. Wrangle/clean for students
NSBpunches <- punches
NSBpunches$Date <- as.Date(NSBpunches$Date, "%m/%d/%y")
str(punches)
NSBpunches$Date

NSBpunches$ready.time <-  strftime(NSBpunches$ready.time,"%H:%M")
