#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("reshape2")

setwd("C:/Users/RMH/Dropbox/BehavEphyRNAseq/RNAseqSamples")
punches<-read.csv("punches_060915.csv", header=TRUE)
summary(punches)

#mice: make new datset with just "Mouse" "Date" "RNAisolationdate" "Slice.collector", "location" , order by date
#CA1toDG: make a wide dataset that counts number of punches, delete irrlevant columns, order by date
punches %>%
  select(Mouse, Slice, Punch, notes.mouse, photos, Date, RNAisolationdate, Slice.collector) %>%
  group_by(Mouse, Punch, notes.mouse, photos, Date, Slice.collector) %>%
  summarize(num.slice=length(Slice)) %>% 
  as.data.frame()->mice
CA1toDG<-spread(mice, Punch, num.slice)
CA1toDG[is.na(CA1toDG)]<-0
CA1toDG$"CA1?" <- NULL
CA1toDG$"other" <- NULL
CA1toDG$"Pr" <- NULL
CA1toDG <- CA1toDG[order(as.Date(CA1toDG$Date, format="%m/%d/%Y")),]
mice <- mice[order(as.Date(mice$Date, format="%m/%d/%Y")),]

#write to file
write.csv(mice, file="mice.csv")
write.csv(CA1toDG, file="CA1toDG.csv")


## Some summary figures to look at all samples
counts <- table(punches$Slice)
barplot(counts, main="Total Slices", 
        xlab="Number of Slices per mouse")
counts <- table(punches$Punch)
barplot(counts, main="Total Slices", 
        xlab="Number of Regions per mouse")
counts <- table(punches$Mouse)
barplot(counts, main="Total Slices", 
        xlab="Mouse")
counts <- table(punches$Slice, punches$Punch)
barplot(counts, main="Regions per Slice",
        xlab="Region", 
        legend = rownames(counts))
