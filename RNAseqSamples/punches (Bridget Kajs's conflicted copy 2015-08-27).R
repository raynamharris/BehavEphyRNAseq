#analysis of samples collected
setwd("C:/Users/RMH/Dropbox/BehavEphyRNAseq/RNAseqSamples")
punches<-read.csv("punches_060915.csv", header=TRUE)
summary(punches)

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

#make new datset with just "Mouse" "Date" "RNAisolationdate" "Slice.collector", "location"  
#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("reshape2")

punches %>%
  select(Mouse, Slice, notes.mouse, Date, RNAisolationdate, Slice.collector) %>%
  group_by(Mouse, notes.mouse, Date, Slice.collector) %>%
  summarize(num.slice=max(Slice)) %>% 
  as.data.frame()->mice
write.csv(mice, file="mice.csv")

punches %>%
  select(Mouse, Slice, Punch, notes.mouse, Date, RNAisolationdate, Slice.collector) %>%
  group_by(Mouse, Punch, notes.mouse, Date, Slice.collector) %>%
  summarize(num.slice=length(Slice)) %>% 
  as.data.frame()->CA1toDG

CA1toDGWIDE<-spread(CA1toDG, Punch, num.slice)
CA1toDGWIDE[is.na(CA1toDGWIDE)]<-0
CA1toDGWIDE$"CA1?"<-NULL
CA1toDGWIDE$"other"<-NULL
CA1toDGWIDE$"Pr"<-NULL
write.csv(CA1toDGWIDE, file="CA1toDGWIDE.csv")
