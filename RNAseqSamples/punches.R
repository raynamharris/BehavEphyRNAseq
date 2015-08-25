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

#make new datset with just "Mouse" "Date" "Slice.collector", "location"  
