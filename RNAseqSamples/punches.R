#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("reshape2")

setwd("~/Github/BehavEphyRNAseq/RNAseqSamples")
punches<-read.csv("punches_060915.csv", header=TRUE)
str(punches)
animals<-read.csv("animals_maddy.csv", header=TRUE)
str(animals)
CA1toDGphotos <- read.csv("CA1toDGphotos.csv", header=TRUE, sep=",")

#combine mine and maddy's notes
all <- full_join(animals, punches)
CA1toDGphotos <- full_join(animals, CA1toDGphotos)

#group and count number of animals
animals_good <- animals %>% 
  group_by(Paradigm, Behavior, E.phy) %>%
  filter(Behavior == 'Good', E.phy == 'Yes')%>%
  tally()
head(animals_good)

#group and count number of punches for "good" animals
num_punches <- all %>%
  select(Mouse, Paradigm, APA, Conflict, Behavior, E.phy, Slice, Punch) %>%
  filter(Punch %in% c("DG", "CA1", "CA2", "CA3", "CA4" )) %>%
  filter(Behavior == 'Good', E.phy == 'Yes') %>%
  group_by(Punch, Paradigm) %>%
  summarize(num_punches=length(Slice)) 
num_punches <- spread(num_punches, Punch, num_punches) 

#show punch location for "good" animals
location <- CA1toDGphotos %>%
  select(Mouse, Paradigm, Behavior, E.phy, L1, L2, L3, L4) %>%
  filter(Behavior == 'Good', E.phy == 'Yes', L1 > 1)



write.csv(animals_good, "animals_good.csv", row.names=FALSE)
write.csv(num_punches, "num_punches.csv", row.names=FALSE)
write.csv(all, "all.csv", row.names=FALSE)
write.csv(location, "location.csv", row.names=FALSE)

