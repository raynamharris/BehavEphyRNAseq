#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("plyr")
library("reshape2")

setwd("~/Github/BehavEphyRNAseq/bin")

#read raw data -----
punches <- read.csv("../data/rnaseq/punches.csv", header=TRUE)
animals <- read.csv("../data/rnaseq/animals.csv", header=TRUE)

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

## for DissociationTest project
DissociationTest <- full %>%
  filter(Mouse %in% c("15-100", "15-101")) %>%
  filter(jobnumber == "JA16444")  %>% droplevels()
DissociationTest <- DissociationTest[c(1:3, 10:14,16,29:32)] 
DissociationTest <- DissociationTest[c(11,4,8,7,1:3,5:6,9:10)]
names(DissociationTest)[names(DissociationTest)=="notes"] <- "Method"
DissociationTest$Method <- revalue(DissociationTest$Method, c("maddy FACS dissociate" = "dissociated")) 
DissociationTest$Method <- revalue(DissociationTest$Method, c("maddy punch" = "homogenized")) 
str(DissociationTest)
#write.csv(DissociationTest, "/Users/raynamharris/Github/DissociationTest/data/sampleinfo.csv", row.names=FALSE)

## for Bigger DissociationTest project
## includes some animals for which we have behavior data
DissociationTest2 <- full %>%
  filter(grepl("15-100|15-101|15-145|15-146|15-147|15-148", Mouse)) %>%
  filter(RNAseqID != "146A-CA3-2") %>%
  filter(jobnumber == "JA16444")  %>% droplevels()

DissociationTest2$Method <- ifelse(grepl("maddy FACS dissociate", DissociationTest2$notes), "Dissociated", 
                                     ifelse(grepl("maddy punch", DissociationTest2$notes), "Homogenized",
                                            ifelse(grepl("WT Conflict Trained", DissociationTest2$Group), "Trained",
                                                   ifelse(grepl("WT Conflict Yoked", DissociationTest2$Group), "Yoked", 
                                                          ifelse(grepl("WT NoConflict Trained", DissociationTest2$Group), "Trained",
                                                                 ifelse(grepl("WT NoConflict Yoked", DissociationTest2$Group), "Yoked",
                                                                        ifelse(grepl("WT NA NA", DissociationTest2$Group), "HomeCage",NA)))))))
DissociationTest2 <- DissociationTest2[c(1:3,12:14,16,21,29:32,38)] 
DissociationTest2 <- DissociationTest2[c(10,13,6,5,1:4,7:9)]
str(DissociationTest2)
#write.csv(DissociationTest2, "/Users/raynamharris/Github/DissociationTest/data/sampleinfo2.csv", row.names=FALSE)

## for WT with mutliple jobs
WT2015samples <- full %>%
  filter(jobnumber %in% c("JA16268", "JA16444")) %>%
  distinct(RNAseqID, Tube, Mouse, Genotype, Conflict, APA, Group, Behavior, E.phy, Punch, Slice, Date, jobnumber) %>% droplevels()
str(WT2015samples)
tail(WT2015samples)
#write.csv(WT2015samples, "WT2015samples.csv", row.names=FALSE)



## for homecage animals
homecage <- full %>%
  filter(jobnumber == "JA16444")  %>% droplevels()
homecage$Group <- as.character(homecage$Group)
homecage$notes <- as.character(homecage$notes)
homecage$Group <- ifelse(grepl("maddy punch", homecage$notes), "homogenized", 
                         ifelse(grepl("maddy FACS", homecage$notes), "dissociated", 
                                ifelse(grepl("Cage", homecage$Behavior), "homecage",
                                ifelse(grepl("Yoked", homecage$Group), "yoked", homecage$Group))))
homecage <- homecage[c(1,2,9,14,30)] 
homecage <- homecage %>%
  filter(grepl("yoked|homecage", Group)) %>% 
  filter(grepl("15-145|15-146|15-147|15-148", Mouse)) %>%   droplevels()

#write.csv(homecage, "/Users/raynamharris/Github/DissociationTest/data/homecage.csv", row.names=FALSE)


############# Summer 2016 Samples for RNAseq
## first, make a spreadsheet for recording photo analysis 
summer2016photos <- full %>%
  filter(Year == "Summer2016", Purpose == "Collaboratorium") %>%
  distinct(Mouse, Slice) %>% droplevels()
# This output was saved as a file and used for doing the photo analysis

## then, use the allen brain atlast to annotate the photos
## then read in photo results 
summer2016photos <- read.csv("summer2016photos.csv", header=T, stringsAsFactors = F)

## subset the full punch dataset and keep the relevant columns for the 2016 project
summer2016forRNAseq <- full %>%
  filter(Year == "Summer2016", Purpose == "Collaboratorium") %>%
  distinct(Tube, Mouse, Genotype, Group, Punch, Slice, Date, storagebox) %>% droplevels()

## merge the photo and the punch data
summer2016forRNAseq <- join(summer2016forRNAseq, summer2016photos, by=c("Mouse","Slice"), type = "full", match = "all")
summer2016forRNAseq$Slice <- as.integer(summer2016forRNAseq$Slice) # first make both integers
str(summer2016forRNAseq)
rm(summer2016photos) #don't need this anymore

# reset the Groups to be Yoked, Same, and Conflict
summer2016forRNAseq <- rename(summer2016forRNAseq, c("Group"="APA"))
summer2016forRNAseq$APA <- ifelse(grepl("NoConflict Trained", summer2016forRNAseq$APA), "Same", 
                                  ifelse(grepl("Yoked", summer2016forRNAseq$APA), "Yoked", "Conflict"))
summer2016forRNAseq$APA


## create a T/F column to say if the the same is from the optimal slice or not. 
summer2016forRNAseq$isbest <- ifelse(summer2016forRNAseq$Slice == summer2016forRNAseq$sliceforRNAseq, T,F) 

## filter data to only the BEST slices and CA1, CA3, and DG Samples
summer2016forRNAseq <- summer2016forRNAseq %>%
  filter(isbest == "TRUE") %>% filter(Punch %in% c("CA1")) %>% 
  filter(APA == "Yoked") %>% 
  arrange(Date) %>% droplevels()

## create an RNAseqID (can't have dashes, must be less than 10 characters)
head(summer2016forRNAseq)
summer2016forRNAseq$idealRNAseqID <- as.factor(paste(summer2016forRNAseq$Mouse,summer2016forRNAseq$Slice,sep="_"))
summer2016forRNAseq$idealRNAseqID <- gsub("-", "_", summer2016forRNAseq$idealRNAseqID, fixed = TRUE)
summer2016forRNAseq$RNAseqID <- summer2016forRNAseq$Mouse
#write.csv(summer2016forRNAseq, "~/Github/BehavEphyRNAseq/data/rnaseq/summer2016forRNAseq.csv", row.names=FALSE)


### calculate sample sizes
summer2016forRNAseqtotals <- select(summer2016forRNAseq, Genotype, APA, Punch)
summer2016forRNAseqtotals <- count(summer2016forRNAseqtotals, c('Genotype','APA', "Punch"))
summer2016forRNAseqtotals <- dcast(summer2016forRNAseqtotals, Genotype + APA ~ Punch, value.var = "freq")
head(summer2016forRNAseqtotals)
rm(summer2016forRNAseqtotals)
