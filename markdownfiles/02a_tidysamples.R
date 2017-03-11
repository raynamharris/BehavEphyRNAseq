#install.packages("tidyr", dependencies=TRUE)
library("tidyr")
library("dplyr")
library("plyr")
library("reshape2")

setwd("~/Github/BehavEphyRNAseq/markdownfiles")

#read raw data -----
punches <- read.csv("../data/rnaseq/punches.csv", header=TRUE)
animals <- read.csv("../data/rnaseq/animals.csv", header=TRUE)

# combine mine and maddy's notes -----
full <- join(animals, punches, by = "Mouse", type = "full", match = "all")
full <- full %>% # filter bad samples
  filter(Mouse != "", Mouse != "???") %>% droplevels()
str(full)


## used to create the RNAseq printed program ----
forRNAseq <- full %>%
  filter(Behavior %in% c("Good","Cage", "None"), best.location >= 1, 
         Punch %in% c("CA1","CA3", "DG")) %>%
  select(Tube, RNAseqID, Group, Mouse, Behavior, E.phy, Punch, 
         Slice, best.location, jobnumber, Random) %>% 
  arrange(Random) 
# write.csv(forRNAseq, "../data/rnaseq/JA16444/forRNAseq.csv", row.names=FALSE)

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
# write.csv(RNAseqAnimals, "../data/rnaseq/JA16444/RNAseqAnimals.csv", row.names=FALSE)

JA16444samples <- full %>%
  filter(jobnumber == "JA16444") %>%
  distinct(RNAseqID, Tube, Mouse, Genotype, Conflict, APA, Group, Behavior, E.phy, Punch, Slice, Date, jobnumber)
str(JA16444samples)
tail(JA16444samples)
# write.csv(JA16444samples, "../data/rnaseq/JA16444/JA16444samples.csv", row.names=FALSE)

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
#write.csv(DissociationTest, "../data/rnaseq/JA16444/sampleinfo.csv", row.names=FALSE)

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
#write.csv(DissociationTest2, "../data/rnaseq/JA16444/sampleinfo2.csv", row.names=FALSE)

## for WT with mutliple jobs
WT2015samples <- full %>%
  filter(jobnumber %in% c("JA16268", "JA16444")) %>%
  distinct(RNAseqID, Tube, Mouse, Genotype, Conflict, APA, Group, Behavior, E.phy, Punch, Slice, Date, jobnumber) %>% droplevels()
str(WT2015samples)
tail(WT2015samples)
#write.csv(WT2015samples, "../data/rnaseq/JA16444/WT2015samples.csv", row.names=FALSE)



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

#write.csv(homecage, "../data/rnaseq/JA16444/homecage.csv", row.names=FALSE)


########


############# Summer 2016 Samples for RNAseq
## first, make a spreadsheet for recording photo analysis 
#summer2016photos <- full %>%
#  filter(Year == "Summer2016", Purpose == "Collaboratorium") %>%
#  distinct(Mouse, Slice) %>% droplevels()
# This output was saved as a file and used for doing the photo analysis

## then, use the allen brain atlast to annotate the photos
## then read in photo results 
summer2016photos <- read.csv("../data/rnaseq/JA17009/summer2016photos.csv", header=T, stringsAsFactors = F)

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
#write.csv(summer2016forRNAseq, "../data/rnaseq/JA17009/summer2016forRNAseq.csv", row.names=FALSE)


### Adding the job number to the JA17009 samples
# make new vector with all the JA17009 samples
JA17009 <- lapply(summer2016forRNAseq$Tube, as.character) 
# use if else to add the job number for the above sample to the full data frame, but first make it a character
full$jobnumber <- as.character(full$jobnumber)
full$jobnumber <- ifelse(full$Tube %in% JA17009, "JA17009", full$jobnumber)
full$jobnumber <- as.factor(full$jobnumber)

full$RNAseqID <- as.character(full$RNAseqID)
full$Mouse <- as.character(full$Mouse)
full$RNAseqID <- ifelse(full$Tube %in% JA17009, full$Mouse, full$RNAseqID) 
full$RNAseqID <- as.factor(full$RNAseqID)
full$Mouse <- as.factor(full$Mouse)
str(full)


###################
## creating a data frame with just the samples that have been sequeneded that I want to analysize
# subset the data
tidysamples <- full %>%
  filter(grepl("JA16268|JA16444|JA17009", jobnumber))  %>% droplevels()
tidysamples <- tidysamples[c(1:3,5:14,16:22,29:30,34,4)] 
tidysamples <- tidysamples[c(22,1,23,3,21,13,8,4:7,9:12,14:20,24)] 

# make better column descriptions for groups
tidysamples$Group <- as.character(tidysamples$Group)
tidysamples$notes <- as.character(tidysamples$notes)
tidysamples$Group <- as.character(tidysamples$Group)
tidysamples$method <- ifelse(grepl("maddy punch", tidysamples$notes), "control", 
                    ifelse(grepl("maddy FACS", tidysamples$notes), "dissociated", "homogenized"))

tidysamples$Group <- ifelse(grepl("maddy punch", tidysamples$notes), "homecage", 
                            ifelse(grepl("maddy FACS", tidysamples$notes), "homecage", 
                                   ifelse(grepl("Cage", tidysamples$Behavior), "homecage",
                                          ifelse(grepl("Yoked", tidysamples$Group), "control", 
                                                 ifelse(grepl("NoConflict", tidysamples$Group), "consistent",
                                                        ifelse(grepl("Conflict", tidysamples$Group), "conflict",tidysamples$Group))))))


# make better column descriptions for possibly dodgy samples
tidysamples$dodgy <- ifelse(grepl("NG", tidysamples$E.phy), "ephys",
                            ifelse(grepl("CYRSTALS", tidysamples$notes), "slice",
                                   ifelse(grepl("late", tidysamples$notes), "slice",
                                          ifelse(grepl("350", tidysamples$SliceSize), "slice", "allgood"))))

# make better column for time of day collected
str(tidysamples$start.time)
tidysamples$start.time <- as.character(tidysamples$start.time)

tidysamples$daytime <- ifelse(grepl("10:|11:", tidysamples$start.time), "beforenoon",
                            ifelse(grepl("4:30", tidysamples$start.time), "earlyAM",
                                   ifelse(grepl("12:|14:|15:", tidysamples$start.time), "afternoon",
                                          ifelse(grepl("17:|18:|19:", tidysamples$start.time), "evening",
                                                 ifelse(grepl("20:|23:0", tidysamples$start.time), "nighttime", "norecord")))))




## drop some JA16268 samples and rename
tidysamples$RNAseqID <- as.character(tidysamples$RNAseqID)
tidysamples <- tidysamples %>%
  filter(!grepl("DG_D|CA1_D", RNAseqID))  %>% droplevels()
tidysamples$RNAseqID <- ifelse(grepl("142C_DG_S", tidysamples$RNAseqID), "142C_DG",
                              ifelse(grepl("142C_CA1_S", tidysamples$RNAseqID), "142C_CA1",
                                     ifelse(grepl("143C_CA1_S", tidysamples$RNAseqID), "143C_CA1",
                                            ifelse(grepl("143C_DG_S", tidysamples$RNAseqID), "143C_DG",tidysamples$RNAseqID))))
                                     
## correct mouse name for 15-100 and 15-101 samples
tidysamples$Mouse <- as.character(tidysamples$Mouse)
tidysamples$Mouse <- ifelse(grepl("15-101", tidysamples$Mouse), "15-100", tidysamples$Mouse)
             
## final drop columns
tidysamples <- tidysamples[c(1:7,10:11,24:26,15:16)] 


## rename "Punch" to "Region"
tidysamples <- rename(tidysamples, c("Punch"="Region"))
tidysamples <- rename(tidysamples, c("method"="Method"))


## write out clearn data file will all samples
#write.csv(tidysamples, "../data/rnaseq/tidysamples.csv", row.names=FALSE)


