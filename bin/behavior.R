# Part 1: Reading and analyzing beahvior and physiology data

## load libraries -----
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!
library(reshape2) #@ for melting dataframe
library(ggdendro) ## for dendrograms!!
library(magrittr) ## to use the weird pipe
library(gplots)

## read and wrangle the data ----

## read the data 
setwd("~/Github/BehavEphyRNAseq/data/sample_info/")
behav <- read.csv("APA_2013-2016.csv", header=TRUE, stringsAsFactors = FALSE, na.strings = c("", "ND", "N/A"))
str(behav)

## rename columns 
names(behav)[18] <- "NumEntrances"
names(behav)[23] <- "NumShock"
names(behav)[26] <- "Speed_cm_s"
names(behav) # check all good
head(behav)

## remove bad/extra rows
behav <- behav %>% dplyr::filter(Year %in% c("2012", "2013", "2014", "2015", "2016", "2017")) ## remove extra headers
behav <- behav %>% dplyr::filter(TrainSession != "C4?") ## remove ?C in train train animal
behav <- behav %>% dplyr::filter(!is.na(TotalTime.s.)) ## remove file with no data


## make characters either factors or numbers, as appropriate
cols = c(1:4,6:14)
cols2 = c(15:58)
behav[,cols] %<>% lapply(function(x) as.factor(as.character(x)))
behav[,cols2] %<>% lapply(function(x) as.numeric(as.character(x)))
str(behav)
names(behav)

## rename some values -----
behav$Genotype <- revalue(behav$Genotype, c("FMR1" = "FMR1-KO")) 
behav$Genotype <- revalue(behav$Genotype, c("FMR1?WT?" = "FMR1-KO")) 
behav$TrainGroup <- revalue(behav$TrainGroup, c("control" = "untrained")) 
behav$TrainSequence[is.na(behav$TrainSequence)] <- "untrained" ## make NA more meaningful
behav$TrainSequence <- as.factor(behav$TrainSequence)

## create combinatorial factor columns 
behav$genoYear <- as.factor(paste(behav$Genotype,behav$Year,sep="_"))
behav$genoYear <- factor(behav$genoYear, 
                    levels = c("WT_2013", "WT_2014", "FMR1-KO_2014", 
                               "WT_2015", "WT_2016", "FMR1-KO_2016"))

behav$APA <- as.factor(paste(behav$TrainSequence,behav$TrainGroup,sep="_"))
behav$APA <- revalue(behav$APA, c("train_trained" = "trained")) 
behav$APA <- revalue(behav$APA, c("untrained_untrained" = "untrained")) 
behav$APA <- revalue(behav$APA, c("train-conflict_trained" = "trained_conflict")) 
behav$APA <- revalue(behav$APA, c("train-conflict_yoked" = "yoked_conflict")) 
behav$APA <- revalue(behav$APA, c("train-train_trained" = "trained_trained")) 
behav$APA <- revalue(behav$APA, c("train-train_yoked" = "yoked_trained")) 
behav$APA <- factor(behav$APA, 
                        levels = c("untrained", "yoked_trained", 
                                   "yoked_conflict", "trained", 
                                   "trained_trained", "trained_conflict"))

behav$genoAPA <- as.factor(paste(behav$Genotype,behav$APA, sep="_"))
behav$genoAPA <- factor(behav$genoAPA, 
                        levels = c("WT_untrained", "WT_trained", 
                                   "WT_yoked_trained", "WT_trained_trained", 
                                   "WT_yoked_conflict", "WT_trained_conflict", 
                                   "FMR1-KO_untrained", "FMR1-KO_trained",
                                   "FMR1-KO_yoked_trained", "FMR1-KO_trained_trained",
                                   "FMR1-KO_yoked_conflict", "FMR1-KO_trained_conflict"))

behav$genoAPAsession  <- as.factor(paste(behav$genoAPA,behav$TrainSession, sep="_"))
behav$genoAPAsessionDay  <- as.factor(paste(behav$genoAPAsession,behav$Day, sep="_"))
behav$genoAPAsessionDayInd <- as.factor(paste(behav$genoAPAsessionDay,behav$ID, sep="_"))

behav$TrainSessionCombo <- behav$TrainSession
levels(behav$TrainSessionCombo)
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C1" = "T4_C1")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T4" = "T4_C1")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C2" = "T5_C2")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T5" = "T5_C2")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C3" = "T6_C3")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T6" = "T6_C3")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T7" = "T+_C+")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("T8" = "T+_C+"))
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C4" = "T+_C+"))
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C5" = "T+_C+")) 
behav$TrainSessionCombo <- revalue(behav$TrainSessionCombo, c("C6" = "T+_C+")) 
levels(behav$TrainSessionCombo)
behav$TrainSessionCombo <- factor(behav$TrainSessionCombo, 
                        levels = c("Hab", "T1","T2","T3","Retest",
                                   "T4_C1","T5_C2","T6_C3","T+_C+","Retention"))

behav$TrainSessionComboDay <- as.factor(paste(behav$TrainSessionCombo, behav$Day, sep="_"))

behav$genoAPAsessionCombo <- as.factor(paste(behav$genoAPA, behav$TrainSessionCombo, sep="_"))

behav$pair1 <- as.factor(paste(behav$ID,behav$TrainSessionComboDay, sep="_"))
behav$pair2 <- as.factor(paste(behav$PairedPartner,behav$TrainSessionComboDay, sep="_"))

## reorders dataframe 
names(behav)
behav <- behav[c(2,59:69,3:9,1,10:58)]  
names(behav)

## subset the data -----
maddy <- behav %>% filter(Experimenter == "Maddy") 
jma <- behav %>% filter(Experimenter != "Maddy") 
wt <- behav %>% filter(Genotype == "WT")
frmr1 <- behav %>% filter(Genotype != "WT")  

## make a df to look at number of shock actually received by the yoked animals 
names(behav)
yoked <- behav %>% 
  select(pair1, pair2, TrainGroup, Experimenter, NumShock, NumEntrances) %>% 
  filter(Experimenter == "Maddy") %>% droplevels()
yokedpair <- behav %>% 
  select(pair1, pair2, TrainGroup, Experimenter, NumShock, NumEntrances) %>% 
  filter(Experimenter == "Maddy") %>% droplevels()
## rename columns 
names(yoked)[1] <- "train"
names(yoked)[2] <- "yoke"
names(yokedpair)[1] <- "yoke"
names(yokedpair)[2] <- "train"
## join an caluclate some values
yokedyokedpair <- full_join(yoked, yokedpair, by = "yoke")
yokedyokedpair <- yokedyokedpair %>%  
  mutate(shockdiff = NumShock.x - NumShock.y) %>%  
  mutate(entrdiff = NumEntrances.x - NumEntrances.y)

## create the ggplot color palettes ----

## for APA palette
#40004b purple train conflict
#762a83 purple train train
#af8dc3 purpletrain 
#7fbf7b green untrained
#1b7837 greenuntrained trained
#00441b green untrained conflict
APApalette <- (values=c("#7fbf7b","#1b7837","#00441b","#af8dc3", "#762a83","#40004b"))
APApaletteSlim <- (values=c("#1b7837","#00441b","#762a83","#40004b"))

## using the same red black orange grey scale JMA used
JMPalette <- c('black','grey50','red','darkorange')
WTPalette <- c('black','red')
FMR1Palette <- c('grey50','darkorange')

## ggplots across sessions!! -----

## ggpots of TimeTarget 
behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=TimeTarget, color=APA)) + 
  stat_smooth() + theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApalette) + 
  scale_y_continuous(name="Time spent in the shock zone") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "T4/C1",
                              "6" = "T5/C2", "7" = "T6/C3")) +
  facet_wrap(~genoYear) 

behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=TimeTarget, color=APA)) + 
  stat_smooth() + theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApaletteSlim) + 
  scale_y_continuous(name="Time spent in the shock zone") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear) 

## ggpots of time to 2nd entrance 
behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=Time2ndEntr, color=APA)) + 
  stat_smooth() + theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApalette) + 
  scale_y_continuous(name="Time to 2nd Entrance") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "T4/C1",
                              "6" = "T5/C2", "7" = "T6/C3")) +
  facet_wrap(~genoYear) 

behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=Time2ndEntr, color=APA)) + 
  stat_smooth() + theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApaletteSlim) + 
  scale_y_continuous(name="Time to 2nd Entrance") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear)

## num entrances 
behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=NumEntrances, color=APA)) + 
  stat_smooth() + theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApalette) + 
  scale_y_continuous(name="Number of Entrances") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "T4/C1",
                              "6" = "T5/C2", "7" = "T6/C3")) +
  facet_wrap(~ genoYear) 

behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=NumEntrances, color=APA)) + 
  stat_smooth() + theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApaletteSlim) + 
  scale_y_continuous(name="Number of Entrances") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear)

## pTimeOpp 
behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=pTimeOPP, color=APA)) + 
  stat_smooth() + theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApalette) + 
  scale_y_continuous(name="pTimeOPP") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "T4/C1",
                              "6" = "T5/C2", "7" = "T6/C3")) +
  facet_wrap(~ genoYear)

behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=TotalPath.m., color=APA)) + 
  stat_smooth() + theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApaletteSlim) + 
  scale_y_continuous(name="Total Path") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear)

## pTimeOpp 
behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=TotalPath.m., color=APA)) + 
  stat_smooth() + theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApalette) + 
  scale_y_continuous(name="Total Path") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "T4/C1",
                              "6" = "T5/C2", "7" = "T6/C3")) +
  facet_wrap(~ genoYear)

behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=Speed1stEntr.cm.s., color=APA)) + 
  stat_smooth() + theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=APApaletteSlim) + 
  scale_y_continuous(name="Speed to 1st Entrance") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear)

## heatmap of data (not correlations, but raw/scaled data) ----
### melt to make long 
behav_long <- melt(behav, id=c("ID","APA","genoAPA","genoAPAsession","genoAPAsessionDay", "genoYear", "genoAPAsessionCombo",
                               "genoAPAsessionDayInd","TrainSessionCombo" ,"Genotype", "TrainSessionComboDay",
                               "TrainProtocol","TrainSequence","TrainGroup","Day","TrainSession",
                               "ShockOnOff","Year","PairedPartner","Experimenter",
                               "Housing","TestLocation","filename", "pair1", "pair2"))
behav_long <- filter(behav_long, !grepl("TotalTime.s|p.miss", variable )) %>% droplevels()
#behav_long <- filter(behav_long, !grepl("Retention|Retest", genoAPAsession )) %>% droplevels()
str(behav_long)

## now widen then lengthen to get group averages
behav_long_genoAPA <- dcast(behav_long, genoAPA ~ variable, value.var= "value", fun.aggregate=mean)
head(behav_long_genoAPA)

## scale columns
rownames(behav_long_genoAPA) <- behav_long_genoAPA$genoAPA    # set $genoAPAsession as rownames
behav_long_genoAPA[1] <- NULL
behav_long_genoAPA <- scale(behav_long_genoAPA)
head(behav_long_genoAPA)

## heatmap clusterd
heatpalette <- colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 100)
heatmap.2(behav_long_genoAPA,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,13),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          #scale = c("column"),
          RowSideColors = c("#7fbf7b", "#af8dc3", "#1b7837", "#762a83", "#00441b", "#40004b",
                                   "#7fbf7b", "#af8dc3", "#1b7837", "#762a83", "#00441b", "#40004b"))

## correlation matrix and plots ----
ggplot(behav_long_genoAPA, aes(x = genoAPA, y = variable, fill = value)) + 
  geom_tile() + coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## first, make the data a matrix with genoAPAsessionInd as the row names
behav_matrix <- behav   #create new dataframe
rownames(behav_matrix) <- behav_matrix$genoAPAsessionDayInd     # set $genoAPAsessionInd as rownames
names(behav_matrix)
behav_matrix <- behav_matrix[-c(1:27)] #delete all non-numeric columns and p.miss and TotalTime
head(behav_matrix)
str(behav_matrix)

behav_matrix_cormat <- round(cor(behav_matrix),2) # compute correlations

##  some necessary functions
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Get lower and lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

## calculate the corerelation matrix
cormat <- reorder_cormat(behav_matrix_cormat)
lower_tri <- get_lower_tri(cormat)
melted_cormat <- melt(lower_tri, na.rm = TRUE)
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() 


### unused plots ----
behav %>% 
  filter(TrainSessionCombo %in% c("T3"))  %>%  droplevels() %>%
  ggplot(aes((x=TrainSessionCombo), y=pTimeTarget, fill=APA)) + 
  geom_boxplot() + theme_bw() + scale_fill_manual(values=APApalette) +
  facet_wrap(~Genotype+Year) +
  scale_y_continuous(name="Probability of being in the shock zone") 

behav %>% 
  filter(TrainSessionCombo %in% c("T4_C1"))  %>%  droplevels() %>%
  ggplot(aes((x=TrainSessionCombo), y=pTimeTarget, fill=APA)) + 
  geom_boxplot() + theme_bw() + scale_fill_manual(values=APApaletteSlim) +
  facet_wrap(~Genotype+Year) +
  scale_y_continuous(name="Probability of being in the shock zone") 

