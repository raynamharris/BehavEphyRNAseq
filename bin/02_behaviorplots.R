

## load libraries -----
library(ggdendro) ## for dendrograms!!
library(magrittr) ## to use the weird pipe
library(gplots) ##for making awesome plots
library(cowplot) ## for some easy to use themes



## color palettes tips ----

## for APApalette
#40004b purple train conflict
#762a83 purple train train
#af8dc3 purpletrain 

#7fbf7b green untrained
#1b7837 green untrained trained
#00441b green untrained conflict

## for APApalette2

#40004b purple train conflict
#9970ab purple train train
#c2a5cf purpletrain 

#fdb863 orange untrained
#f1a340 orange untrained trained
#b35806 orange untrained conflict

## create the ggplot color palettes ----

APApalette <- (values=c("#7fbf7b","#1b7837","#00441b","#af8dc3", "#762a83","#40004b"))
APApaletteSlim <- (values=c("#1b7837","#00441b","#762a83","#40004b"))

APApalette2 <- (values=c("#fdb863","#f1a340","#b35806","#c2a5cf", "#9970ab","#40004b"))
APApaletteSlim2 <- (values=c("#f1a340","#b35806","#9970ab","#40004b"))

APApaletteSlim3 <- (values=c("#9970ab","#40004b"))

## using the same red black orange grey scale JMA used-----
JMPalette <- c('black','grey50','red','darkorange')
WTPalette <- c('black','red')
FMR1Palette <- c('grey50','darkorange')

## ggplots across sessions!! -----

## creating better facet label names
genoAPAnames <- c(
  "WT_2015" = "WT 2015", 
  "WT_2016" = "WT 2016",
  "FMR1-KO_2016" = "FMR1-KO 2015"
)


## ggpots of TimeTarget - saved as beahv_TimeTarget_6 or _3
ptime <- behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=pTimeTarget, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.1, size=2)  +
  theme_cowplot(font_size = 20, line_size = 1) + 
  theme(strip.background = element_blank()) +  
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  scale_colour_manual(values=APApaletteSlim2,
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                      labels=c("Yoked to Trained", "Yoked to Conflict", "Trained", "Conflict")) + 
  scale_y_continuous(name="Probability of being in the Shock Zone") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_grid(~genoYear, labeller = as_labeller(genoAPAnames))

## ggpots of time to 2nd entrance - saved as beahv_Time2ndEntrance_6 or _3
two <- behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=Time2ndEntr, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.2, size=2)  +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  theme(strip.background = element_blank()) +  
  scale_colour_manual(values=APApaletteSlim2,
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                      labels=c("Yoked to Trained", "Yoked to Conflict", "Trained", "Conflict")) + 
  scale_y_continuous(name="Time to 2nd Entrance (s)", limits = c(0, 600)) + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_grid(~genoYear, labeller = as_labeller(genoAPAnames))

one <- behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=Time1stEntr, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.2, size=2)  +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  theme(strip.background = element_blank()) +  
  scale_colour_manual(values=APApaletteSlim2,
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                      labels=c("Yoked to Trained", "Yoked to Conflict", "Trained", "Conflict")) + 
  scale_y_continuous(name="Time to 1st Entrance (s)", limits = c(0, 600)) + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_grid(~genoYear, labeller = as_labeller(genoAPAnames))

## num entrances - saved as beahv_NumEntrance_6 or _3
entrances <- behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=NumEntrances, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.2, size=2)  +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  theme(strip.background = element_blank()) +  
  scale_colour_manual(values=APApaletteSlim2,
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                      labels=c("Yoked to Trained", "Yoked to Conflict", "Trained", "Conflict")) + 
  scale_y_continuous(name="Total Entrances") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear, labeller = as_labeller(genoAPAnames)) 

## num shocks 
shocks <- behav %>% 
  filter(TrainSessionCombo %in% c("T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest"))  %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(!grepl("yoked", TrainGroup)) %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=NumShock, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.2, size=2)  +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  theme(strip.background = element_blank()) +  
  scale_colour_manual(values=APApaletteSlim3,
                      name="APA Training",
                      breaks=c("trained_trained", "trained_conflict"),
                      labels=c("Trained", "Conflict")) + 
  scale_y_continuous(name="Total Shocks") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels=c("1" = "T1", "2" = "T2", 
                              "3" = "T3", "4" = "Retest", "5" = "T4/C1",
                              "6" = "T5/C2", "7" = "T6/C3")) +
  facet_wrap(~genoYear, labeller = as_labeller(genoAPAnames)) 

## pTimeOpp - saved as beahv_pTimeOpp_6 or _3
pTimeOPP <- behav %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() %>%
  ggplot(aes(as.numeric(x=TrainSessionCombo), y=pTimeOPP, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.2, size=2)  +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  theme(strip.background = element_blank()) +  
  scale_colour_manual(values=APApaletteSlim2,
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                      labels=c("Yoked to Trained", "Yoked to Conflict", "Trained", "Conflict")) + 
  scale_y_continuous(name="pTimeOPP") + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_wrap(~genoYear, labeller = as_labeller(genoAPAnames)) 


## save plots ----
setwd("/Users/raynamharris/Github/BehavEphyRNAseq/results/2016_10_24_sfnposter/")
save_plot("one.pdf", one,
          base_aspect_ratio = 3 # make room for figure legend
)
save_plot("two.pdf", two, base_aspect_ratio = 3 )
save_plot("timespent.pdf", timespent, base_aspect_ratio = 3)
save_plot("entrances.pdf", entrances, base_aspect_ratio = 3)

save_plot("one.png", one, base_aspect_ratio = 3)
save_plot("two.png", two, base_aspect_ratio = 3 )


NumEntrancesTimeTarget <- plot_grid(entrances, timespent, ncol = 1)
save_plot("NumEntrancesTimeTarget.pdf", entrantime, 
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 3)

save_plot("NumEntrancesTimeTarget.png", entrantime, 
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 3)



## heatmap of data (not correlations, but raw/scaled data) ----
### melt to make long 



behav_long <- melt(behav, id=c("ID","APA","genoAPA","genoAPAsession","genoAPAsessionDay", "genoYear", "genoAPAsessionCombo",
                               "genoAPAsessionDayInd","TrainSessionCombo" ,"Genotype", "TrainSessionComboDay", "genoAPAyear",
                               "genoAPAsessionComboInd","TrainProtocol","TrainSequence","TrainGroup","Day","TrainSession",
                               "ShockOnOff","Year","PairedPartner","Experimenter",
                               "Housing","TestLocation","filename", "pair1", "pair2"))

behav_long <- filter(behav_long, !grepl("TotalTime.s|p.miss", variable )) %>% 
  filter(!grepl("16-357A|16-357B|16-357D", ID)) %>% 
  filter(TrainSessionCombo %in% c("Hab", "T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3", "Retest", "Retention"))  %>% 
  droplevels() 


#behav_long <- filter(behav_long, !grepl("Retention|Retest", genoAPAsession )) %>% droplevels()
behav_long$value <- as.numeric(behav_long$value)
str(behav_long)


##widen to get individual values by session
behav_long_genoAPAsessionComboInd <- dcast(behav_long, genoAPAsessionComboInd ~ variable, value.var= "value", fun.aggregate=mean)
names(behav_long_genoAPAsessionComboInd)
head(behav_long_genoAPAsessionComboInd)
behav_long_genoAPAsessionComboInd <- behav_long_genoAPAsessionComboInd[c(1:43)] #drop non-numeric colums 

##widen to get one animal per row!!
behav_long_oneanimal <- behav_long
behav_long_oneanimal$measure <- as.factor(paste(behav_long_oneanimal$variable, behav_long_oneanimal$TrainSessionCombo, sep="_"))
behav_long_oneanimalwide <- dcast(behav_long_oneanimal, ID ~ measure, value.var= "value", fun.aggregate=mean)
str(behav_long_oneanimalwide)
head(behav_long_oneanimalwide)

## now widen then lengthen to get group averages
behav_long_genoAPA <- dcast(behav_long, genoAPA ~ variable, value.var= "value", fun.aggregate=mean)
head(behav_long_genoAPA)

## scale columns
rownames(behav_long_genoAPA) <- behav_long_genoAPA$genoAPA    # set $genoAPAsession as rownames
behav_long_genoAPA <- behav_long_genoAPA[-c(4:5)]

behav_long_genoAPA[1] <- NULL
behav_long_genoAPA <- scale(behav_long_genoAPA)
head(behav_long_genoAPA)

## heatmap clusterd - saved as behav_heatmap.png
heatpalette <- colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 100)
heatmap.2(behav_long_genoAPA,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,25),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          RowSideColors = c("#7fbf7b", "#af8dc3", "#1b7837", "#762a83", "#00441b", "#40004b",
                            "#7fbf7b", "#af8dc3", "#1b7837", "#762a83", "#00441b", "#40004b"),
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 2, cexCol = 1.1)      # font size

## now widen then lengthen to get group averages by year (omit WT 2014)
behav_long_genoAPAyear <- dcast(behav_long, genoAPAyear ~ variable, value.var= "value", fun.aggregate=mean)
behav_long_genoAPAyear <- filter(behav_long_genoAPAyear, genoAPAyear != "WT_trained_2014")
behav_long_genoAPAyear$genoAPAyear

## scale columns
rownames(behav_long_genoAPAyear) <- behav_long_genoAPAyear$genoAPA    # set $genoAPAsession as rownames
behav_long_genoAPAyear <- behav_long_genoAPAyear[-c(4:5)]
behav_long_genoAPAyear[1] <- NULL
behav_long_genoAPAyear <- scale(behav_long_genoAPAyear)
head(behav_long_genoAPAyear)
str(behav_long_genoAPAyear)

## heatmap clusterd - saved as behav_heatmap.png
heatpalette <- colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 100)
heatmap.2(behav_long_genoAPAyear, 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,30),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          RowSideColors = c("#af8dc3", "#40004b", "#762a83", "#7fbf7b", 
                            "#00441b", "#1b7837", "#af8dc3", "#40004b",
                            "#40004b", "#762a83", "#762a83", "#7fbf7b",
                            "#00441b", "#00441b", "#1b7837", "#1b7837"),
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 2, cexCol = 1.1)      # font size


head(behav_long_oneanimalwide)
is.data.frame(behav_long_oneanimalwide)
rownames(behav_long_oneanimalwide) <- behav_long_oneanimalwide$ID    # set $genoAPAsession as rownames
behav_long_oneanimalwide[1] <- NULL
behav_long_oneanimalwide <- scale(behav_long_oneanimalwide)
head(behav_long_oneanimalwide)
behav_long_oneanimalwide <- na.omit(behav_long_oneanimalwide)

heatmap.2(behav_long_oneanimalwide, 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,15),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          RowSideColors = c("#af8dc3", "#40004b", "#762a83", "#7fbf7b", 
                            "#00441b", "#1b7837", "#af8dc3", "#40004b",
                            "#40004b", "#762a83", "#762a83", "#7fbf7b",
                            "#00441b", "#00441b", "#1b7837", "#1b7837"),
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 2, cexCol = 1.1)      # font size


## correlation matrix and plots ----
ggplot(behav_long_genoAPA, aes(x = genoAPA, y = variable, fill = value)) + 
  geom_tile() + coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## first, make the data a matrix with genoAPAsessionInd as the row names
behav_matrix <- behav   #create new dataframe
rownames(behav_matrix) <- behav_matrix$genoAPAsessionDayInd     # set $genoAPAsessionInd as rownames
names(behav_matrix)
behav_matrix <- behav_matrix[-c(1:32)] #delete all non-numeric columns and p.miss and TotalTime
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


