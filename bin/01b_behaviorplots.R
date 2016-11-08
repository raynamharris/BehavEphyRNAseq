setwd("~/Github/BehavEphyRNAseq/data/behavior/")
behav <- read.csv("behav.csv", header = T, check.names = F)
summary(behav)
## load libraries -----
library(ggdendro) ## for dendrograms!!
library(magrittr) ## to use the weird pipe
library(gplots) ##for making awesome plots
library(cowplot) ## for some easy to use themes
library(tidyr)
library(dplyr) ## for filtering and selecting rows
library(reshape2) #@ for melting dataframe
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!

behavMaddy <- behav %>% 
  filter(Experimenter %in% c("Maddy"))  %>%  droplevels() 
behavsofinterest <- c("Path1stEntr","Path2ndEntr", "pTimeOPP", "pTimeTarget")

## need to make TrainSessionComboNum for plotting purposes
levels(behavMaddy$TrainSessionCombo)
behavMaddy$TrainSessionComboNum <- ifelse(grepl("Hab", behavMaddy$TrainSessionCombo), "1", 
                                          ifelse(grepl("T1", behavMaddy$TrainSessionCombo), "2",
                                                 ifelse(grepl("T2", behavMaddy$TrainSessionCombo), "3",
                                                        ifelse(grepl("T3", behavMaddy$TrainSessionCombo), "4", 
                                                               ifelse(grepl("Retest", behavMaddy$TrainSessionCombo), "5",
                                                                      ifelse(grepl("T4_C1", behavMaddy$TrainSessionCombo), "6",
                                                                             ifelse(grepl("T5_C2", behavMaddy$TrainSessionCombo), "7",
                                                                                    ifelse(grepl("T6_C3", behavMaddy$TrainSessionCombo), "8",
                                                                                           ifelse(grepl("Retention", behavMaddy$TrainSessionCombo), "9",NA)))))))))
 
behavMaddy$TrainSessionComboNum <- as.numeric(as.character(behavMaddy$TrainSessionComboNum))
summary(behavMaddy$TrainSessionComboNum) 

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

RedBlack22 <- (values=c("#bababa","#1a1a1a","#f4a582","#b2182b"))

JMPalette <- c('black','grey50','red','darkorange')
WTPalette <- c('black','red')
FMR1Palette <- c('grey50','darkorange')

OrangePurple22 <- (values=c("#f1a340", "#7f3b08", "#9970ab", "#40004b"))

## ggplots across sessions!! -----

## creating better facet label names
genoAPAnames <- c(
  "WT_2015" = "WT 2015", 
  "WT_2016" = "WT 2016",
  "FMR1KO_2016" = "FMR1-KO 2016"
)

genonames <- c(
  "WT" = "WT", 
  "FMR1KO" = "FMR1-KO"
)


# single behavior plot
behavMaddy %>% 
  ggplot(aes(x=TrainSessionComboNum, y=Path2ndEntr, color=APA)) + 
  geom_point(size=2) + geom_jitter() +
  stat_smooth(alpha=0.2, size=2)  +
  theme_cowplot(font_size = 25, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
  #theme(strip.background = element_blank()) +  
  scale_colour_manual(values=c("#f1a340", "#7f3b08", "#9970ab", "#40004b"),
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                      labels=c("Yoked to Trained", "Yoked to Conflict", "Trained Trained", "Trained Conflict")) + 
  scale_y_continuous(name="Path to 2nd Entrance (m)", limits = c(0, 17)) + 
  scale_x_continuous(name =NULL, 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                              "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                              "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
  facet_grid(~Genotype, labeller = as_labeller(genonames))

## behavior plot function!!! with legend
behaviorplot <- function(data, xcol, ycol, yaxislabel, colorcode){
  print(ycol)
  print(yaxislabel)
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=2, aes(shape = Genotype)) + geom_jitter() +
    stat_smooth(alpha=0.2, size=2)  +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
    scale_colour_manual(values=c("#f1a340", "#7f3b08", "#9970ab", "#40004b"),
                        name="APA Training",
                        breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                        labels=c("Yoked to Trained", "Yoked to Conflict", "Trained Trained", "Trained Conflict")) + 
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name =NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
    facet_grid(~Genotype, labeller = as_labeller(genonames))
  print(plot)
  return(plot)
}

A4 <- behaviorplot(data=behavMaddy, xcol="TrainSessionComboNum", ycol="Path2ndEntr", yaxislabel="Path to 2nd Entrance (m)", colorcode="APA")
B4 <- behaviorplot(data=behavMaddy, xcol="TrainSessionComboNum", ycol="Speed2ndEntr", yaxislabel='Speed to 2nd Entrance (cm/s)' , colorcode="APA")


## behavior plot function without legend
behaviorplotnolegend <- function(data, xcol, ycol, yaxislabel, colorcode){
  print(ycol)
  print(yaxislabel)
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=2, aes(shape = Genotype)) + geom_jitter() +
    stat_smooth(alpha=0.2, size=2)  +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    theme(axis.text.x = element_text(angle=70, vjust=0.5)) +
    scale_colour_manual(values=c("#f1a340", "#7f3b08", "#9970ab", "#40004b"),
                        name="APA Training",
                        breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                        labels=c("Yoked to Trained", "Yoked to Conflict", "Trained Trained", "Trained Conflict")) +
    theme(legend.position="none") +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name =NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
    facet_grid(~Genotype, labeller = as_labeller(genonames))
  print(plot)
  return(plot)
}

names(behavMaddy)
A3 <- behaviorplotnolegend(data=behavMaddy, xcol="TrainSessionComboNum", ycol="Path1stEntr", yaxislabel="Path to 1st Entrance (m)", colorcode="APA")
A2 <- behaviorplotnolegend(data=behavMaddy, xcol="TrainSessionComboNum", ycol="pTimeTarget", yaxislabel='Time in the "Shock Zone" (s)' , colorcode="APA")
B1 <- behaviorplotnolegend(data=behavMaddy, xcol="TrainSessionComboNum", ycol="RayleigLength", yaxislabel='Rayleig Length' , colorcode="APA")
B3 <- behaviorplotnolegend(data=behavMaddy, xcol="TrainSessionComboNum", ycol="Speed1stEntr.cm.s.", yaxislabel='Speed to 1st Entrance (cm/s)' , colorcode="APA")
A1 <- behaviorplotnolegend(data=behavMaddy, xcol="TrainSessionComboNum", ycol="TotalPathArena.m", yaxislabel='Total Path (m)' , colorcode="APA")
B2 <- behaviorplotnolegend(data=behavMaddy, xcol="TrainSessionComboNum", ycol="AnnularMaxVal", yaxislabel='Maximum Annular Value' , colorcode="APA")

## make a grid plot

plot_grid(A1, A2, A3, B1, B2, B3, nrow = 2)
plot_grid(A4, B4, nrow = 2)


## save plots ----
setwd("/Users/raynamharris/Github/BehavEphyRNAseq/results/2016-11-01_sfnposter/")
save_plot("timespent.pdf", timespent, base_aspect_ratio = 3)
save_plot("entrances.pdf", entrances, base_aspect_ratio = 3)

save_plot("one.png", one, base_aspect_ratio = 3)
save_plot("two.png", two, base_aspect_ratio = 3 )
save_plot("pTimeTarget.png", pTimeTarget, base_aspect_ratio = 3 )
save_plot("pTimeOPP.png", pTimeOPP, base_aspect_ratio = 3 )
save_plot("RayleigLength.png", RayleigLength, base_aspect_ratio = 3 )





## heatmap of data (not correlations, but raw/scaled data) ----

#melt the data
behav_long  <- melt(behav, id = c(1:21))
behav_long <- behav_long %>% drop_na()
## now widen then lengthen to get group averages
behav_long_genoAPA <- dcast(behav_long, genoAPA ~ variable, value.var= "value", fun.aggregate=mean)
head(behav_long_genoAPA)
## scale columns
rownames(behav_long_genoAPA) <- behav_long_genoAPA$genoAPA    # set $genoAPAsession as rownames
behav_long_genoAPA[1] <- NULL
behav_long_genoAPA <- scale(behav_long_genoAPA)
head(behav_long_genoAPA)
## heatmap clusterd - saved as behav_heatmap.png
heatpalette <- colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(n = 100)
heatmap.2(behav_long_genoAPA, 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,18),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          #RowSideColors = c("#7fbf7b", "#af8dc3", "#1b7837", "#762a83", "#00441b", "#40004b",
          #                  "#7fbf7b", "#af8dc3", "#1b7837", "#762a83", "#00441b", "#40004b"),
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 1.5, cexCol = 1)      # font size


#melt the data
behav_long_year  <- melt(behav, id = c(1:21))
behav_long_year <- behav_long_year %>% drop_na()
## now widen then lengthen to get group averages
behav_long_year$genoAPAyear <- as.factor(paste(behav_long_year$APA,behav_long_year$Genotype, behav_long_year$Year, sep="_"))
behav_long_year <- dcast(behav_long_year, genoAPAyear ~ variable, value.var= "value", fun.aggregate=mean)
head(behav_long_year)
## scale columns
rownames(behav_long_year) <- behav_long_year$genoAPA    # set $genoAPAsession as rownames
behav_long_year[1] <- NULL
behav_long_year <- scale(behav_long_year)
head(behav_long_year)
## heatmap clusterd - saved as behav_heatmap.png
heatpalette <- colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(n = 100)
heatmap.2(behav_long_year, 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,22),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          RowSideColors = c("#40004b", "#40004b", "#40004b", "#9970ab", "#9970ab", "#9970ab",
                            "#b35806", "#b35806", "#b35806", "#f1a340", "#f1a340", "#f1a340"),
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 1.5, cexCol = 1)      # font size


## for APApalette2

#40004b purple train conflict
#9970ab purple train train
#c2a5cf purpletrain 

#fdb863 orange untrained
#f1a340 orange untrained trained
#b35806 orange untrained conflict



## behavbysessionHabT3
behavbysessionHabT3 <- behav %>% filter(Experimenter %in% c("Maddy")) %>%
                                          filter(TrainSessionCombo %in% c("T1","T2","T3","Hab")) %>%  droplevels()

behavbysessionHabT3 <- melt(behavbysessionHabT3, id = c(1:21))
behavbysessionHabT3 <- behavbysessionHabT3 %>% drop_na()
## create the bysession column and wide
names(behavbysessionHabT3)
head(behavbysessionHabT3)
behavbysessionHabT3$bysession <- as.factor(paste(behavbysessionHabT3$ShockOnOff, behavbysessionHabT3$TrainSessionCombo, behavbysessionHabT3$Genotype, behavbysessionHabT3$TrainGroup, sep="_"))
behavbysessionHabT3 <- dcast(behavbysessionHabT3, bysession   ~ variable, value.var= "value", fun.aggregate = mean)
summary(behavbysessionHabT3) 
head(behavbysessionHabT3)
## makin' it a matrix
is.data.frame(behavbysessionHabT3)
rownames(behavbysessionHabT3) <- behavbysessionHabT3$bysession    # set $genoAPAsession as rownames
names(behavbysessionHabT3)
behavbysessionHabT3 <- behavbysessionHabT3[-c(1)]
## send the factors to 
behavbysessionHabT3 <- scale(behavbysessionHabT3)
head(behavbysessionHabT3)
behavbysessionHabT3 <- na.omit(behavbysessionHabT3)
#dev.off()
heatmap.2(behavbysessionHabT3, 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 0.8, cexCol = 0.6)      # font size
# saved as behavbysessionHabT3.png


## behavbysessionT1C4T5T6
behavbysessionT1C4T5T6 <- melt(behav, id = c(1:21))
behavbysessionT1C4T5T6 <- behavbysessionT1C4T5T6 %>% 
  filter(Experimenter %in% c("Maddy")) %>%
  filter(Year == 2015)  %>%
  filter(TrainSessionCombo %in% c("T1","T2","T3","T4_C1", "T5_C2", "T6_C3"))  %>%  droplevels() 
## create the bysession column and wide
names(behavbysessionT1C4T5T6)
head(behavbysessionT1C4T5T6)
behavbysessionT1C4T5T6$bysession <- as.factor(paste(behavbysessionT1C4T5T6$TrainSessionCombo, behavbysessionT1C4T5T6$Genotype, behavbysessionT1C4T5T6$TrainSequence, sep="_"))
behavbysessionT1C4T5T6 <- dcast(behavbysessionT1C4T5T6, bysession   ~ variable, value.var= "value", fun.aggregate = mean)
summary(behavbysessionT1C4T5T6) 
head(behavbysessionT1C4T5T6)
## makin' it a matrix
is.data.frame(behavbysessionT1C4T5T6)
rownames(behavbysessionT1C4T5T6) <- behavbysessionT1C4T5T6$bysession    # set $genoAPAsession as rownames
names(behavbysessionT1C4T5T6)
behavbysessionT1C4T5T6 <- behavbysessionT1C4T5T6[-c(1)]
## send the factors to 
behavbysessionT1C4T5T6 <- scale(behavbysessionT1C4T5T6)
head(behavbysessionT1C4T5T6)
behavbysessionT1C4T5T6 <- na.omit(behavbysessionT1C4T5T6)
#dev.off()
heatmap.2(behavbysessionT1C4T5T6, 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          col=heatpalette,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          srtCol=45,  adjCol = c(1,1), # angled column label
          cexRow = 0.8, cexCol = 0.6)      # font size
# saved as behavbysessionT1C4T5T6.png

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


