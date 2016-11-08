library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(reshape2) #@ for melting dataframe
# Create the dataset WT trained only 

forpca <- melt(behav, id = c(1:21))
forpca <- forpca %>% 
  filter(Experimenter == "Maddy") %>%  
  #filter(!grepl("Bin", variable)) %>%  
  droplevels()
forpca$value <- as.numeric(forpca$value)
forpca$bysession <- as.factor(paste(forpca$TrainSessionCombo, forpca$variable, sep="_"))
forpca <- dcast(forpca, ID + APA + Genotype + TrainProtocol + TrainSequence + TrainGroup + 
                          PairedPartner + Experimenter + Housing + TestLocation + genoYear + APA + 
                          genoAPA + pair1 + pair2 ~ bysession, value.var= "value", fun.aggregate = mean)
summary(forpca)
nacols<-which(colSums(is.na(forpca))>=1)   ## removing cols with NA
forpca<-forpca[,-c(nacols)]

# Pick out the numerical columns from the data set (run ?iris for further info)
Z = forpca[,16:171]

Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
summary(Z)

# Run PCA
pc = prcomp(Z, scale.=TRUE)

# Look at the summary methods
pc
summary(pc)
head(pc)

behav2015pca <-  as.data.frame(pc$x)

#write.csv(pc$rotation, "pca.csv")

## a bunch of different plots
plot(pc)
#biplot(pc)


### quick pca plot
pc$rotation
loadings = pc$rotation
scores = pc$x
head(scores)
qplot(scores[,1], scores[,2], color=forpca$APA, xlab='Component 1', ylab='Component 2', size=2)
qplot(scores[,1], scores[,3], color=forpca$APA, xlab='Component 1', ylab='Component 3', size=2)
qplot(scores[,2], scores[,3], color=forpca$APA, xlab='Component 2', ylab='Component 3', size=2)

qplot(scores[,1], scores[,2], color=forpca$Genotype, xlab='Component 1', ylab='Component 2', size=2)
qplot(scores[,1], scores[,3], color=forpca$Genotype, xlab='Component 1', ylab='Component 3', size=2)
qplot(scores[,2], scores[,3], color=forpca$Genotype, xlab='Component 2', ylab='Component 3', size=2)



## exporting the data
scoresdf <- as.data.frame(scores)
scoresdf$ID <-  forpca$ID
scoresdf$APA <- forpca$APA
scoresdf$Genotype <- forpca$Genotype
scoresdf$TrainGroup <- forpca$TrainGroup
scoresdf$genoYear <- forpca$genoYear


### ggplotified
ggplot(scoresdf, aes(PC1, PC2, colour=APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  #theme(strip.background = element_blank()) +
  facet_grid(~Genotype, labeller = as_labeller(genonames)) + 
  scale_colour_manual(name="APA Training",
                      values=c("#f1a340", "#7f3b08", "#9970ab","#40004b")) 


genoAPAnames <- c(
  "WT_2015" = "WT 2015", 
  "WT_2016" = "WT 2016",
  "FMR1KO_2016" = "FMR1-KO 2016"
)


## pca plot function!!! 
pcaplotgenotypeyear <- function(data, xcol, ycol, colorcode){
  print(ycol)
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size = 8, aes(shape = Genotype)) +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    #theme(strip.background = element_blank()) +
    facet_grid(~genoYear, labeller = as_labeller(genoAPAnames)) + 
    scale_colour_manual(name="APA Training",
                        values=c("#f1a340", "#7f3b08", "#9970ab","#40004b")) 
  return(plot)
}

pcaplotgenotypeyear(data=scoresdf, xcol="PC1", ycol="PC2", colorcode="APA")
pcaplotgenotypeyear(data=scoresdf, xcol="PC1", ycol="PC5", colorcode="APA")


## pca plot function!!! 
pcaplotgenotype <- function(data, xcol, ycol, colorcode){
  print(ycol)
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size = 8, aes(shape = Genotype)) +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    #theme(strip.background = element_blank()) +
    scale_colour_manual(name="APA Training",
                        values=c("#f1a340", "#7f3b08", "#9970ab","#40004b")) 
  return(plot)
}

pcaplotgenotype(data=scoresdf, xcol="PC1", ycol="PC2", colorcode="APA")
pcaplotgenotype(data=scoresdf, xcol="PC1", ycol="PC5", colorcode="APA")





## other pca plot found online http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
#library("devtools")
#install_github("kassambara/factoextra")
library("factoextra")

res.pca <- prcomp(Z,  scale = TRUE)
fviz_pca_biplot(res.pca, habillage = forpca$APA, addEllipses = FALSE, pointshape = 19, pointsize = 5,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                select.var = list(contrib = 6)) +
                scale_colour_manual(values=c("#9970ab","#40004b"), #762a83 train #40004b conflict
                                      name="APA Training",
                                      breaks=c("trained_trained", "trained_conflict"),
                                      labels=c("Trained", "Trained with Conflict")) 

fviz_pca_biplot(res.pca, geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5, 
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 4), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("trained_trained", "trained_conflict"),
                      labels=c("Trained", "Trained with Conflict"))  + 
  scale_x_continuous(limits = c(-80, 20))


fviz_pca_biplot(pc, axes = c(1, 2), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 2), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("trained_trained", "trained_conflict"),
                      labels=c("Trained", "Trained with Conflict")) + 
  scale_x_continuous(limits = c(-20, 20))














# Create the dataset WT trained only 

forpca <- melt(behav, id = c(1:21))
forpca <- forpca %>% 
  filter(genoYear == "WT_2015") %>%  
  filter(!grepl("Bin|Shock", variable)) %>%  
  filter(!grepl("Retest|Hab|Reten", TrainSessionCombo)) %>%  
  droplevels()
forpca$value <- as.numeric(forpca$value)
forpca$bysession <- as.factor(paste(forpca$TrainSessionCombo, forpca$variable, sep="_"))
forpca <- dcast(forpca, ID + APA + Genotype + TrainProtocol + TrainSequence + TrainGroup + 
                  PairedPartner + Experimenter + Housing + TestLocation + genoYear + APA + 
                  genoAPA + pair1 + pair2 ~ bysession, value.var= "value", fun.aggregate = mean)
summary(forpca) 
names(forpca)


# Pick out the numerical columns from the data set (run ?iris for further info)
Z = forpca[,16:207]

Z <- Z[!colSums(Z)%in%NA]   ## removing cols with NA
Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]

# Run PCA
pc = prcomp(Z, scale.=TRUE)

# Look at the summary methods
pc
summary(pc)


## a bunch of different plots
plot(pc)
#biplot(pc)

### quick pca plot
pc$rotation
loadings = pc$rotation
scores = pc$x
qplot(scores[,1], scores[,3], color=forpca$APA, xlab='Component 1', ylab='Component 3') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))
qplot(scores[,1], scores[,2], color=forpca$APA, xlab='Component 1', ylab='Component 2') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))
qplot(scores[,2], scores[,3], color=forpca$APA, xlab='Component 2', ylab='Component 3') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))
qplot(scores[,4], scores[,3], color=forpca$APA, xlab='Component 2', ylab='Component 3') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))


## exporting the data
scoresdf <- as.data.frame(scores)
scoresdf$ID <-  forpca$ID
scoresdf$APA <- forpca$APA


### ggplotified
ggplot(scoresdf, aes(PC1, PC3, colour=forpca$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank()) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict","trained_trained", "trained_conflict"),
                      labels=c("Yoked", "Yoked with Conflict", "Trained", "Trained with Conflict")
  ) 

ggplot(scoresdf, aes(PC1, PC2, colour=forpca$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank()) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict","trained_trained", "trained_conflict"),
                      labels=c("Yoked", "Yoked with Conflict", "Trained", "Trained with Conflict")
  ) 




# now plot, using geom_segment() for arrows and geom_text for labels
rotation_data <- data.frame(pc$rotation, variable=row.names(pc$rotation))
rotation_data <- subset(rotation_data, PC1 > 0.119 | PC1 < -0.12, select=c(PC1, PC2))
PC2 <- subset(rotation_data, PC2 > 0.112 | PC2 < -0.18, select=c(PC2))


arrow_style <- arrow(length = unit(0.05, "inches"),
                     type = "closed")

ggplot(rotation_data) + 
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow_style) + 
  geom_text(aes(x=PC1, y=PC2, label=variable), hjust=0, size=3, color='red') + 
  xlim(-1.,1.25) + 
  ylim(-1.,1.) +
  coord_fixed()


## other pca plot found online http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
#library("devtools")
#install_github("kassambara/factoextra")
library("factoextra")

fviz_pca_biplot(pc, axes = c(1, 2),  geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 4), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) + 
  scale_x_continuous(limits = c(-18, 18)) 


fviz_pca_biplot(pc, axes = c(1, 3), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 4), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) + 
  scale_x_continuous(limits = c(-18, 18))



























# Create the dataset WT trained only 

forpca <- melt(behav, id = c(1:21))
forpca <- forpca %>% 
  filter(genoYear %in% c("WT_2015", "FMR1KO_2016")) %>%  
  filter(!grepl("Bin|Shock", variable)) %>%  
  filter(!grepl("Retest|Hab|Reten", TrainSessionCombo)) %>%  
    droplevels()
  
forpca$value <- as.numeric(forpca$value)

forpca$bysession <- as.factor(paste(forpca$TrainSessionCombo, forpca$variable, sep="_"))
forpca <- dcast(forpca, ID + APA + Genotype + TrainProtocol + TrainSequence + TrainGroup + 
                  PairedPartner + Experimenter + Housing + TestLocation + genoYear + APA + 
                  genoAPA + pair1 + pair2 ~ bysession, value.var= "value", fun.aggregate = mean)
summary(forpca) 
names(forpca)


# Pick out the numerical columns from the data set (run ?iris for further info)
Z = forpca[,16:207]

Z <- Z[!colSums(Z)%in%NA]   ## removing cols with NA
Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]

# Run PCA
pc = prcomp(Z, scale.=TRUE)

# Look at the summary methods
pc
summary(pc)

behavWT15pca <-  pc$rotation

#write.csv(pc$rotation, "pca.csv")

## a bunch of different plots
plot(pc)
#biplot(pc)

### quick pca plot
pc$rotation
loadings = pc$rotation
scores = pc$x
qplot(scores[,1], scores[,3], color=forpca$APA, xlab='Component 1', ylab='Component 3') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))
qplot(scores[,1], scores[,2], color=forpca$APA, xlab='Component 1', ylab='Component 2') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))
qplot(scores[,2], scores[,3], color=forpca$APA, xlab='Component 2', ylab='Component 3') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))
qplot(scores[,4], scores[,3], color=forpca$APA, xlab='Component 2', ylab='Component 3') + 
  geom_point(size = 8) + scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"))


## exporting the data
scoresdf <- as.data.frame(scores)
scoresdf$ID <-  forpca$ID
scoresdf$APA <- forpca$APA

### ggplotified
ggplot(scoresdf, aes(PC1, PC3, colour=forpca$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank()) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict","trained_trained", "trained_conflict"),
                      labels=c("Yoked", "Yoked with Conflict", "Trained", "Trained with Conflict")
  ) 

ggplot(scoresdf, aes(PC1, PC2, colour=forpca$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank()) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("yoked_trained", "yoked_conflict","trained_trained", "trained_conflict"),
                      labels=c("Yoked", "Yoked with Conflict", "Trained", "Trained with Conflict")
  ) 




# now plot, using geom_segment() for arrows and geom_text for labels
rotation_data <- data.frame(pc$rotation, variable=row.names(pc$rotation))
rotation_data <- subset(rotation_data, PC1 > 0.119 | PC1 < -0.12, select=c(PC1, PC2,varoable))
PC2 <- subset(rotation_data, PC2 > 0.112 | PC2 < -0.18, select=c(PC2))


arrow_style <- arrow(length = unit(0.05, "inches"),
                     type = "closed")

ggplot(rotation_data) + 
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow_style) + 
  geom_text(aes(x=PC1, y=PC2, label=variable), hjust=0, size=3, color='red') + 
  xlim(-1.,1.25) + 
  ylim(-1.,1.) +
  coord_fixed()


## other pca plot found online http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
#library("devtools")
#install_github("kassambara/factoextra")
library("factoextra")

p12 <- fviz_pca_biplot(pc, axes = c(1, 2),  geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 5), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) + 
  scale_x_continuous(limits = c(-18, 18)) 


p34 <- fviz_pca_biplot(pc, axes = c(3, 4), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 5), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) 

p56 <- fviz_pca_biplot(pc, axes = c(5, 6), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                select.var = list(contrib = 5), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) 


p78 <- fviz_pca_biplot(pc, axes = c(7, 8), geom = c("point"),
                       label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                       habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                       col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                       col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                       #col.circle = "grey70", 
                       select.var = list(contrib = 5), 
                       select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                       jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#f1a340","#b35806","#9970ab","#40004b")) 

setwd("~/Github/BehavEphyRNAseq/results/2016-11-01_sfnposter")
save_plot("p12.png", p12, base_aspect_ratio = 2)
save_plot("p34.png", p34, base_aspect_ratio = 2)
save_plot("p56.png", p56, base_aspect_ratio = 2)
save_plot("p78.png", p78, base_aspect_ratio = 2)
