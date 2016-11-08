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
    scale_colour_manual(values=c("#f1a340", "#7f3b08", "#9970ab", "#40004b"),
                        name="APA Training",
                        breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                        labels=c("Yoked to Trained", "Yoked to Conflict", "Trained Trained", "Trained Conflict")) 
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
    scale_colour_manual(values=c("#f1a340", "#7f3b08", "#9970ab", "#40004b"),
                        name="APA Training",
                        breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                        labels=c("Yoked to Trained", "Yoked to Conflict", "Trained Trained", "Trained Conflict")) 
    return(plot)
}

pcaplotgenotype(data=scoresdf, xcol="PC1", ycol="PC2", colorcode="APA")
pcaplotgenotype(data=scoresdf, xcol="PC1", ycol="PC5", colorcode="APA")

#save.image("~/Github/BehavEphyRNAseq/results/2016-11-07_sfnposter/Nov7.Rdata")