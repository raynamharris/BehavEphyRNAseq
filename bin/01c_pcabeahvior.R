library(ggplot2)
library(cowplot)
library(dplyr)

# Create the dataset
forpca <- melt(behav, id=c("ID","APA","genoAPA","genoAPAsession","genoAPAsessionDay", "genoYear", "genoAPAsessionCombo",
                                   "genoAPAsessionDayInd","TrainSessionCombo" ,"Genotype", "TrainSessionComboDay", "genoAPAyear",
                                   "genoAPAsessionComboInd","TrainProtocol","TrainSequence","TrainGroup","Day","TrainSession",
                                   "ShockOnOff","Year","PairedPartner","Experimenter",
                                   "Housing","TestLocation", "filename", "pair1", "pair2"))
forpca <- filter(forpca, !grepl("TotalTime.s|p.miss", variable )) %>% 
  filter(!grepl("15141A|15141B", ID)) %>%  # removed two animals with no retention values
  filter(genoYear == "WT_2015", TrainGroup == "trained") %>% 
  filter(TrainSessionCombo %in% c("T1","T2","T3","T4_C1", 
                                  "T5_C2", "T6_C3"))  %>%  droplevels() 
## create the bysession column and wide
forpca$bysession <- as.factor(paste(forpca$TrainSessionCombo, forpca$variable, sep="_"))
forpca$value <- as.numeric(forpca$value)
forpca <- dcast(forpca, ID + APA + Genotype + TrainSequence + TrainGroup ~ bysession, value.var= "value", fun.aggregate = mean)
summary(forpca) 
head(forpca)
names(forpca)


# Pick out the numerical columns from the data set (run ?iris for further info)
Z = forpca[,6:257]

# Run PCA
pc = prcomp(Z, scale.=TRUE)

# Look at the summary methods
pc
summary(pc)

behavWT15pca <-  pc$rotation

#write.csv(pc$rotation, "pca.csv")

## a bunch of different plots
plot(pc)
biplot(pc)


### quick pca plot
pc$rotation
loadings = pc$rotation
scores = pc$x
qplot(scores[,1], scores[,3], color=forpca$APA, xlab='Component 1', ylab='Component 3', size=2)


## exporting the data
scoresdf <- as.data.frame(scores)
scoresdf$ID <-  forpca$ID
scoresdf$APA <- forpca$APA

### ggplotified
ggplot(scoresdf, aes(PC2, PC3, colour=forpca$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank()) +
  scale_colour_manual(values=c("#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("trained_trained", "trained_conflict"),
                      labels=c("Trained", "Trained with Conflict"))   


## other pca plot found online http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
library("devtools")
install_github("kassambara/factoextra")
library("factoextra")

fviz_pca_biplot(pc, axes = c(2, 3), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 5,
                habillage = forpca$APA , addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 3,
                col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
                #col.circle = "grey70", 
                #select.var = list(contrib = 3), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_colour_manual(values=c("#9970ab","#40004b"), #762a83 train #40004b conflict
                      name="APA Training",
                      breaks=c("trained_trained", "trained_conflict"),
                      labels=c("Trained", "Trained with Conflict"))  

  
#wilke cow plot way
rotation_data <- data.frame(pc$rotation, variable=row.names(pc$rotation))
arrow_style <- arrow(length = unit(0.05, "inches"),
                     type = "closed")
ggplot(rotation_data) + 
  geom_segment(aes(xend=rotation_data$PC2, yend=rotation_data$PC3), x=0, y=0, arrow=arrow_style) + 
  geom_text(aes(x=PC1, y=PC2, label=variable), hjust=0, size=3, color='red')

## Rpubs 
install.packages("ggfortify")
library(ggfortify)
df <- forpca[,6:257]
autoplot(prcomp(df))
autoplot(prcomp(df), data = forpca, colour = 'APA', 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


