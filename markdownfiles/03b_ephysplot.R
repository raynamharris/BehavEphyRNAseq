source("figureoptions.R")
levels(ephys2$APA2) <- c("yoked-consistent","yoked-conflict","consistent", "conflict")
levels(ephys2_long$APA2) <- c("yoked-consistent","yoked-conflict","consistent", "conflict")

head(ephys2)
ephys2summaryNum <- dplyr::summarise(group_by(ephys2, Genotype, APA2, genoAPA), m = mean(min), se = sd(min)/sqrt(length(min)), len = length(min))
ephys2summaryNum <- as.data.frame(ephys2summaryNum)
levels(ephys2summaryNum$Genotype) <- c("WT","FMR1KO")
levels(ephys2summaryNum$APA2) <- c("yoked-consistent","yoked-conflict", "consistent", "conflict")


ephys2$APA2 <- factor(ephys2$APA2, levels = c("yoked-consistent",  "consistent", "yoked-conflict", "conflict"))

## Stimulus strength on x fepsp slope on y
ephys2_long %>% 
  ggplot(aes(x=variablenumeric, y=value, color=APA2 )) + 
  stat_smooth(alpha=0.2, size=2) +
  geom_jitter(size=2, width = 0.5) +
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  facet_grid(~Genotype) +   
  scale_y_continuous(trans = "reverse",
                     limits = c(0.005, -0.010)) + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("0", "10", "15", "20", "30", "40", "50")) +
  scale_color_manual(values = colorvalAPA2) + 
  labs(x = "Stimulus Strenght (V)", y = "fEPSP Slope") + 
  theme(legend.position="none")

ephys2_long %>% 
  ggplot(aes(x=variablenumeric, y=value, color=Genotype )) + 
  stat_smooth(alpha=0.2, size=2) +
  geom_jitter(size=2, width = 0.5) +
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  facet_wrap(~APA2,nrow=2) +   
  scale_y_continuous(trans = "reverse") + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("0", "10", "15", "20", "30", "40", "50")) +
  scale_color_manual(values = colorvalGenotype) + 
  labs(x = "Stimulus Strenght (V)", y = "fEPSP Slope") + 
  theme(legend.position="none")

ephys2 %>% 
  ggplot(aes(x=APA2, y=min, fill=APA2)) +
  geom_boxplot() +
  scale_y_continuous(trans = "reverse") + 
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  scale_fill_manual(values = colorvalAPA2) +
  labs(y = "Maximum fEPSP Slope", x = NULL) +
  facet_grid(~Genotype)

myggplot <- ephys2 %>% 
  ggplot(aes(x=Genotype, y=min, fill=Genotype)) +
  geom_boxplot() +
  scale_y_continuous(trans = "reverse") + 
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 8, line_size = 1) + 
  scale_fill_manual(values = colorvalGenotype) +
  labs(y = "Maximum fEPSP Slope", x = NULL) +
  facet_grid(~APA2)  +
  theme(legend.position="none")
myggplot

pdf(file="~/Github/FMR1rnaseqCA1/figures/04_ephys/ephys1.pdf", width=3.8, height=2.5)
plot(myggplot)
dev.off()


myggplot <- ephys2 %>% 
  ggplot(aes(x=Genotype, y=min, fill=Genotype)) +
  geom_boxplot() +
  scale_y_continuous(trans = "reverse") + 
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 8, line_size = 0.5) + 
  scale_fill_manual(values = colorvalGenotype) +
  labs(y = "Maximum fEPSP Slope", x = NULL) +
  facet_wrap(~APA2, nrow=2) +
  theme(legend.position="none")

pdf(file="~/Github/FMR1rnaseqCA1/figures/04_ephys/ephys2.pdf", width=2.5, height=2.5)
plot(myggplot)
dev.off()


aov1 <- aov(min ~ Genotype, data=ephys2)
summary(aov1) 
aov1 <- aov(min ~ Genotype * APA2, data=ephys2)
summary(aov1) 

aov1 <- aov(max ~ Genotype, data=ephys2)
summary(aov1) 
aov1 <- aov(max ~ Genotype * APA2, data=ephys2)
summary(aov1)

summary(ephys2$min)
  

#### pca ----

W = ephys2[,6:19]

W <- W[!colSums(W)%in%NA]   ## removing cols with NA
W <- W[,apply(W, 2, var, na.rm=TRUE) != 0]

# Run PCA
W
pcW = prcomp(W, scale.=TRUE)

# Look at the summary methods
pcW
summary(pcW)

## a bunch of different plots
plot(pcW)
#biplot(pc)


### quick pca plot
pcW$rotation
loadingsW = pcW$rotation
scoresW = pcW$x
qplot(scoresW[,1], scoresW[,2], color=ephys2$APA, xlab='Component 1', ylab='Component 3')
qplot(scoresW[,1], scoresW[,2], color=ephys2$Genotype, xlab='Component 1', ylab='Component 3')
qplot(scoresW[,1], scoresW[,3], color=ephys2$APA, xlab='Component 1', ylab='Component 3')
qplot(scoresW[,2], scoresW[,3], color=ephys2$APA, xlab='Component 1', ylab='Component 3')


## exporting the data
scoresWdf <- as.data.frame(scoresW)
scoresWdf$ID <-  ephys2$ID
scoresWdf$APA <- ephys2$APA

# capture the rotation matrix in a data frame
rotation_data <- data.frame(pcW$rotation, variable=row.names(pcW$rotation))
# define a pleasing arrow style
arrow_style <- arrow(length = unit(0.05, "inches"),
                     type = "closed")

ggplot(rotation_data) + 
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow_style) + 
  geom_text(aes(x=PC1, y=PC2, label=variable), hjust=0, size=5, color='red') 

ggplot(rotation_data) + 
  geom_segment(aes(xend=PC2, yend=PC3), x=0, y=0, arrow=arrow_style) + 
  geom_text(aes(x=PC2, y=PC3, label=variable), hjust=0, size=5, color='red') 



### ggplotified
ggplot(scoresWdf, aes(PC1, PC2, colour=ephys2$Genotype)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC1, PC2, colour=ephys2$genoAPA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())   

ggplot(scoresWdf, aes(PC2, PC3, colour=ephys2$genoAPA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC3, PC4, colour=ephys2$genoAPA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC3, PC4, colour=ephys2$Genotype)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC3, PC4, colour=ephys2$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  


