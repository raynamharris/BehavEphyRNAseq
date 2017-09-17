source("figureoptions.R")


## Stimulus strength on x fepsp slope on y
ephys2_long %>% 
  ggplot(aes(x=variablenumeric, y=value, color=APA )) + 
  stat_smooth(alpha=0.2, size=2) +
  geom_jitter(size=2, width = 0.5) +
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  facet_grid(~Genotype) +   
  scale_y_continuous(trans = "reverse") + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("0", "10", "15", "20", "30", "40", "50")) +
  scale_color_manual(values = colorvalAPA) + 
  labs(x = "Stimulus Strenght (V)", y = "fEPSP Slope") + 
  theme(legend.position="none")


ephys2 %>% 
  ggplot(aes(x=APA, y=min, fill=APA)) +
  geom_boxplot() +
  scale_y_continuous(trans = "reverse") + 
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  scale_fill_manual(values = colorvalAPA) +
  labs(y = "Maximum fEPSP Slope", x = NULL) +
  facet_grid(~Genotype) 
  

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


