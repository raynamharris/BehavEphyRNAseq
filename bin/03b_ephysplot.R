

ephys_long$APA <- ifelse(grepl("train_trained", ephys_long$APA), "trained", 
                    ifelse(grepl("trained_conflict", ephys_long$APA), "trained",
                           ifelse(grepl("yoked_trained", ephys_long$APA), "yoked",
                                  ifelse(grepl("yoked_conflict", ephys_long$APA), "yoked", NA))))

ephys$APA <- ifelse(grepl("train_trained", ephys$APA), "trained", 
                         ifelse(grepl("trained_conflict", ephys$APA), "trained",
                                ifelse(grepl("yoked_trained", ephys$APA), "yoked",
                                       ifelse(grepl("yoked_conflict", ephys$APA), "yoked", NA))))

  

## plots  
ephys_long %>% 
  ggplot(aes(x=as.factor(variable), y=value, color=ID)) + 
  geom_point(size = 4) + 
  facet_grid(~APA + Genotype) +  
  scale_y_continuous(trans = "reverse")

ephys_long %>% 
  filter(Year == "2015") %>% 
  ggplot(aes(x=as.factor(variable), y=value, color=ID)) + 
  geom_point(size = 4) + 
  #facet_grid(~APA + Genotype) +  
  scale_y_continuous(trans = "reverse")









#### pca ----

W = ephys[,6:19]

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
qplot(scoresW[,1], scoresW[,2], color=ephys$APA, xlab='Component 1', ylab='Component 3')
qplot(scoresW[,1], scoresW[,2], color=ephys$Genotype, xlab='Component 1', ylab='Component 3')
qplot(scoresW[,1], scoresW[,3], color=ephys$APA, xlab='Component 1', ylab='Component 3')
qplot(scoresW[,2], scoresW[,3], color=ephys$APA, xlab='Component 1', ylab='Component 3')


## exporting the data
scoresWdf <- as.data.frame(scoresW)
scoresWdf$ID <-  ephys$ID
scoresWdf$APA <- ephys$APA

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
ggplot(scoresWdf, aes(PC1, PC2, colour=ephys$Genotype)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC1, PC2, colour=ephys$genoAPA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())   

ggplot(scoresWdf, aes(PC2, PC3, colour=ephys$genoAPA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC3, PC4, colour=ephys$genoAPA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC3, PC4, colour=ephys$Genotype)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  

ggplot(scoresWdf, aes(PC3, PC4, colour=ephys$APA)) + 
  geom_point(size = 8) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  theme(strip.background = element_blank())  


