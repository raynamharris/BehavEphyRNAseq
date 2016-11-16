
## make long and tidy
ephys_long <- melt(ephysnoNA, id = c(1:5))
iomax_long <- melt(iomax, id = c(12:14))
ephysnoNA_long <- melt(ephysnoNA, id = c(1:5))
summary(ephysnoNA_long)
head(ephys_long$APA)


## add numeric for stat smooth
ephys_long$variablenumeric <- ifelse((ephys_long$variable == "V0"), "1", 
                                     ifelse(grepl("V10", ephys_long$variable ), "2",
                                            ifelse(grepl("V15", ephys_long$variable ), "3",
                                                   ifelse(grepl("V20", ephys_long$variable), "4", 
                                                          ifelse(grepl("V30", ephys_long$variable), "5",
                                                                 ifelse(grepl("V40", ephys_long$variable), "6",
                                                                        ifelse(grepl("V50", ephys_long$variable), "7", NA)))))))
ephys_long <- ephys_long %>% drop_na()

levels(ephys_long$variable)
ephys_long$variablenumeric <- as.numeric(ephys_long$variablenumeric)

ephys_long$Genotype <- factor(ephys_long$Genotype, 
                    levels = c("WT", "FMR1KO"))
ephys_long$APA <- factor(ephys_long$APA, 
                              levels = c("Yoked", "Trained", "Conflict"))

levels(behav$APA)


## plots  
iobygeneo <- ephys_long %>% 
  ggplot(aes(x=variablenumeric, y=value, color=APA )) + 
  stat_smooth(alpha=0.2, size=2) +
  geom_jitter(size=2, width = 0.5) +
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  facet_grid(~Genotype) +   
  scale_y_continuous(trans = "reverse") + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("0", "10", "15", "20", "30", "40", "50")) +
  scale_colour_manual(name="APA Training",
                      values=c("#f1a340", "#9970ab","#40004b"))+ 
  labs(x = "Stimulus Strenght (V)", y = "fEPSP Slope")

iobygeneonolegend <- ephys_long %>% 
  ggplot(aes(x=variablenumeric, y=value, color=APA )) + 
  stat_smooth(alpha=0.2, size=2) +
  geom_jitter(size=2, width = 0.5) +
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  facet_grid(~Genotype) +   
  scale_y_continuous(trans = "reverse") + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("0", "10", "15", "20", "30", "40", "50")) +
  scale_colour_manual(name="APA Training",
                      values=c("#f1a340", "#9970ab","#40004b"))+ 
  labs(x = "Stimulus Strenght (V)", y = "fEPSP Slope") + 
  theme(legend.position="none")


save_plot("ephysiobygeneotype.pdf", iobygeneo, base_height = 7, base_aspect_ratio = 1.5)
save_plot("ephysiobygeneotype.png", iobygeneo, base_height = 7, base_aspect_ratio = 1.5)


## plot just max response
ephysnoNA$Genotype <- factor(ephysnoNA$Genotype, 
                              levels = c("WT", "FMR1KO"))
ephysnoNA$APA <- factor(ephysnoNA$APA, 
                         levels = c("Yoked", "Trained", "Conflict"))


names(ephysnoNA)
slopmaxgeno <- ephysnoNA %>% 
  ggplot(aes(x=APA, y=min, fill=APA)) +
  geom_boxplot() +
  scale_y_continuous(trans = "reverse") + 
  background_grid(major = "xy", minor = "none") + 
  theme_cowplot(font_size = 20, line_size = 1) + 
  scale_fill_manual(name="APA Training",
                    values=c("#f1a340", "#9970ab","#40004b"))+ 
  labs(y = "Maximum fEPSP Slope", x = NULL) +
  facet_grid(~Genotype) 
  

save_plot("ephysslopmaxgeno.pdf", slopmaxgeno, base_aspect_ratio = 1.2, base_height = 8)
save_plot("ephysslopmaxgeno.png", slopmaxgeno, base_aspect_ratio = 1.2, base_height = 8)


ephsy_io_box_combo <-  plot_grid(iobygeneonolegend, slopmaxgeno, rel_widths=c(1,1.1))
save_plot("ephsy_io_box_combo.png", ephsy_io_box_combo, base_aspect_ratio = 2.5, base_height = 8)
save_plot("ephsy_io_box_combo.pdf", ephsy_io_box_combo, base_aspect_ratio = 2.5, base_height = 8)


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


