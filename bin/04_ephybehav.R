library(reshape2)
scoresdf$ID
ephys$ID

## was having problems with joining, so had to trim white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
ephys$ID <- trim(ephys$ID)
scoresdf$ID <- trim(scoresdf$ID)

scoresdf$ID <- as.factor(as.character(scoresdf$ID))
ephys$ID <- as.factor(as.character(ephys$ID))

scoresdf$ID
ephys$ID

scoresdf_ephys <- inner_join(scoresdf, ephys, by="ID")
str(scoresdf_ephys)
head(scoresdf_ephys$APA.x)

scoresdf_ephys$training <- ifelse(grepl("trained_conflict", scoresdf_ephys$APA.x), "trained", 
                    ifelse(grepl("trained_trained", scoresdf_ephys$APA.x), "trained",
                           ifelse(grepl("yoked_conflict", scoresdf_ephys$APA.x), "yoked",
                                  ifelse(grepl("yoked_trained", scoresdf_ephys$APA.x), "yoked", NA))))

names(scoresdf_ephys)[names(scoresdf_ephys)=="Genotype.y"] <- "Genotype"
levels(scoresdf_ephys$Genotype)
scoresdf_ephys$Genotype <- factor(scoresdf_ephys$Genotype, levels = c("WT", "FMR1KO"))

levels(scoresdf_ephys$genoYear)
scoresdf_ephys_wt2015 <- scoresdf_ephys %>%
  filter(genoYear == "WT_2015")


## behavpca ephs plot
behavephyscorrlationPC2 <- ggplot(scoresdf_ephys, aes(PC2, min, color=APA.x)) + 
  geom_smooth(method = "lm", se = TRUE, alpha=0.2) + 
  geom_point(size = 8, aes(shape=Genotype)) + 
  scale_y_continuous(trans = "reverse") +
  facet_wrap(~APA.x) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "none", minor = "none") + 
  #theme(strip.background = element_blank()) +
  scale_colour_manual(values=c( "#f1a340", "#9970ab","#40004b"),
                      name="APA Training") + 
  labs(x = "Behavior PC2", y = "Maximum fEPSE Slope")

save_plot("behavephyscorrlationPC2.pdf", behavephyscorrlationPC2, base_aspect_ratio = 1.3, base_height = 8)


library(plyr)
ddply(scoresdf_ephys, .(APA.x), summarise, "corr" = cor.test(PC1, min, method = "spearman"))
ddply(scoresdf_ephys, .(APA.x), summarise, "corr" = cor.test(PC2, min, method = "spearman"))
ddply(scoresdf_ephys, .(APA.x), summarise, "corr" = cor.test(PC5, min, method = "spearman"))


## function
beahvepcaphys <- function(data, xcol, ycol, colorcode){
  print(ycol)
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size = 8, aes(shape=Genotype)) + 
    geom_smooth(method = "lm", se = FALSE) + 
    scale_y_continuous(trans = "reverse") +
    ylab("Input Output Slope Maximum") +
    #facet_wrap(~Genotype) +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    #theme(strip.background = element_blank()) +
    scale_colour_manual(values=c( "#f1a340", "#7f3b08", "#9970ab","#40004b"),
                        name="APA Training",
                        breaks=c("yoked_trained", "yoked_conflict", "trained_trained", "trained_conflict"),
                        labels=c("Yoked to Trained", "Yoked to Conflict", "Trained Trained", "Trained Conflict")) 
  return(plot)
}

beahvepcaphys(data=scoresdf_ephys, xcol="PC1", ycol="min", colorcode="APA.x")

#write.csv(loadings, "loadings.csv", row.names = T)


