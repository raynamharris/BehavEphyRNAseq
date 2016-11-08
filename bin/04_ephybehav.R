## integrating behavior and ephys
## 

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
ggplot(scoresdf_ephys, aes(PC1, min, color=APA.x)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point(size = 8, aes(shape=Genotype)) + 
  scale_y_continuous(trans = "reverse") +
  #facet_wrap(~Genotype) +
  theme_cowplot(font_size = 20, line_size = 1) + 
  background_grid(major = "xy", minor = "none") + 
  #theme(strip.background = element_blank()) +
  scale_colour_manual(name="APA Training",
                      breaks=c("trained_conflict", "trained_trained", "yoked_conflict", "yoked_trained"),
                      labels=c("Trained Conflict", "Trained Trained", "Yoked to Conflict", "Yoked to Trained"),
                      values=c("#7f3b08", "#f1a340", "#40004b","#9970ab"))  

## function
beahvepcaphys <- function(data, xcol, ycol, colorcode){
  print(ycol)
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size = 8, aes(shape=Genotype)) + 
    geom_smooth(method = "lm", se = FALSE) + 
    scale_y_continuous(trans = "reverse") +
    ylab("Input Output Slope Maximum") +
    facet_wrap(~Genotype) +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    #theme(strip.background = element_blank()) +
    scale_colour_manual(name="APA Training",
                        breaks=c("trained_conflict", "trained_trained", "yoked_conflict", "yoked_trained"),
                        labels=c("Trained Conflict", "Trained Trained", "Yoked to Conflict", "Yoked to Trained"),
                        values=c("#7f3b08", "#f1a340", "#40004b","#9970ab"))  
  return(plot)
}

beahvepcaphys(data=scoresdf_ephys, xcol="PC1", ycol="min", colorcode="APA.x")
beahvepcaphys(data=scoresdf_ephys_wt2015, xcol="PC1", ycol="min", colorcode="APA.x")
beahvepcaphys(data=scoresdf_ephys_wt2015, xcol="PC2", ycol="min", colorcode="APA.x")
beahvepcaphys(data=scoresdf_ephys_wt2015, xcol="PC3", ycol="min", colorcode="APA.x")
beahvepcaphys(data=scoresdf_ephys_wt2015, xcol="PC4", ycol="min", colorcode="APA.x")

#write.csv(loadings, "loadings.csv", row.names = T)



lm1 <- lm(min ~ APA.x + PC1 + PC2 + PC3 + PC4 + PC5, data = scoresdf_ephys)
summary(lm1)

## subset to just numeric columsn for corrlation
names(scoresdf_ephys)
#scoresdf_ephys_numeric <- scoresdf_ephys[c(1:10,41:44,46,48,50,54)]
scoresdf_ephys_numeric <- scoresdf_ephys[c(1:5,50,54,55)]
scoresdf_ephys_numeric <- scale(scoresdf_ephys_numeric)

rownames(scoresdf_ephys_numeric) <- scoresdf_ephys$ID

cormat <- round(cor(scoresdf_ephys_numeric),2)
melted_cormat <- melt(cormat)
head(melted_cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
