## integrating behavior and ephys
#load("~/Github/BehavEphyRNAseq/bin/rnaseqpca.Rdata")

## behav ephys data
head(scoresdf)  # behavior pcas
head(ephys) # beahvior pcas with ephys
head(pcadata16)  # rnaseq pcas

## slim behavior ephy to top 10 pcs
scoresdfslim <- scoresdf[(c(1:10,79:83))]
ephysslim <- ephys[c(-10,-12,-14,-16, -17, -18)]

## rename the columsn
names(scoresdfslim)[names(scoresdfslim)=="Beahvior_PC1"] <- "Behavior_PC1"

names(scoresdfslim)[names(scoresdfslim)=="PC1"] <- "Behavior_PC1"
names(scoresdfslim)[names(scoresdfslim)=="PC2"] <- "Behavior_PC2"
names(scoresdfslim)[names(scoresdfslim)=="PC3"] <- "Behavior_PC3"
names(scoresdfslim)[names(scoresdfslim)=="PC4"] <- "Behavior_PC4"
names(scoresdfslim)[names(scoresdfslim)=="PC5"] <- "Behavior_PC5"
names(scoresdfslim)[names(scoresdfslim)=="PC6"] <- "Behavior_PC6"
names(scoresdfslim)[names(scoresdfslim)=="PC7"] <- "Behavior_PC7"
names(scoresdfslim)[names(scoresdfslim)=="PC8"] <- "Behavior_PC8"
names(scoresdfslim)[names(scoresdfslim)=="PC9"] <- "Behavior_PC9"
names(scoresdfslim)[names(scoresdfslim)=="PC10"] <- "Behavior_PC10"

names(ephysslim)[names(ephysslim)=="min"] <- "Max_Slope"
names(ephysslim)[names(ephysslim)=="V0"] <- "Ephys_10"
names(ephysslim)[names(ephysslim)=="V10"] <- "Ephys_V10"
names(ephysslim)[names(ephysslim)=="V20"] <- "Ephys_V20"
names(ephysslim)[names(ephysslim)=="V30"] <- "Ephys_V30"
names(ephysslim)[names(ephysslim)=="V40"] <- "Ephys_V40"
names(ephysslim)[names(ephysslim)=="V50"] <- "Ephys_V50"

head(scoresdfslim)  # behavior pcas
head(ephysslim) # beahvior pcas with ephys

## remove fmr1 ephys
ephysslim <- filter(ephysslim, Genotype == "WT") %>% droplevels()
ephysslim$genoAPA

## revalue genoAPA
scoresdfslim$APA <- revalue(scoresdfslim$APA, c("trained_trained" = "Trained")) 
scoresdfslim$APA <- revalue(scoresdfslim$APA, c("trained_conflict" = "Conflict")) 
scoresdfslim$APA <- revalue(scoresdfslim$APA, c("yoked_trained" = "Yoked")) 
scoresdfslim$APA <- revalue(scoresdfslim$APA, c("yoked_conflict" = "Yoked")) 
scoresdfslim$APA <- factor(scoresdfslim$APA, 
                         levels = c("Yoked", "Trained", "Conflict"))

tail(ephysslim)

## wrangle the rnaseq pca ----
rnaseqpcabyregion <- melt(pcadata16, by=7:10)
head(rnaseqpcabyregion)
rnaseqpcabyregion$measure <- as.factor(paste("RNAseq", rnaseqpcabyregion$variable, rnaseqpcabyregion$Punch, sep="_"))

# create a mouse id
rnaseqpcabyregion$ID <- as.factor(paste("15", rnaseqpcabyregion$name, sep=""))
rnaseqpcabyregion$ID <- sapply(strsplit(as.character(rnaseqpcabyregion$ID),'\\-'), "[", 1)
head(rnaseqpcabyregion)
rnaseqpcabyregion$ID <- as.factor(rnaseqpcabyregion$ID)


# not necessary
rnaseqpcabyregion$group <- NULL
rnaseqpcabyregion$Punch <- NULL
rnaseqpcabyregion$variable <- NULL
rnaseqpcabyregion$APA <- NULL
rnaseqpcabyregion$name <- NULL
head(rnaseqpcabyregion)

#subset by region
rnaseqpcaCA1 <- rnaseqpcabyregion %>%
  filter(grepl("CA1", measure))
rnaseqpcaCA3 <- rnaseqpcabyregion %>%
  filter(grepl("CA3", measure))
rnaseqpcaDG <- rnaseqpcabyregion %>%
  filter(grepl("DG", measure))


# spread  -- by region keeps more data
rnaseqpcabyregion <- dcast(rnaseqpcabyregion, ID ~ measure, value.var = "value")

## a full one but long
scoresdf_ephys_rnaseqpca <- inner_join(rnaseqpcabyregion, scoresdfslim, by="ID")
scoresdf_ephys_rnaseqpca <- left_join(scoresdf_ephys_rnaseqpca, ephysslim, by="ID")
names(scoresdf_ephys_rnaseqpca)
str(scoresdf_ephys_rnaseqpca)

## making a matrix for a heatmap
scoresdf_ephys_rnaseqpca_factors <- scoresdf_ephys_rnaseqpca[c(1,30:37)]
str(scoresdf_ephys_rnaseqpca_factors)
row.names(scoresdf_ephys_rnaseqpca_factors) <- scoresdf_ephys_rnaseqpca_factors$ID
names(scoresdf_ephys_rnaseqpca_factors)[names(scoresdf_ephys_rnaseqpca_factors)=="APA.x"] <- "APA"


scoresdf_ephys_rnaseqpca_matrix <- scoresdf_ephys_rnaseqpca[c(1,2:29,38:45)]
row.names(scoresdf_ephys_rnaseqpca_matrix) <- scoresdf_ephys_rnaseqpca_matrix$ID
scoresdf_ephys_rnaseqpca_matrix$ID <- NULL

scoresdf_ephys_rnaseqpca_matrix <- scale(scoresdf_ephys_rnaseqpca_matrix)
scoresdf_ephys_rnaseqpca_matrix <- t(scoresdf_ephys_rnaseqpca_matrix)


names(APA_annot)
## heatmap
APA_annot <- scoresdf_ephys_rnaseqpca_factors %>% select(APA)
APA_colors = list(
  APA =  c(Yoked = (values=c("#f1a340")), Trained = (values=c("#9970ab")), 
           Conflict = (values=c("#40004b"))))

quartz()
pheatmap(scoresdf_ephys_rnaseqpca_matrix, 
         show_colnames=FALSE, show_rownames=TRUE,
         annotation_col = APA_annot, 
         annotation_colors = APA_colors,
         fontsize = 20, fontsize_row = 10, 
         cellwidth=10, cellheight=10,
         width = 8,
         border_color = "grey60", filename = "pheatbehavephy.pdf"
)



## gg plots geom points

allplot <- function(data, xcol, ycol){
  print(xcol)
  plot <- data %>% 
    ggplot(aes_string(x = xcol, y = ycol)) +
    geom_point(size=5, aes(color=APA.x)) +
    stat_smooth(method=lm) +
    theme_cowplot(font_size = 20, line_size = 1) + 
    background_grid(major = "xy", minor = "none") + 
    #theme(strip.background = element_blank()) +
    scale_colour_manual(values=c( "#f1a340", "#9970ab","#40004b"),
                        name="APA Training")
  return(plot)
}
colnames(scoresdf_ephys_rnaseqpca)
allplot(data=scoresdf_ephys_rnaseqpca, xcol="Behavior_PC4", ycol="RNAseq_PC4_CA1")

#save_plot("rnaseqpca24.pdf", rnaseqpca24, base_aspect_ratio = 2)
#save_plot("rnaseqpca24.png", rnaseqpca24, base_aspect_ratio = 2)


## cormatrix
names(scoresdf_ephys_rnaseqpca)
rownames(scoresdf_ephys_rnaseqpca) <- scoresdf_ephys_rnaseqpca$ID
cormat <- scoresdf_ephys_rnaseqpca
cormat <- cormat[c(2:29,38:45)]
cormat <- cor(cormat, use = "complete.obs")

library(Hmisc)
rcorr(cormat, type = c("pearson"))
corematres <- cor(cormat, use = "complete.obs")
cormatres2 <- rcorr(as.matrix(cormat))
cormatres2

cormatres2$r  # get the R2 values
cormatres2$P        # pvalues

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
cormatres2flat <-  flattenCorrMatrix(cormatres2$r, cormatres2$P)
head(cormatres2flat)
min(cormatres2flat$cor)

library(corrplot)
quartz()
corrplot(cormat, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

library("PerformanceAnalytics")
quartz()
chart.Correlation(cormat, histogram=TRUE, pch=19)

col<- colorRampPalette(c("blue", "white", "red"))(20)
quartz()
heatmap(x = cormat, col = col, symm = TRUE)
view(cormat)


cormatdf <- as.data.frame(cormat)
rownames(cormatdf)
head(cormatdf)

cormatdf$rowcolors <- rownames(cormatdf)
cormatdf$rowcolors <- sapply(strsplit(as.character(cormatdf$rowcolors),'\\_'), "[", 1)
cormatdf$rowcolors <- ifelse(grepl("RNAseq", cormatdf$rowcolors), "green", 
                                    ifelse(grepl("Behavior", cormatdf$rowcolors), "blue", "red"))
            
cormatrowsidecols <- as.vector(cormatdf$rowcolors)
quartz()
beahv_ephys_rna_cormat <- heatmap(x = cormat, col = col, symm = TRUE, 
        RowSideColors = cormatrowsidecols, na.rm = T,
        labCol = F)

save_plot("beahv_ephys_rna_cormat.pdf", beahv_ephys_rna_cormat)

