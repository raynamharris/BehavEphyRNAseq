setwd("~/Github/BehavEphyRNAseq/data/ephys")
library(magrittr)  # for function function "%<>%"
#source("figureoptions.R")

## read and transform data

ephys2 <- read.csv("Ephy_forDataAnalysisV2.csv", header=T, na.strings=c(""), stringsAsFactors = FALSE)

## rename columns 
names(ephys2)[names(ephys2)=="Peak..V."] <- "Peak V"
names(ephys2)[names(ephys2)=="training"] <- "group"

## new cols
ephys2$APA <- ifelse(grepl("TTT", ephys2$group), "consistent", 
                     ifelse(grepl("TCT", ephys2$group), "conflict",
                            ifelse(grepl("TTY", ephys2$group), "control",
                                   ifelse(grepl("TCY", ephys2$group), "control", NA))))
     
ephys2$APA2 <- ifelse(grepl("TTT", ephys2$group), "consistent", 
                     ifelse(grepl("TCT", ephys2$group), "conflict",
                            ifelse(grepl("TTY", ephys2$group), "control-consistent",
                                   ifelse(grepl("TCY", ephys2$group), "control-conflict", NA))))

ephys2$Year <- ifelse(grepl("151", ephys2$ID), "2015", 
                    ifelse(grepl("16-", ephys2$ID), "2016", NA))

ephys2$genoAPA <- as.factor(paste(ephys2$Genotype,ephys2$APA,sep="_"))



## reorder, drop old group coloumn
ephys2 <- ephys2[c(2,15:19,3:14)]

### Drop the voltages with very little recording and the Vmax
ephys2 <- ephys2[c(-11,-13, -15, -17, -18)]


## make characters either factors or numbers, as appropriate
str(ephys2)
colsfactor = c(1:6)
#colsnumeric = c(18)
ephys2[,colsfactor] %<>% lapply(function(x) as.factor(as.factor(x)))
#ephys2[,colsnumeric] %<>% lapply(function(x) as.numeric(as.integer(x)))

## make with max column with max i/o
ephys2$min <- apply(ephys2[c(7:13)], 1, min,na.rm=TRUE)
ephys2$max <- apply(ephys2[c(7:13)], 1, max,na.rm=TRUE)


## prep levels for visualization
ephys2$Genotype <- factor(ephys2$Genotype, 
                              levels = c("WT", "FMR1KO"))
ephys2$APA <- factor(ephys2$APA, 
                         levels = c("control", "consistent", "conflict"))




## make long and tidy
ephys2_long <- melt(ephys2, id = c(1:6))

## add numeric for stat smooth
ephys2_long$variablenumeric <- ifelse((ephys2_long$variable == "V0"), "1", 
                                      ifelse(grepl("V10", ephys2_long$variable ), "2",
                                             ifelse(grepl("V15", ephys2_long$variable ), "3",
                                                    ifelse(grepl("V20", ephys2_long$variable), "4", 
                                                           ifelse(grepl("V30", ephys2_long$variable), "5",
                                                                  ifelse(grepl("V40", ephys2_long$variable), "6",
                                                                         ifelse(grepl("V50", ephys2_long$variable), "7", NA)))))))
ephys2_long <- ephys2_long %>% drop_na()

levels(ephys2_long$variable)
ephys2_long$variablenumeric <- as.numeric(ephys2_long$variablenumeric)

ephys2_long$Genotype <- factor(ephys2_long$Genotype, 
                               levels = c("WT", "FMR1KO"))
ephys2_long$APA <- factor(ephys2_long$APA, 
                          levels = c("control", "consistent", "conflict"))
######################################
## make a WT 2016 only
ephys2015 <- ephys2 %>%
  filter(Year %in% c("2015"))
ephys2015_long <- ephys2_long %>%
  filter(Year %in% c("2015"))


write.csv(ephys2015, "~/Github/IntegrativeProjectWT2015/data/03_ephys2.csv", row.names = F)
write.csv(ephys2015_long, "~/Github/IntegrativeProjectWT2015/data/03_ephys3.csv", row.names = F)

