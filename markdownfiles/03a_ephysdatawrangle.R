setwd("~/Github/BehavEphyRNAseq/data/ephys")

## read and transform data

ephys <- read.csv("Ephy_forDataAnalysisV2.csv", header=T, na.strings=c(""), stringsAsFactors = FALSE)
str(ephys)

## rename columns 
names(ephys)[names(ephys)=="Peak..V."] <- "Peak V"
names(ephys)[names(ephys)=="training"] <- "group"
str(ephys)
head(ephys)

## new cols
ephys$APA <- ifelse(grepl("TTT", ephys$group), "Trained", 
                     ifelse(grepl("TCT", ephys$group), "Conflict",
                            ifelse(grepl("TTY", ephys$group), "Yoked",
                                   ifelse(grepl("TCY", ephys$group), "Yoked", NA))))
                                                                             
ephys$Year <- ifelse(grepl("151", ephys$ID), "2015", 
                    ifelse(grepl("16-", ephys$ID), "2016", NA))

ephys$genoAPA <- as.factor(paste(ephys$Genotype,ephys$APA,sep="_"))


## reorder
head(ephys)
names(ephys)
ephys <- ephys[c(2,15:18,3:14)]
str(ephys)

## make characters either factors or numbers, as appropriate
colsfactor = c(1:4)
colsnumeric = c(17)
ephys[,colsfactor] %<>% lapply(function(x) as.factor(as.factor(x)))
ephys[,colsnumeric] %<>% lapply(function(x) as.numeric(as.integer(x)))

str(ephys)
summary(ephys)

## make with max column with max i/o
ephys$max <- apply(ephys[c(6:16)], 1, max,na.rm=TRUE)
ephys$min <- apply(ephys[c(6:16)], 1, min,na.rm=TRUE)
summary(ephys)

## make new iomax  df
iomax <- ephys[c(6:18,20)]
summary(iomax)
iomax <- iomax[,-c(1,13)]/iomax[,14]
iomax$ID <- ephys$ID
iomax$APA <- ephys$APA
iomax$Year <- ephys$Year
iomax$Genotype <- ephys$Genotype
iomax$genoAPA <- ephys$genoAPA
iomax$min <- NULL

## checking nas
summary(ephys)
ephysnoNA <- ephys[c(-10,-12, -14, -16, -17)]
head(ephysnoNA)


