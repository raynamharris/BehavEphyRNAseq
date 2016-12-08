library(dplyr)
library(reshape2)

# merge tpm and gene id dataframe
tpmbygene <-  full_join(geneids, tpm)
str(tpmbygene)

countbygene <- full_join(geneids, count)
str(countbygene)

## remove unnecesary columns (aka, keep gene name and counts for samples)
tpmbygene <- tpmbygene[-c(1:6,8:12)]  
countbygene <- countbygene[-c(1:6,8:12)]  


## lenghten 
tpmbygene <- melt(tpmbygene, id=c("gene"))
head(tpmbygene)

countbygene <- melt(countbygene, id=c("gene"))
head(countbygene)

## remove 100 and 100 animals because they aren't APA trained
## remove mice 147D_CA1_1 and 145B_CA3_1 because they produced almost 0 zeros
## remove 147 and 148 because they are home cage animals
tpmbygene <- tpmbygene %>% dplyr::filter(!grepl("100-|101-|147D-CA1-1|145B-CA3-1|147-|148-", variable)) 
countbygene <- countbygene %>% dplyr::filter(!grepl("100-|101-|147D-CA1-1|145B-CA3-1|147-|148-", variable)) 


## create three separate files for CA1 CA3 and DG
head(tpmbygene)
tail(tpmbygene)
tpmbygeneCA1 <- filter(tpmbygene, grepl("CA1", variable))
tpmbygeneCA3 <- filter(tpmbygene, grepl("CA3", variable))
tpmbygeneDG <- filter(tpmbygene, grepl("DG", variable))
tpmbygeneCA1DG <- filter(tpmbygene, grepl("DG|CA1", variable))
tpmbygenetrained <- filter(tpmbygene, grepl("B|D", variable))
tpmbygeneyoked <- filter(tpmbygene, grepl("A|C", variable))

countbygeneCA1DG <- filter(countbygene, grepl("DG|CA1", variable))


#then widen by sum
tpmbygeneCA1 <- dcast(tpmbygeneCA1, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygeneCA3 <- dcast(tpmbygeneCA3, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygeneDG <- dcast(tpmbygeneDG, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygenetrained <- dcast(tpmbygenetrained, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygeneyoked <- dcast(tpmbygeneyoked, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygeneCA1DG <- dcast(tpmbygeneCA1DG, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygene <- dcast(tpmbygene, gene ~ variable, value.var= "value", fun.aggregate=mean)

countbygeneCA1DG <- dcast(countbygeneCA1DG, gene ~ variable, value.var= "value", fun.aggregate=mean)
countbygene  <- dcast(countbygene, gene ~ variable, value.var= "value", fun.aggregate=mean)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygeneCA1) <- tpmbygeneCA1$gene
tpmbygeneCA1[1] <- NULL
tpmbygeneCA1 <- round(tpmbygeneCA1)
summary(tpmbygeneCA1)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygeneCA3) <- tpmbygeneCA3$gene
tpmbygeneCA3[1] <- NULL
tpmbygeneCA3 <- round(tpmbygeneCA3)
summary(tpmbygeneCA3)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygeneDG) <- tpmbygeneDG$gene
tpmbygeneDG[1] <- NULL
tpmbygeneDG <- round(tpmbygeneDG)
summary(tpmbygeneDG)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygeneCA1DG) <- tpmbygeneCA1DG$gene
tpmbygeneCA1DG[1] <- NULL
tpmbygeneCA1DG <- round(tpmbygeneCA1DG)
summary(tpmbygeneCA1DG)


## make gene the row name then round all value to nearest 1s place
row.names(tpmbygene) <- tpmbygene$gene
tpmbygene[1] <- NULL
tpmbygene <- round(tpmbygene)
summary(tpmbygene)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygenetrained) <- tpmbygenetrained$gene
tpmbygenetrained[1] <- NULL
tpmbygenetrained <- round(tpmbygenetrained)
summary(tpmbygenetrained)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygeneyoked) <- tpmbygeneyoked$gene
tpmbygeneyoked[1] <- NULL
tpmbygeneyoked <- round(tpmbygeneyoked)
summary(tpmbygeneyoked)

## make gene the row name then round all value to nearest 1s place
row.names(countbygeneCA1DG) <- countbygeneCA1DG$gene
countbygeneCA1DG[1] <- NULL
countbygeneCA1DG <- round(countbygeneCA1DG)
summary(countbygeneCA1DG)

## make gene the row name then round all value to nearest 1s place
row.names(countbygene) <- countbygene$gene
countbygene[1] <- NULL
countbygene <- round(countbygene)
summary(countbygene)
