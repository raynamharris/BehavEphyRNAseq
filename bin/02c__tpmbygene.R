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

## remove 147 and 148 because they are home cage animals
tpmbygene <- tpmbygene %>% dplyr::filter(!grepl("147-|148-", variable)) 
countbygene <- countbygene %>% dplyr::filter(!grepl("147-|148-", variable)) 

#then widen by sum
tpmbygene <- dcast(tpmbygene, gene ~ variable, value.var= "value", fun.aggregate=mean)
countbygene  <- dcast(countbygene, gene ~ variable, value.var= "value", fun.aggregate=mean)

## make gene the row name then round all value to nearest 1s place
row.names(tpmbygene) <- tpmbygene$gene
tpmbygene[1] <- NULL
tpmbygene <- round(tpmbygene)
summary(tpmbygene)

row.names(countbygene) <- countbygene$gene
countbygene[1] <- NULL
countbygene <- round(countbygene)
summary(countbygene)
