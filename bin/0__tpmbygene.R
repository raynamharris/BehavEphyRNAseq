# merge tpm and gene id dataframe
tpmbygene <-  full_join(geneids, tpm)
str(tpmbygene)

## remove unnecesary columns (aka, keep gene name and counts for samples)
tpmbygene <- tpmbygene[-c(1:6,8:12)]  

## lenghten 
tpmbygene <- melt(tpmbygene, id=c("gene"))
head(tpmbygene)

## remove 100 and 100 animals because they aren't APA trained
## remove mice 147D_CA1_1 and 145B_CA3_1 because they produced almost 0 zeros
## remove 147 and 148 because they are home cage animals
tpmbygene <- tpmbygene %>% dplyr::filter(!grepl("100-|101-|147D-CA1-1|145B-CA3-1|147-|148-", variable)) 

## create three separate files for CA1 CA3 and DG
head(tpmbygene)
tpmbygeneCA1 <- filter(tpmbygene, grepl("CA1", variable))
tpmbygeneCA3 <- filter(tpmbygene, grepl("CA3", variable))
tpmbygeneDG <- filter(tpmbygene, grepl("DG", variable))

#then widen by sum
tpmbygeneCA1 <- dcast(tpmbygeneCA1, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygeneCA3 <- dcast(tpmbygeneCA3, gene ~ variable, value.var= "value", fun.aggregate=mean)
tpmbygeneDG <- dcast(tpmbygeneDG, gene ~ variable, value.var= "value", fun.aggregate=mean)


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
