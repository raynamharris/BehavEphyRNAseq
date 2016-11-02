
Traits <- read.csv("JA16444samples.csv", sep=",", header = TRUE, stringsAsFactors=FALSE, na.string = "NA")
rownames(Traits) <- Traits$RNAseqID    # set $genoAPAsessionInd as rownames
names(Traits)

#keeping informative volumns
Traits <- Traits[c(1,3,5:6,10:11)] 
tail(Traits)
str(Traits)

## remove 100 and 100 animals because they aren't APA trained
## remove mice 147D_CA1_1 and 145B_CA3_1 because they produced almost 0 zeros
## remove 147 and 148 because they are home cage animals
Traits <- Traits %>% dplyr::filter(!grepl("100-|101-|147D-CA1-1|145B-CA3-1|147-|148-", RNAseqID)) 


## adding combinatorial traits
Traits$APAconflict <- as.factor(paste(Traits$APA, Traits$Conflict, sep="_")) 
summary(Traits)

## making it a numeric
Traits$Mouse <- as.integer(factor(Traits$Mouse))
Traits$Conflict <- as.integer(factor(Traits$Conflict))
Traits$APA <- as.integer(factor(Traits$APA))
Traits$Slice <- as.integer(factor(Traits$Slice))
Traits$APAconflict <- as.integer(factor(Traits$APAconflict))
head(Traits)
str(Traits)
summary(Traits)

## Subset by punch
TraitsCA1 <- Traits %>% filter(Punch == "CA1") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict)
TraitsCA3 <- Traits %>% filter(Punch == "CA3") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict)
TraitsDG <- Traits %>% filter(Punch == "DG") %>% select (RNAseqID, Mouse, Conflict, APA, Slice, APAconflict)


## make gene the row name then round all value to nearest 1s place
row.names(TraitsCA1) <- TraitsCA1$RNAseqID
TraitsCA1[1] <- NULL
row.names(TraitsCA3) <- TraitsCA3$RNAseqID
TraitsCA3[1] <- NULL
row.names(TraitsDG) <- TraitsDG$RNAseqID
TraitsDG[1] <- NULL


