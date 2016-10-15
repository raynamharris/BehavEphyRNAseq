#Sleuth from https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html

#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#biocLite("biomaRt")
#install.packages('devtools')
#devtools::install_github('pachterlab/sleuth')
library("sleuth")
library("dplyr")

## set paths to sample info and kallisto data
setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444")
base_dir <- "~/Github/BehavEphyRNAseq/TACC-copy/JA16444"
sample_id <- dir(file.path(base_dir,"03_kallistoquantcandidategenes/"))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "03_kallistoquantcandidategenes/", id))
kal_dirs

## read in the sample info, rename and revalue things
s2c <- read.csv("JA16444samples.csv", sep=",", header = TRUE, stringsAsFactors=FALSE, na.string = "NA")
s2c <- dplyr::rename(s2c, sample = RNAseqID) ## unique identifier must be called sample
s2c$APA[is.na(s2c$APA)] <- "noAPA" ## make NA more meaningful
s2c$Behavior[is.na(s2c$Behavior)] <- "noAPA" ## make NA more meaningful
s2c$E.phy[is.na(s2c$E.phy)] <- "noAPA" ## make NA more meaningful
s2c$Conflict[is.na(s2c$Conflict)] <- "noAPA" ## make NA more meaningful

### for full transcriptome analysis only
## Not all samples had mapped reads....created a files.csv with only mapped samples, used this to filter
#files <- read.csv("files.csv",  sep=",", header = TRUE, stringsAsFactors=FALSE)
#s2c <- inner_join(s2c, files)
#s2c

## add more grouping columns
s2c$APAconflict <- as.factor(paste(s2c$APA, s2c$Conflict, sep="_")) 
s2c$APApunch <- as.factor(paste(s2c$APA, s2c$Punch, sep="_")) 
s2c$ConflictPunch <- as.factor(paste(s2c$Conflict, s2c$Punch, sep="_")) 
s2c$APAconflictPunch <- as.factor(paste(s2c$APA, s2c$Conflict, s2c$Punch, sep="_")) 
s2c
## only APA trained or yoked animals

## add the path to the file with mutate
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
names(s2c)
## remove bad files
s2c <- s2c %>% filter(!grepl("147D-CA1-1", path)) 
s2c <- s2c %>% filter(!grepl("145B-CA3-1", path)) 


## subset the data 
noAPA <- filter(s2c, APA == "noAPA")
withAPA <- filter(s2c, APA != "noAPA")
CA1 <- filter(s2c, Punch == "CA1")
DG <- filter(s2c, Punch == "DG")
CA3 <- filter(s2c, Punch == "CA3")
conflict <- filter(s2c, Conflict == "Conflict")


## looking at just the tissue dilution samples

## Now the “sleuth object” can be constructed. 
## This requires three commands that 
## (1) load the kallisto processed data into the object 
## (2) estimate parameters for the sleuth response error measurement model and 
## (3) perform differential analysis (testing). On a laptop the three steps should take about 2 minutes altogether.
so <- sleuth_prep(s2c, ~ APA)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'APAYoked') 
models(so)
sleuth_live(so)
results_table <- sleuth_results(so, 'APAYoked')
