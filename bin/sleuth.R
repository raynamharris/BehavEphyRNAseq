#Sleuth from https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html

source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
biocLite("biomaRt")
install.packages('devtools')
devtools::install_github('pachterlab/sleuth')
library("sleuth")

## set paths to sample info and kallisto data
setwd("~/Github/BehavEphyRNAseq/TACC-copy/JA16444")
base_dir <- "~/Github/BehavEphyRNAseq/TACC-copy/JA16444"
sample_id <- dir(file.path(base_dir,"02_kallistoquant"))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "02_kallistoquant", id))
kal_dirs

## read in the sample info
s2c <- read.csv("JA16444samples.csv", sep=",", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::rename(s2c, sample = RNAseqID)
s2c


## remove non APA animals
s2c <- filter(s2c, APA != "NA")

## Not all samples had mapped reads....created a files.csv with only mapped samples, used this to filter
files <- read.csv("files.csv",  sep=",", header = TRUE, stringsAsFactors=FALSE)
s2c <- inner_join(s2c, files)
s2c

## add the path to the file with mutate
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

## subset the data to remove sample without


## looking at just the tissue dilution samples

## Now the “sleuth object” can be constructed. 
## This requires three commands that 
## (1) load the kallisto processed data into the object 
## (2) estimate parameters for the sleuth response error measurement model and 
## (3) perform differential analysis (testing). On a laptop the three steps should take about 2 minutes altogether.
so <- sleuth_prep(s2c, ~  APA + Conflict + Punch)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'APAYoked') 
models(so)
sleuth_live(so)
results_table <- sleuth_results(so, 'APAYoked')

sp <- sleuth_wt(so, 'PunchCA3') 
models(sp)
sleuth_live(sp)


## plotting densities by group
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$APAYoked),
                                                     "sample"), offset = 1)


setwd("~/Github/SingleNeuronSeq/results")
library(dplyr)
ERCCs <- read.csv("processed_data_table_ERCC_2016-03-25.csv", sep = ",", header = TRUE)
mus <- read.csv("processed_data_table_Mus_tissue_neuron.csv", sep = ",", header = TRUE)
combo <- dplyr::full_join(mus, ERCCs, by = "sample")
combo <- combo[c(-1, -6, -7,-8,-9,-10,-11,-12,-13,-22)]
combo <- dplyr::mutate(combo, total_mapped = frac_mapped.x + frac_mapped.y)

library(plyr)
combo <- rename(combo, c("frac_mapped.x" = "frac_mapped_mus"))
combo <- rename(combo, c("frac_mapped.y" = "frac_mapped_ERCC"))
plot(combo$frac_mapped.mus, combo$frac_mapped.ERCC)

library(ggplot2)
ggplot(combo, aes(frac_mapped_ERCC, frac_mapped_mus), size = 3) + geom_point(size= 3) + geom_point(aes(colour = factor(RNAseqBatch.y, size = 3)))
theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
#write.csv(combo, "combo.csv")


## now plot some ERCC stuff

library(dplyr)
library(ggplot2)
library(reshape2)
setwd("~/Github/SingleNeuronSeq/results")
ERCC <- read.csv("kallisto_table_ERCC.csv", sep=",", header=TRUE)
#plot tpm versus counts
a <- ggplot(ERCC, aes(x = target_id, y = est_counts))
a + geom_point(aes(color=factor(batch))) + theme(legend.position="none")

## make the data wide so I can look at correlations across samples
ERCC_counts <- select(ERCC,target_id, sample, est_counts)
ERCC_counts <- dcast(ERCC_counts, target_id ~ sample)
ERCC_counts <- data.frame(ERCC_counts[,-1], row.names=ERCC_counts[,1])
ERCC_tpm <- select(ERCC,target_id, sample, tpm)
ERCC_tpm <- dcast(ERCC_tpm, target_id ~ sample)
ERCC_tpm <- data.frame(ERCC_tpm[,-1], row.names=ERCC_tpm[,1])

## make correlation matrix
library(Hmisc)
ERCC_counts <- log(ERCC_counts + 1) 
ERCC_counts_corr <- cor(ERCC_counts, use="complete.obs", method="kendall") 
install.packages("corrgram")
library(corrgram)
corrgram(ERCC_counts, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="Car Milage Data (unsorted)")


require(lattice)
levelplot(ERCC_counts_corr)


install.packages("qtlcharts")
library(qtlcharts)
iplotCorr(mat=ERCC_counts, group=ERCC_counts$batch, reorder=TRUE)

ERCC_counts_corr <- round(cor(ERCC_counts),2)
ERCC_counts_corr_melt <- melt(ERCC_counts_corr)
ggplot(data = ERCC_counts_corr_melt, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = ERCC_counts_corr, col = col, symm = TRUE)
