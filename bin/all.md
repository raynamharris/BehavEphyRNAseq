The Project
-----------

The research project was designed to understand how experience shapes
the brain. In particular, we are looking at learned avoidance behavhior.
How do animals change their behavior to avoidance an unpleasant
experience?

RNAseq samples
--------------

In the summers of 2015 and 2016, I processed a bunch of hippocampal
tissue samples from WT and FMR1-KO mice. Most mice were trained in an
active place avoidance task or used as yoked controls; however, a few
animals were taken straight from the home cage.

This data has been cleaned using a different script `02a_punches.R`.

THis

    ##       RNAseqID      Mouse      year    Genotype    jobnumber  Punch   
    ##  100-CA1-1: 1   15-100 :14   2015:71   FMR1: 9   JA16268: 4   CA1:43  
    ##  100-CA1-2: 1   15-143C: 3   2016:17   WT  :79   JA16444:67   CA3:21  
    ##  100-CA1-3: 1   15-144A: 3                       JA17009:17   DG :24  
    ##  100-CA3-1: 1   15-144C: 3                                            
    ##  100-CA3-4: 1   15-145A: 3                                            
    ##  100-DG-2 : 1   15-145B: 3                                            
    ##  (Other)  :82   (Other):59                                            
    ##         Group          Conflict       APA             method  
    ##  conflict  :14   Conflict  :36   Trained:28   dissociated: 7  
    ##  consistent:14   NoConflict:32   Yoked  :40   homogenized:81  
    ##  control   :40   NA's      :20   NA's   :20                   
    ##  homecage  :20                                                
    ##                                                               
    ##                                                               
    ##                                                               
    ##      dodgy          daytime   Slice       Date   
    ##  allgood:77   afternoon : 3   1:22   9/28/15:14  
    ##  ephys  : 2   beforenoon:15   2:31   7/18/15: 9  
    ##  slice  : 9   earlyAM   : 6   3:22   7/23/15: 6  
    ##               evening   : 5   4:13   7/24/15: 6  
    ##               nighttime : 3          7/29/15: 6  
    ##               norecord  :56          7/30/15: 6  
    ##                                      (Other):41

Kallisto Gather
---------------

The kallisto output gives you read counts for sample in an abundance
file for every single sample. This portion of the code goes through and
finds each samples' abundance.tsv file, extracts the data, and combines
it all into a dataframe. The "counts" file is unnormalized, but the
"tpm" is the data after being normalized by transcripts per million. I
also use some string splitting to take the very long transcript
identifying and create a "geneids" file that has all the database
identifiers for each transcript.

(P.S. Unfortunately, I have no idea how to do this next part without
changing directories.)

\`\`\`{r kallisto gather}
=========================

setwd("../data/JA16444/") \#\# this will create lists of all the samples
kallistoDirs = dir(".") kallistoDirs =
kallistoDirs\[!grepl("\\.(R|py|pl|sh|xlsx?|txt|tsv|csv|org|md|obo|png|jpg|pdf)$",
kallistoDirs, ignore.case=TRUE)\]

kallistoFiles = paste0(kallistoDirs, "/abundance.tsv")
names(kallistoFiles) = kallistoDirs if(file.exists(kallistoFiles))
kallistoData = lapply( kallistoFiles, read.table, sep = "", row.names =
1, header = TRUE )

this for loop uses the reduce function to make two data frame with counts or tpm from all the samples
-----------------------------------------------------------------------------------------------------

ids = Reduce(f=union, x=lapply(kallistoData, rownames)) if
(all(sapply(kallistoData, function(x) {all(rownames(x)==ids)}))) { count
= data.frame( id = ids, sapply(kallistoData, function(x)
{x$est\_counts}), check.names = FALSE, stringsAsFactors = FALSE ) tpm = data.frame( id = ids, sapply(kallistoData, function(x) {x$tpm}),
check.names = FALSE, stringsAsFactors = FALSE ) }

make a dataframe with the parts of the gene id as columns
---------------------------------------------------------

geneids &lt;- count\[c(1)\]
geneids*E**N**S**M**U**S**T* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 1)
geneids*E**N**S**M**U**S**G* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 2)
geneids*O**T**T**M**U**S**G* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 3)
geneids*O**T**T**M**U**S**T* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 4)
geneids*t**r**a**n**s**c**r**i**p**t* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 5)
geneids*g**e**n**e* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 6)
geneids*l**e**n**g**t**h* &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 7)
geneids*s**t**r**u**c**t**u**r**e*1 &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 8)
geneids*s**t**r**u**c**t**u**r**e*2 &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 9)
geneids*s**t**r**u**c**t**u**r**e*3 &lt; −*s**a**p**p**l**y*(*s**t**r**s**p**l**i**t*(*a**s*.*c**h**a**r**a**c**t**e**r*(*g**e**n**e**i**d**s*id),'\\|'),
"\[", 10)
geneids*t**r**a**n**s**c**r**i**p**t*<sub>*l*</sub>*e**n**g**h**t* &lt; −*a**s*.*f**a**c**t**o**r*(*p**a**s**t**e*(*g**e**n**e**i**d**s*transcript,
geneids$length, sep="\_"))

prep data for wgcna
-------------------

countswgcna &lt;- count row.names(countswgcna) &lt;-
geneids$transcript\_lenght countswgcna\[1\] &lt;- NULL countswgcna &lt;-
round(countswgcna) summary(countswgcna)

prep data for wgcna
-------------------

tpmswgcna &lt;- tpm row.names(tpmswgcna) &lt;-
geneids$transcript\_lenght tpmswgcna\[1\] &lt;- NULL tpmswgcna &lt;-
round(tpmswgcna) summary(tpmswgcna) \#\`\`\`

Merge transcipts counts to gene counts
--------------------------------------

Kallisto is cool because it does give you transcript level stuff, but
right now I think I have more power if I just look at gene level stuff.
I'll go back and look at transcripts if I want to.

\`\`\`{r tpmbygene}
===================

merge tpm and gene id dataframe
===============================

tpmbygene &lt;- full\_join(geneids, tpm) str(tpmbygene) head(tpmbygene)

countbygene &lt;- full\_join(geneids, count) str(countbygene)

remove unnecesary columns (aka, keep gene name and counts for samples)
----------------------------------------------------------------------

tpmbygene &lt;- tpmbygene\[-c(1:6,8:12)\]  
countbygene &lt;- countbygene\[-c(1:6,8:12)\]

lenghten
--------

tpmbygene &lt;- melt(tpmbygene, id=c("gene")) head(tpmbygene)

countbygene &lt;- melt(countbygene, id=c("gene")) head(countbygene)

then widen by sum
=================

tpmbygene &lt;- dcast(tpmbygene, gene ~ variable, value.var= "value",
fun.aggregate=mean) countbygene &lt;- dcast(countbygene, gene ~
variable, value.var= "value", fun.aggregate=mean)

make gene the row name then round all value to nearest 1s place
---------------------------------------------------------------

row.names(tpmbygene) &lt;- tpmbygene$gene tpmbygene\[1\] &lt;- NULL
tpmbygene &lt;- round(tpmbygene) summary(tpmbygene) head(tpmbygene)

row.names(countbygene) &lt;- countbygene$gene countbygene\[1\] &lt;-
NULL countbygene &lt;- round(countbygene) summary(countbygene)

\`\`\`
======

DESeq Analysis
--------------

Now, I'll look for differential gene expression between the FMR1-KO and
WT mice. This analysis was developed by reading the DESEq manual. In
many place, I try to provide the chapter where these steps are described
in more details.

\`\`\`{r DESeq}
===============

1.3.3 Count matrix input ----
=============================

countData &lt;- countbygene colData &lt;- Traits %&gt;%
arrange(RNAseqID) \# needs to be in same order a countData
head(countData) head(colData)

making sure colData and countData have the same number of rows
--------------------------------------------------------------

savecols &lt;- as.character(colData$RNAseqID) \#select the sample name
column that corresponds to row names savecols &lt;- as.vector(savecols)
\# make it a vector countData &lt;- countData %&gt;%
select(one\_of(savecols)) \# select just the columns that match the
samples in colData

remove genes with total counts across all samples &lt; 2
--------------------------------------------------------

countData\[countData &lt; 2\] &lt;- 0

differential gene expression
----------------------------

dds &lt;- DESeqDataSetFromMatrix(countData = countData, colData =
colData, design = ~ Group + Punch) dds

1.3.6 Pre-filtering
-------------------

dds &lt;- dds\[ rowSums(counts(dds)) &gt; 1, \]

1.3.7 Note on factor levels
---------------------------

dds*G**r**o**u**p* &lt; −*f**a**c**t**o**r*(*d**d**s*Group,
levels=c("homecage","yoked"))
dds*P**u**n**c**h* &lt; −*f**a**c**t**o**r*(*d**d**s*Punch,
levels=c("DG","CA1", "CA3"))

1.4 Differential expression analysi
-----------------------------------

dds &lt;- DESeq(dds)

general deseq
=============

res &lt;- results(dds, independentFiltering = F) resOrdered &lt;-
res\[order(res$padj),\] summary(res) head(resOrdered,10)
sum(res*p**a**d**j* &lt; 0.1, *n**a*.*r**m* = *T**R**U**E*)*r**e**s*05 &lt; −*r**e**s**u**l**t**s*(*d**d**s*, *a**l**p**h**a* = 0.05)*s**u**m**m**a**r**y*(*r**e**s*05)*t**a**b**l**e*(*r**e**s*05padj
&lt; .05) sum(res05$padj &lt; 0.05, na.rm=TRUE)

1.5 exploring and reporting results
-----------------------------------

plotMA(res, main="plotMA")

resMLE &lt;- results(dds) head(resMLE, 4)

hist(res*p**v**a**l**u**e*\[*r**e**s*baseMean &gt; 1\], breaks=0:20/20,
col="grey50", border="white")

plotCounts(dds,
gene=which.min(res$padj), intgroup="Group") plotCounts(dds, gene=which.min(res$padj),
intgroup="Punch")

respadj &lt;- as.data.frame(res$padj) head(respadj)

1.5 more info
-------------

mcols(res)$description

for variance stablized gene expression and log transformed data
---------------------------------------------------------------

rld &lt;- rlog(dds, blind=FALSE) vsd &lt;-
varianceStabilizingTransformation(dds, blind=FALSE) vsd.fast &lt;-
vst(dds, blind=FALSE) head(assay(rld), 3)

\`\`\`
======

pca plot
--------

\`\`\`{r pca}
=============

pcaData &lt;- plotPCA(rld, intgroup = c( "Group", "Punch"),
returnData=TRUE) pcaData percentVar &lt;- round(100 \* attr(pcaData,
"percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Group, shape = Punch)) +
geom\_point(size=3) + xlab(paste0("PC1: ",percentVar\[1\],"% variance"))
+ ylab(paste0("PC2: ",percentVar\[2\],"% variance")) + coord\_fixed()
\#\`\`\`

`{r heatmap} library("genefilter") library("pheatmap") topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25) mat <- assay(rld)[ topVarGenes, ] mat <- mat - rowMeans(mat) df <- as.data.frame(colData(rld)[,c("Group", "Punch")]) pheatmap(mat) pheatmap(mat, show_colnames=F, show_rownames = T, annotation_col=df) #`
====================================================================================================================================================================================================================================================================================================================================

Session Info
------------

    sessionInfo()

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: OS X 10.10.5 (Yosemite)
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] DESeq2_1.14.0              SummarizedExperiment_1.4.0
    ##  [3] Biobase_2.34.0             GenomicRanges_1.26.1      
    ##  [5] GenomeInfoDb_1.10.1        IRanges_2.8.0             
    ##  [7] S4Vectors_0.12.0           BiocGenerics_0.20.0       
    ##  [9] cowplot_0.7.0              gplots_3.0.1              
    ## [11] magrittr_1.5               ggplot2_2.1.0             
    ## [13] reshape2_1.4.2             plyr_1.8.4                
    ## [15] dplyr_0.5.0                tidyr_0.6.0               
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] locfit_1.5-9.1       Rcpp_0.12.7          lattice_0.20-34     
    ##  [4] gtools_3.5.0         assertthat_0.1       rprojroot_1.2       
    ##  [7] digest_0.6.11        R6_2.2.0             chron_2.3-47        
    ## [10] backports_1.0.5      acepack_1.4.1        RSQLite_1.0.0       
    ## [13] evaluate_0.10        zlibbioc_1.20.0      data.table_1.9.6    
    ## [16] annotate_1.52.0      gdata_2.17.0         rpart_4.1-10        
    ## [19] Matrix_1.2-7.1       rmarkdown_1.3        splines_3.3.1       
    ## [22] BiocParallel_1.8.1   geneplotter_1.52.0   stringr_1.1.0       
    ## [25] foreign_0.8-67       RCurl_1.95-4.8       munsell_0.4.3       
    ## [28] htmltools_0.3.5      nnet_7.3-12          tibble_1.2          
    ## [31] gridExtra_2.2.1      htmlTable_1.7        Hmisc_4.0-0         
    ## [34] XML_3.98-1.4         bitops_1.0-6         grid_3.3.1          
    ## [37] xtable_1.8-2         gtable_0.2.0         DBI_0.5-1           
    ## [40] scales_0.4.0         KernSmooth_2.23-15   stringi_1.1.2       
    ## [43] XVector_0.14.0       genefilter_1.56.0    latticeExtra_0.6-28 
    ## [46] Formula_1.2-1        RColorBrewer_1.1-2   tools_3.3.1         
    ## [49] survival_2.40-1      yaml_2.1.14          AnnotationDbi_1.36.0
    ## [52] colorspace_1.2-7     cluster_2.0.5        caTools_1.17.1      
    ## [55] knitr_1.15.1
