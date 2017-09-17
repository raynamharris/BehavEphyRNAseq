Loading the data
----------------

In the summer of 2016, I processed a bunch of hippocampal tissue samples
from WT and FMR1-KO mice that were trained in an active place avoidance
task.

This data was added the the epically large sample collection database
contained in the two files "animals.csv" and "punches.csv" which
provided a detailed account of all animals processed and all tissue
samples collected. Then, I tidy the dataframe a little bit to get it
prepared for RNAseq sample submision.

    #install.packages("tidyr", dependencies=TRUE)
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("DESeq2")
    library("tidyr") 
    library("dplyr") ## for filtering and selecting rows
    library("plyr")  ## for renmaing factors
    library("reshape2") ##  for melting dataframe
    library("ggplot2") ## for awesome plots!
    library("magrittr") ## to use the weird pipe
    library("gplots") ##for making awesome plots
    library("cowplot") ## for some easy to use themes
    library("DESeq2") ## for differnetial gene expression profiling
    library("ggrepel") ## for labeling volcano plot

    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../results/fmr1/')

    source("functions_behavior.R")
    source("figureoptions.R")

Read the data

    # this starts with data genearated from code described in 01a_behavdatawrangling.R 
    behavior <- read.csv('../data/behavior/fmr1.csv' , header = T)

    # this starts with data genearated from code described in 02a_tidysamples.R and 02b_KallistoGather.Rmd
    colData <- read.csv('../data/rnaseq/fmr1ColData.csv')
    rownames(colData) <- colData$RNAseqID
    countData <-  read.csv('../data/rnaseq/fmr1CountData.csv', check.names = F, row.names = 1)

Behavior analysis
-----------------

    ## sert and revalue factors
    behavior$APA <- factor(behavior$APA, levels = c("control", "consistent", "conflict"))
    behavior$Genotype <- factor(behavior$Genotype, levels = c("WT", "FMR1KO"))
    behavior$Time1stEntrLog <- log(behavior$Time1stEntr) 

    ## behavior summary time
    behavior <- na.omit(behavior)
    behaviorsummaryNum <- dplyr::summarise(group_by(behavior, APA, Genotype, TrainSessionComboNum), m = mean(NumEntrances), se = sd(NumEntrances)/sqrt(length(NumEntrances)))
    behaviorsummaryNum

    ## Source: local data frame [54 x 5]
    ## Groups: APA, Genotype [?]
    ## 
    ##        APA Genotype TrainSessionComboNum        m        se
    ##     <fctr>   <fctr>                <int>    <dbl>     <dbl>
    ## 1  control       WT                    1 26.37500 2.0610114
    ## 2  control       WT                    2 12.00000 1.5118579
    ## 3  control       WT                    3 13.44444 0.8012336
    ## 4  control       WT                    4 17.57143 1.2883571
    ## 5  control       WT                    5 17.83333 2.8333333
    ## 6  control       WT                    6 17.22222 1.4979410
    ## 7  control       WT                    7 15.77778 1.4021588
    ## 8  control       WT                    8 13.75000 1.5323884
    ## 9  control       WT                    9 20.58333 1.9087452
    ## 10 control   FMR1KO                    1 28.57143 2.8440175
    ## # ... with 44 more rows

    numentrance <- ggplot(behaviorsummaryNum, aes(x=, TrainSessionComboNum, y=m, color=Genotype)) + 
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=Genotype), width=.1) +
        geom_point(size = 2) +
        geom_line() +
        facet_wrap(~APA) +
        scale_y_continuous(name="Number of entrances") +
        scale_x_continuous(name = NULL, 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels=c("1" = "Hab.", "2" = "T1", "3" = "T2", 
                                    "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                    "7" = "T5/C2", "8" = "T6/C3", "9"= "Reten.")) +
      #theme_cowplot(font_size = 14, line_size = 1) +
      background_grid(major = "y", minor = "y") +
      #scale_color_manual(values = colorvalAPA) + 
      theme(legend.position="top") + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    numentrance

![](../results/fmr1/behaviordatawrangle-1.png)

    pdf(file="../results/fmr1/behaviordatawrangle-1.pdf", width=7.5, height=3)
    plot(numentrance)
    dev.off()

    ## quartz_off_screen 
    ##                 2

RNAsequencing analysis
----------------------

    ## remove outliers
    colData <- colData %>% filter(RNAseqID != "16-125B", RNAseqID != "16-123B") #set the coldata to be the countbygene df

    ## colData and countData must contain the exact same sample. I'll use the next three lines to make that happen
    savecols <- as.character(colData$RNAseqID) #select the sample name column that corresponds to row names
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% dplyr::select(one_of(savecols)) # select just the columns that match the samples in colData

DESeq Analysis
--------------

Now, I'll look for differential gene expression between the FMR1-KO and
WT mice. This analysis was developed by reading the DESEq manual. In
many place, I try to provide the chapter where these steps are described
in more details.

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 134 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 
    ## out of 16871 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 2, 0.012% 
    ## LFC < 0 (down)   : 9, 0.053% 
    ## outliers [1]     : 45, 0.27% 
    ## low counts [2]   : 0, 0% 
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    ## log2 fold change (MAP): Genotype FMR1 vs WT 
    ## Wald test p-value: Genotype FMR1 vs WT 
    ## DataFrame with 10 rows and 6 columns
    ##            baseMean log2FoldChange      lfcSE       stat       pvalue
    ##           <numeric>      <numeric>  <numeric>  <numeric>    <numeric>
    ## Ccnd2      46.67797     -1.5818867 0.14281653 -11.076356 1.633885e-28
    ## Fmr1       93.47210     -1.0773829 0.13916618  -7.741701 9.809573e-15
    ## Arel1     886.18923      0.3785902 0.07530677   5.027306 4.974185e-07
    ## Cry2      508.25058      0.4323033 0.09101104   4.750010 2.034066e-06
    ## Kcnt1      72.55067     -0.6511921 0.14092433  -4.620864 3.821461e-06
    ## Cacna1g   189.89551     -0.5805430 0.13527607  -4.291543 1.774361e-05
    ## Mtus1      18.57150     -0.5788129 0.13638252  -4.244040 2.195308e-05
    ## Serpina3n 133.33547     -0.4908262 0.11598815  -4.231693 2.319392e-05
    ## Car4       51.78762     -0.5900355 0.14407291  -4.095395 4.214490e-05
    ## Ephx1      39.42961     -0.5889576 0.14354308  -4.103002 4.078234e-05
    ##                   padj
    ##              <numeric>
    ## Ccnd2     2.752932e-24
    ## Fmr1      8.264075e-11
    ## Arel1     2.793668e-03
    ## Cry2      8.567995e-03
    ## Kcnt1     1.287756e-02
    ## Cacna1g   4.884930e-02
    ## Mtus1     4.884930e-02
    ## Serpina3n 4.884930e-02
    ## Car4      7.100994e-02
    ## Ephx1     7.100994e-02

    ##                16-116B  16-117D  16-118B  16-118D  16-119B  16-119D
    ## 0610007P14Rik 5.386474 5.537126 5.120084 5.377359 5.387834 5.288434
    ## 0610009B22Rik 4.292298 4.325869 4.097429 4.305376 4.210168 4.168866
    ## 0610009L18Rik 2.038900 1.956492 1.974084 2.106994 1.976300 1.962694
    ##                16-120B  16-120D  16-122B  16-122D  16-123D  16-124D
    ## 0610007P14Rik 5.176035 5.357994 5.412306 5.444757 5.349111 5.301910
    ## 0610009B22Rik 4.172616 4.014141 4.495170 4.280107 4.219794 4.293576
    ## 0610009L18Rik 1.987914 1.949666 2.051807 2.016388 2.021595 1.992984
    ##                16-125D  16-126B
    ## 0610007P14Rik 5.284383 5.336251
    ## 0610009B22Rik 4.139385 4.231975
    ## 0610009L18Rik 1.950443 2.173058

Data viz
--------

![](../results/fmr1/plots-1.png)![](../results/fmr1/plots-2.png)![](../results/fmr1/plots-3.png)

Volcano plot
------------

    res <- results(dds, contrast =c("Genotype", "FMR1", "WT"), independentFiltering = F)
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
    with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
    with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))

![](../results/fmr1/volcanoplot-1.png)

    resOrdered <- res[order(res$padj),]
    head(resOrdered)

    ## log2 fold change (MAP): Genotype FMR1 vs WT 
    ## Wald test p-value: Genotype FMR1 vs WT 
    ## DataFrame with 6 rows and 6 columns
    ##          baseMean log2FoldChange      lfcSE       stat       pvalue
    ##         <numeric>      <numeric>  <numeric>  <numeric>    <numeric>
    ## Ccnd2    46.67797     -1.5818867 0.14281653 -11.076356 1.633885e-28
    ## Fmr1     93.47210     -1.0773829 0.13916618  -7.741701 9.809573e-15
    ## Arel1   886.18923      0.3785902 0.07530677   5.027306 4.974185e-07
    ## Cry2    508.25058      0.4323033 0.09101104   4.750010 2.034066e-06
    ## Kcnt1    72.55067     -0.6511921 0.14092433  -4.620864 3.821461e-06
    ## Cacna1g 189.89551     -0.5805430 0.13527607  -4.291543 1.774361e-05
    ##                 padj
    ##            <numeric>
    ## Ccnd2   2.752932e-24
    ## Fmr1    8.264075e-11
    ## Arel1   2.793668e-03
    ## Cry2    8.567995e-03
    ## Kcnt1   1.287756e-02
    ## Cacna1g 4.884930e-02

    data <- data.frame(gene = row.names(res), pvalue = -log10(res$padj), lfc = res$log2FoldChange)
    data <- na.omit(data)
    head(data)

    ##            gene       pvalue         lfc
    ## 1 0610007P14Rik 0.0012632616 -0.04815661
    ## 2 0610009B22Rik 0.0012632616 -0.16802970
    ## 3 0610009L18Rik 0.0012632616 -0.10012651
    ## 4 0610009O20Rik 0.0008721801  0.02837217
    ## 5 0610010F05Rik 0.0012632616  0.08569528
    ## 6 0610010K14Rik 0.0012632616 -0.10993716

    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                            yes = "WT", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                        yes = "FMR1", 
                                        no = "none")))
    top_labelled <- top_n(data, n = 9, wt = pvalue)

    # Color corresponds to fold change directionality
    colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend 
      scale_color_manual(values = c("FMR1" = "#7570b3",
                                    "WT" = "#d95f02", 
                                    "none" = "#bdbdbd")) + theme(panel.grid.minor=element_blank(),
               panel.grid.major=element_blank()) + 
      theme(axis.title.x = element_blank())+ 
      theme(axis.title.y = element_blank()) + 
      geom_text_repel(data = top_labelled, 
                              mapping = aes(label = gene), 
                              size = 3,
                              fontface = 'bold', 
                              color = 'black',
                              box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.5, "lines"))

    colored

![](../results/fmr1/volcanoplot-2.png)

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
    ##  [1] ggrepel_0.6.5              DESeq2_1.14.1             
    ##  [3] SummarizedExperiment_1.4.0 Biobase_2.34.0            
    ##  [5] GenomicRanges_1.26.3       GenomeInfoDb_1.10.3       
    ##  [7] IRanges_2.8.1              S4Vectors_0.12.1          
    ##  [9] BiocGenerics_0.20.0        cowplot_0.7.0             
    ## [11] gplots_3.0.1               magrittr_1.5              
    ## [13] ggplot2_2.2.1              reshape2_1.4.2            
    ## [15] plyr_1.8.4                 dplyr_0.5.0               
    ## [17] tidyr_0.6.1               
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] splines_3.3.1        gtools_3.5.0         Formula_1.2-1       
    ##  [4] assertthat_0.1       latticeExtra_0.6-28  yaml_2.1.14         
    ##  [7] RSQLite_1.1-2        backports_1.0.5      lattice_0.20-34     
    ## [10] digest_0.6.12        RColorBrewer_1.1-2   XVector_0.14.0      
    ## [13] checkmate_1.8.2      colorspace_1.3-2     htmltools_0.3.5     
    ## [16] Matrix_1.2-8         XML_3.98-1.5         genefilter_1.56.0   
    ## [19] zlibbioc_1.20.0      xtable_1.8-2         scales_0.4.1        
    ## [22] gdata_2.17.0         BiocParallel_1.8.1   htmlTable_1.9       
    ## [25] tibble_1.2           annotate_1.52.1      nnet_7.3-12         
    ## [28] lazyeval_0.2.0       survival_2.40-1      memoise_1.0.0       
    ## [31] evaluate_0.10        foreign_0.8-67       tools_3.3.1         
    ## [34] data.table_1.10.0    stringr_1.2.0        munsell_0.4.3       
    ## [37] locfit_1.5-9.1       cluster_2.0.5        AnnotationDbi_1.36.2
    ## [40] caTools_1.17.1       grid_3.3.1           RCurl_1.95-4.8      
    ## [43] htmlwidgets_0.8      bitops_1.0-6         base64enc_0.1-3     
    ## [46] labeling_0.3         rmarkdown_1.3        gtable_0.2.0        
    ## [49] DBI_0.6              R6_2.2.0             gridExtra_2.2.1     
    ## [52] knitr_1.15.1         Hmisc_4.0-2          rprojroot_1.2       
    ## [55] KernSmooth_2.23-15   stringi_1.1.2        Rcpp_0.12.9         
    ## [58] geneplotter_1.52.0   rpart_4.1-10         acepack_1.4.1
