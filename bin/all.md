The Overall Research Question
-----------------------------

The research project was designed to understand how experience shapes
the brain. In particular, we are looking at learned avoidance behavhior.
How do animals change their behavior to avoidance an unpleasant
experience?

RNAseq samples
--------------

In the summers of 2015 and 2016, I processed a bunch of hippocampal
tissue samples from 59 mice. Most mice were trained in an active place
avoidance task or used as yoked controls; however, a few animals were
taken straight from the home cage. This data has been cleaned using a
different script `02a_punches.R`.

THis output provides a summary of the samples and the various factors
that descibe them.

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

Merge transcipts counts to gene counts
--------------------------------------

Kallisto is cool because it does give you transcript level stuff, but
right now I think I have more power if I just look at gene level stuff.
I'll go back and look at transcripts if I want to.

    ## Joining, by = "id"
    ## Joining, by = "id"

    ##     142C_CA1           142C_DG          143A-CA3-1         143A-DG-1      
    ##  Min.   :    0.00   Min.   :   0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:   0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    2.00   Median :   2.00   Median :    2.00   Median :    2.0  
    ##  Mean   :   24.02   Mean   :  21.98   Mean   :   23.66   Mean   :   25.9  
    ##  3rd Qu.:   14.00   3rd Qu.:  15.00   3rd Qu.:   12.00   3rd Qu.:   12.0  
    ##  Max.   :19994.00   Max.   :9036.00   Max.   :13308.00   Max.   :25544.0  
    ##    143B-CA1-1        143B-DG-1           143C_CA1          143C_DG        
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    2.0   Median :    2.00   Median :    2.0   Median :    2.00  
    ##  Mean   :   26.4   Mean   :   24.46   Mean   :   23.4   Mean   :   23.57  
    ##  3rd Qu.:   11.0   3rd Qu.:   13.00   3rd Qu.:   14.0   3rd Qu.:   15.00  
    ##  Max.   :27782.0   Max.   :18912.00   Max.   :13133.0   Max.   :18685.00  
    ##    143C-CA1-1         143D-CA1-3         143D-DG-3       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    1.00   Median :    2.00  
    ##  Mean   :   24.73   Mean   :   24.64   Mean   :   23.75  
    ##  3rd Qu.:   12.00   3rd Qu.:   12.00   3rd Qu.:   14.00  
    ##  Max.   :18125.00   Max.   :19703.00   Max.   :16880.00  
    ##    144A-CA1-2         144A-CA3-2         144A-DG-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    1.00   Median :    2.00  
    ##  Mean   :   26.07   Mean   :   28.67   Mean   :   24.35  
    ##  3rd Qu.:   11.00   3rd Qu.:    8.00   3rd Qu.:   14.00  
    ##  Max.   :24613.00   Max.   :38260.00   Max.   :19211.00  
    ##    144B-CA1-1         144B-CA3-1         144C-CA1-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    1.00   Median :    2.00  
    ##  Mean   :   25.35   Mean   :   24.93   Mean   :   24.81  
    ##  3rd Qu.:   12.00   3rd Qu.:   11.00   3rd Qu.:   12.00  
    ##  Max.   :23010.00   Max.   :20084.00   Max.   :19363.00  
    ##    144C-CA3-2         144C-DG-2          144D-CA3-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    1.00   Median :    2.00   Median :    2.00  
    ##  Mean   :   27.33   Mean   :   23.82   Mean   :   27.66  
    ##  3rd Qu.:    9.00   3rd Qu.:   14.00   3rd Qu.:    9.00  
    ##  Max.   :30938.00   Max.   :13913.00   Max.   :35044.00  
    ##    144D-DG-2          145A-CA1-2         145A-CA3-2     
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    2.00   Median :    2.00   Median :    0.0  
    ##  Mean   :   23.88   Mean   :   26.57   Mean   :   26.3  
    ##  3rd Qu.:   13.00   3rd Qu.:   11.00   3rd Qu.:    9.0  
    ##  Max.   :15967.00   Max.   :27499.00   Max.   :27390.0  
    ##    145A-DG-2          145B-CA1-1         145B-DG-1       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    1.00   Median :    2.00  
    ##  Mean   :   23.62   Mean   :   26.78   Mean   :   23.51  
    ##  3rd Qu.:   14.00   3rd Qu.:   11.00   3rd Qu.:   14.00  
    ##  Max.   :14783.00   Max.   :31778.00   Max.   :18049.00  
    ##    146A-CA1-2         146A-CA3-2         146A-DG-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    1.00   Median :    1.00   Median :    2.00  
    ##  Mean   :   26.42   Mean   :   29.41   Mean   :   24.69  
    ##  3rd Qu.:   10.00   3rd Qu.:    8.00   3rd Qu.:   13.00  
    ##  Max.   :25334.00   Max.   :43024.00   Max.   :19323.00  
    ##    146B-CA1-2         146B-CA3-2         146B-DG-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    2.00   Median :    0.00  
    ##  Mean   :   23.34   Mean   :   24.81   Mean   :   20.54  
    ##  3rd Qu.:   14.00   3rd Qu.:   11.00   3rd Qu.:    8.00  
    ##  Max.   :12770.00   Max.   :17434.00   Max.   :30815.00  
    ##    146C-CA1-4         146C-CA3-4         146C-DG-4       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    0.00   Median :    1.00  
    ##  Mean   :   24.66   Mean   :   27.42   Mean   :   22.69  
    ##  3rd Qu.:   13.00   3rd Qu.:    9.00   3rd Qu.:   14.00  
    ##  Max.   :20521.00   Max.   :26486.00   Max.   :10508.00  
    ##    146D-CA1-3         146D-CA3-3         146D-DG-3      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    0.00   Median :    2.00   Median :    0.0  
    ##  Mean   :   23.19   Mean   :   24.67   Mean   :   20.1  
    ##  3rd Qu.:    8.00   3rd Qu.:   12.00   3rd Qu.:    4.0  
    ##  Max.   :22343.00   Max.   :21103.00   Max.   :29508.0  
    ##    147-CA1-4          147-CA3-4           147-DG-4       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   22.36   Mean   :   23.63   Mean   :   20.58  
    ##  3rd Qu.:    4.00   3rd Qu.:    8.00   3rd Qu.:    2.00  
    ##  Max.   :25916.00   Max.   :18064.00   Max.   :21507.00  
    ##    147C-CA1-3         147C-CA3-3         147C-DG-3       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    2.00   Median :    3.00  
    ##  Mean   :   25.15   Mean   :   28.39   Mean   :   24.47  
    ##  3rd Qu.:   13.00   3rd Qu.:   10.00   3rd Qu.:   14.00  
    ##  Max.   :21638.00   Max.   :38426.00   Max.   :18962.00  
    ##    147D-CA3-1         147D-DG-1          148-CA1-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    2.00   Median :    1.00  
    ##  Mean   :   27.19   Mean   :   23.61   Mean   :   24.83  
    ##  3rd Qu.:   10.00   3rd Qu.:   14.00   3rd Qu.:   12.00  
    ##  Max.   :30897.00   Max.   :15216.00   Max.   :22177.00  
    ##    148-CA3-2          148-DG-2          148A-CA1-3      
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    1.0   Median :    1.00   Median :    2.00  
    ##  Mean   :   24.3   Mean   :   23.05   Mean   :   24.93  
    ##  3rd Qu.:   11.0   3rd Qu.:   14.00   3rd Qu.:   12.00  
    ##  Max.   :16496.0   Max.   :12329.00   Max.   :18801.00  
    ##    148A-CA3-3         148A-DG-3          148B-CA1-4      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    1.00   Median :    3.00   Median :    0.00  
    ##  Mean   :   25.43   Mean   :   23.54   Mean   :   25.25  
    ##  3rd Qu.:   10.00   3rd Qu.:   15.00   3rd Qu.:    8.00  
    ##  Max.   :19835.00   Max.   :15700.00   Max.   :54348.00  
    ##    148B-CA3-4         148B-DG-4           16-116B        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    2.00   Median :    2.00  
    ##  Mean   :   25.15   Mean   :   25.11   Mean   :   22.18  
    ##  3rd Qu.:   11.00   3rd Qu.:   12.00   3rd Qu.:   13.00  
    ##  Max.   :19770.00   Max.   :21584.00   Max.   :11884.00  
    ##     16-116D            16-117D            16-118B        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    2.00   Median :    2.00  
    ##  Mean   :   20.66   Mean   :   22.48   Mean   :   22.37  
    ##  3rd Qu.:    3.00   3rd Qu.:   13.00   3rd Qu.:   13.00  
    ##  Max.   :89831.00   Max.   :11514.00   Max.   :10517.00  
    ##     16-118D            16-119B            16-119D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    2.00   Median :    2.00  
    ##  Mean   :   23.38   Mean   :   24.11   Mean   :   22.43  
    ##  3rd Qu.:   13.00   3rd Qu.:   12.00   3rd Qu.:   13.00  
    ##  Max.   :12942.00   Max.   :14617.00   Max.   :10258.00  
    ##     16-120B           16-120D            16-122B        
    ##  Min.   :   0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:   0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :   2.00   Median :    2.00   Median :    2.00  
    ##  Mean   :  22.06   Mean   :   22.69   Mean   :   23.94  
    ##  3rd Qu.:  13.00   3rd Qu.:   13.00   3rd Qu.:   12.00  
    ##  Max.   :9640.00   Max.   :12280.00   Max.   :14795.00  
    ##     16-122D            16-123B            16-123D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    0.00   Median :    2.00  
    ##  Mean   :   23.63   Mean   :   21.45   Mean   :   22.99  
    ##  3rd Qu.:   13.00   3rd Qu.:   10.00   3rd Qu.:   13.00  
    ##  Max.   :13117.00   Max.   :24154.00   Max.   :12331.00  
    ##     16-124D            16-125B            16-125D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    2.00   Median :    1.00   Median :    2.00  
    ##  Mean   :   22.54   Mean   :   14.37   Mean   :   22.79  
    ##  3rd Qu.:   13.00   3rd Qu.:    5.00   3rd Qu.:   13.00  
    ##  Max.   :10451.00   Max.   :78181.00   Max.   :11716.00  
    ##     16-126B        
    ##  Min.   :    0.00  
    ##  1st Qu.:    0.00  
    ##  Median :    2.00  
    ##  Mean   :   24.34  
    ##  3rd Qu.:   12.00  
    ##  Max.   :16647.00

DESeq Analysis
--------------

Now, I'll look for differential gene expression between of the various
factors. This analysis was developed by reading the DESEq manual. In
many place, I try to provide the chapter where these steps are described
in more details.

    # 1.3.3 Count matrix input ----
    countData <- countbygene 
    colData <- tidysamples %>%
      arrange(RNAseqID) # needs to be in same order a countData
    head(countData)

    ##               142C_CA1 142C_DG 143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1
    ## 0610007P14Rik       60      44         42        56         30        24
    ## 0610009B22Rik       31      14         12        17         10         5
    ## 0610009L18Rik        0       0          4         9         10         8
    ## 0610009O20Rik      128     112         85       185         44        72
    ## 0610010F05Rik       14      16         18        19          7        15
    ## 0610010K14Rik        2       1          2         7          1         2
    ##               143C_CA1 143C_DG 143C-CA1-1 143D-CA1-3 143D-DG-3 144A-CA1-2
    ## 0610007P14Rik       41      16         19         14        22         40
    ## 0610009B22Rik       22      10         10          0         0         15
    ## 0610009L18Rik        6       0          2          0         2          9
    ## 0610009O20Rik      170     109         76         25        38         95
    ## 0610010F05Rik       43      16          7          5         8         15
    ## 0610010K14Rik        6       3          2          2         1          3
    ##               144A-CA3-2 144A-DG-2 144B-CA1-1 144B-CA3-1 144C-CA1-2
    ## 0610007P14Rik         10        40         36         17         32
    ## 0610009B22Rik          4         4          7          4          8
    ## 0610009L18Rik          5         0          1          4          2
    ## 0610009O20Rik         20        82         49         24         89
    ## 0610010F05Rik          2        17          9          3         12
    ## 0610010K14Rik          1         2          4          1          3
    ##               144C-CA3-2 144C-DG-2 144D-CA3-2 144D-DG-2 145A-CA1-2
    ## 0610007P14Rik         14        25         22        75         66
    ## 0610009B22Rik         12         8          6        12         18
    ## 0610009L18Rik          4         6          9        13         21
    ## 0610009O20Rik         48        85         46       151         96
    ## 0610010F05Rik          7        12         14        18         22
    ## 0610010K14Rik          1         2          2         6          3
    ##               145A-CA3-2 145A-DG-2 145B-CA1-1 145B-DG-1 146A-CA1-2
    ## 0610007P14Rik          2        20         26        10         14
    ## 0610009B22Rik          3         7          8         5          8
    ## 0610009L18Rik          0         1          3         0          8
    ## 0610009O20Rik         12        52        124        46         68
    ## 0610010F05Rik          1         6          4         5          9
    ## 0610010K14Rik          0         2          2         2          3
    ##               146A-CA3-2 146A-DG-2 146B-CA1-2 146B-CA3-2 146B-DG-2
    ## 0610007P14Rik         44        12          6         16         0
    ## 0610009B22Rik          4         3          3         22         0
    ## 0610009L18Rik          9         6          0          2         0
    ## 0610009O20Rik        102        44         16         46         2
    ## 0610010F05Rik         10         6          7         14         1
    ## 0610010K14Rik          1         2          1          1         1
    ##               146C-CA1-4 146C-CA3-4 146C-DG-4 146D-CA1-3 146D-CA3-3
    ## 0610007P14Rik         19          2        11          6         44
    ## 0610009B22Rik          4          4         3          2         12
    ## 0610009L18Rik          9          0         0          0          7
    ## 0610009O20Rik         31          3        10          1         83
    ## 0610010F05Rik          5          1         4          1         18
    ## 0610010K14Rik          1          0         0          1          1
    ##               146D-DG-3 147-CA1-4 147-CA3-4 147-DG-4 147C-CA1-3 147C-CA3-3
    ## 0610007P14Rik         3         1        46        0         34         82
    ## 0610009B22Rik         0         1        10        4          8         15
    ## 0610009L18Rik         0         0         0        0          2         11
    ## 0610009O20Rik         3         0        29        0         58        191
    ## 0610010F05Rik         0         2         0        2         19         41
    ## 0610010K14Rik         0         0         1        0          4          3
    ##               147C-DG-3 147D-CA3-1 147D-DG-1 148-CA1-2 148-CA3-2 148-DG-2
    ## 0610007P14Rik        41         40       152        24        29       21
    ## 0610009B22Rik        20         20        52         3        23        8
    ## 0610009L18Rik         3          9        67        13         0        0
    ## 0610009O20Rik       145         88       377        27        91      164
    ## 0610010F05Rik        23         28        56         9        19        7
    ## 0610010K14Rik         2          0        16         0         2        1
    ##               148A-CA1-3 148A-CA3-3 148A-DG-3 148B-CA1-4 148B-CA3-4
    ## 0610007P14Rik         68         28        52          2         61
    ## 0610009B22Rik         30         12         8          0         22
    ## 0610009L18Rik         16          7        11          0         11
    ## 0610009O20Rik        162         65       226          0         70
    ## 0610010F05Rik         25         19        22          0         22
    ## 0610010K14Rik          5          4         2          0          4
    ##               148B-DG-4 16-116B 16-116D 16-117D 16-118B 16-118D 16-119B
    ## 0610007P14Rik         8      19       2      19      13      21      25
    ## 0610009B22Rik         1      10       0       8       6      12      10
    ## 0610009L18Rik         1       5       0       0       2      10       2
    ## 0610009O20Rik        18      40       3      35      86      28      62
    ## 0610010F05Rik         3       8       0       9      16      10      10
    ## 0610010K14Rik         1       3       0       1       1       1       3
    ##               16-119D 16-120B 16-120D 16-122B 16-122D 16-123B 16-123D
    ## 0610007P14Rik      19      16      19      28      35       4      22
    ## 0610009B22Rik       8       8       2      25      16       0       9
    ## 0610009L18Rik       1       3       0       8       6       0       5
    ## 0610009O20Rik      60      93      94      69      89       4      62
    ## 0610010F05Rik      12      13      11      14      15       3      12
    ## 0610010K14Rik       2       2       0       2       3       0       1
    ##               16-124D 16-125B 16-125D 16-126B
    ## 0610007P14Rik      20       2      15      24
    ## 0610009B22Rik      12       0       5      11
    ## 0610009L18Rik       3       0       0      18
    ## 0610009O20Rik      40      13      61      92
    ## 0610010F05Rik      12       1       9       9
    ## 0610010K14Rik       3       0       1       4

    head(colData)

    ##    RNAseqID  Mouse year Genotype jobnumber Punch    Group Conflict  APA
    ## 1 100-CA1-1 15-100 2015       WT   JA16444   CA1 homecage     <NA> <NA>
    ## 2 100-CA1-2 15-100 2015       WT   JA16444   CA1 homecage     <NA> <NA>
    ## 3 100-CA1-3 15-100 2015       WT   JA16444   CA1 homecage     <NA> <NA>
    ## 4 100-CA3-1 15-100 2015       WT   JA16444   CA3 homecage     <NA> <NA>
    ## 5 100-CA3-4 15-100 2015       WT   JA16444   CA3 homecage     <NA> <NA>
    ## 6  100-DG-2 15-100 2015       WT   JA16444    DG homecage     <NA> <NA>
    ##        method   dodgy  daytime Slice    Date
    ## 1 homogenized allgood norecord     1 9/28/15
    ## 2 homogenized allgood norecord     2 9/28/15
    ## 3 homogenized allgood norecord     3 9/28/15
    ## 4 homogenized allgood norecord     1 9/28/15
    ## 5 homogenized allgood norecord     4 9/28/15
    ## 6 homogenized allgood norecord     2 9/28/15

    ## eh, I have a few too many samples in the countDataFrame
    savecols <- colnames(countData)
    colData <- colData %>%
      filter(RNAseqID %in% savecols) %>% droplevels()

    ## remove outliers as per 1st pca

    ## removeing outliers
    colData <- colData %>%
      filter(RNAseqID != "146C-CA3-4", RNAseqID != "16-116D") # needs to be in same order a countData
    savecols <- as.character(colData$RNAseqID) #select the sample name column that corresponds to row names
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% select(one_of(savecols)) # select just the columns that match the samples in colData

    ## remove genes with total counts across all samples < 2
    countData[countData < 2] <- 0

    ## write files for easy import
    write.csv(countData, "../data/rnaseq/countData.csv", row.names=FALSE)
    write.csv(colData, "../data/rnaseq/colData.csv", row.names=FALSE)

\`\`\`{r DESeq1}
================

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
levels=c("homecage","control","consistent", "conflict"))
dds*P**u**n**c**h* &lt; −*f**a**c**t**o**r*(*d**d**s*Punch,
levels=c("DG","CA1", "CA3"))

1.4 Differential expression analysi
-----------------------------------

dds &lt;- DESeq(dds)

general deseq
=============

res &lt;- results(dds, independentFiltering = F) resOrdered &lt;-
res\[order(res$padj),\] head(resOrdered,10)

sum(res*p**a**d**j* &lt; 0.1, *n**a*.*r**m* = *T**R**U**E*)*r**e**s*05 &lt; −*r**e**s**u**l**t**s*(*d**d**s*, *a**l**p**h**a* = 0.05)*t**a**b**l**e*(*r**e**s*05padj
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
vst(dds, blind=FALSE) head(assay(rld), 3) \`\`\`

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
    ## [13] evaluate_0.10        zlibbioc_1.20.0      lazyeval_0.2.0      
    ## [16] data.table_1.9.6     annotate_1.52.0      gdata_2.17.0        
    ## [19] rpart_4.1-10         Matrix_1.2-7.1       rmarkdown_1.3       
    ## [22] splines_3.3.1        BiocParallel_1.8.1   geneplotter_1.52.0  
    ## [25] stringr_1.1.0        foreign_0.8-67       RCurl_1.95-4.8      
    ## [28] munsell_0.4.3        htmltools_0.3.5      nnet_7.3-12         
    ## [31] tibble_1.2           gridExtra_2.2.1      htmlTable_1.7       
    ## [34] Hmisc_4.0-0          XML_3.98-1.4         bitops_1.0-6        
    ## [37] grid_3.3.1           xtable_1.8-2         gtable_0.2.0        
    ## [40] DBI_0.5-1            scales_0.4.0         KernSmooth_2.23-15  
    ## [43] stringi_1.1.2        XVector_0.14.0       genefilter_1.56.0   
    ## [46] latticeExtra_0.6-28  Formula_1.2-1        RColorBrewer_1.1-2  
    ## [49] tools_3.3.1          survival_2.40-1      yaml_2.1.14         
    ## [52] AnnotationDbi_1.36.0 colorspace_1.2-7     cluster_2.0.5       
    ## [55] caTools_1.17.1       knitr_1.15.1
