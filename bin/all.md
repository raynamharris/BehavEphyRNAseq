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

    ## remove genes with total counts across all samples < 2
    countData[countData < 2] <- 0

    ## differential gene expression
    dds <- DESeqDataSetFromMatrix(countData = countData,
    colData = colData,
    design = ~ Group + Punch)

    ## converting counts to integer mode

    dds

    ## class: DESeqDataSet 
    ## dim: 22485 72 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(72): 142C_CA1 142C_DG ... 16-125D 16-126B
    ## colData names(14): RNAseqID Mouse ... Slice Date

    ## 1.3.6 Pre-filtering
    dds <- dds[ rowSums(counts(dds)) > 1, ]

    ## 1.3.7 Note on factor levels
    dds$Group <- factor(dds$Group, levels=c("homecage","control","consistent", "conflict"))
    dds$Punch <- factor(dds$Punch, levels=c("DG","CA1", "CA3"))

    ## 1.4  Differential expression analysi
    dds <- DESeq(dds)

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 858 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    # general deseq
    res <- results(dds, independentFiltering = F)
    resOrdered <- res[order(res$padj),]
    summary(res)

    ## 
    ## out of 17958 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 3213, 18% 
    ## LFC < 0 (down)   : 4095, 23% 
    ## outliers [1]     : 276, 1.5% 
    ## low counts [2]   : 0, 0% 
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head(resOrdered,10)

    ## log2 fold change (MAP): Punch CA3 vs DG 
    ## Wald test p-value: Punch CA3 vs DG 
    ## DataFrame with 10 rows and 6 columns
    ##          baseMean log2FoldChange      lfcSE      stat        pvalue
    ##         <numeric>      <numeric>  <numeric> <numeric>     <numeric>
    ## Doc2b    389.1356      -6.478901 0.18814262 -34.43612 7.266237e-260
    ## Pitpnm2  132.4383      -3.279066 0.10397585 -31.53681 2.719445e-218
    ## C1ql3    248.7928      -6.992053 0.22526521 -31.03920 1.595686e-211
    ## Gnao1    161.1137       1.415813 0.05151084  27.48573 2.600102e-166
    ## Adcy1   2515.2328      -4.008200 0.15001460 -26.71873 2.851800e-157
    ## Lynx1    261.8223       2.578374 0.09677866  26.64197 2.217611e-156
    ## Fam163b  475.1924      -5.343691 0.20108893 -26.57377 1.364687e-155
    ## Syn2     321.9665       1.560201 0.05909040  26.40363 1.244809e-153
    ## Syngr1   138.9374       1.842224 0.07288544  25.27561 5.923733e-141
    ## Prkcg    573.9046       2.246629 0.08929613  25.15931 1.117839e-139
    ##                  padj
    ##             <numeric>
    ## Doc2b   1.287287e-255
    ## Pitpnm2 2.408884e-214
    ## C1ql3   9.423057e-208
    ## Gnao1   1.151585e-162
    ## Adcy1   1.010450e-153
    ## Lynx1   6.547865e-153
    ## Fam163b 3.453827e-152
    ## Syn2    2.756629e-150
    ## Syngr1  1.166054e-137
    ## Prkcg   1.980364e-136

    sum(res$padj < 0.1, na.rm = TRUE) 

    ## [1] 7308

    res05 <- results(dds, alpha=0.05)
    summary(res05) 

    ## 
    ## out of 17958 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)     : 3004, 17% 
    ## LFC < 0 (down)   : 3827, 21% 
    ## outliers [1]     : 276, 1.5% 
    ## low counts [2]   : 3689, 21% 
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    table(res05$padj < .05)

    ## 
    ## FALSE  TRUE 
    ##  7196  6831

    sum(res05$padj < 0.05, na.rm=TRUE)

    ## [1] 6831

    ## 1.5 exploring and reporting results

    plotMA(res, main="plotMA")

![](../results/all/DESeq-1.png)

    resMLE <- results(dds)
    head(resMLE, 4)

    ## log2 fold change (MAP): Punch CA3 vs DG 
    ## Wald test p-value: Punch CA3 vs DG 
    ## DataFrame with 4 rows and 6 columns
    ##                baseMean log2FoldChange     lfcSE       stat      pvalue
    ##               <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ## 0610007P14Rik 19.294373      0.3554797 0.1956735  1.8166984 0.069263309
    ## 0610009B22Rik  7.132269      0.8739272 0.2683269  3.2569501 0.001126162
    ## 0610009L18Rik  2.921439      0.5003623 0.5243815  0.9541953 0.339984797
    ## 0610009O20Rik 44.797784     -0.2198744 0.1831514 -1.2005062 0.229942816
    ##                      padj
    ##                 <numeric>
    ## 0610007P14Rik 0.126354933
    ## 0610009B22Rik 0.003309787
    ## 0610009L18Rik 0.463152266
    ## 0610009O20Rik 0.339638786

    hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")

![](../results/all/DESeq-2.png)

    plotCounts(dds, gene=which.min(res$padj), intgroup="Group")

![](../results/all/DESeq-3.png)

    plotCounts(dds, gene=which.min(res$padj), intgroup="Punch")

![](../results/all/DESeq-4.png)

    respadj <- as.data.frame(res$padj)
    head(respadj)

    ##      res$padj
    ## 1 0.155739152
    ## 2 0.004084153
    ## 3 0.565343595
    ## 4 0.416445197
    ## 5 0.109977342
    ## 6 0.283418977

    ## 1.5 more info
    mcols(res)$description

    ## [1] "mean of normalized counts for all samples"
    ## [2] "log2 fold change (MAP): Punch CA3 vs DG"  
    ## [3] "standard error: Punch CA3 vs DG"          
    ## [4] "Wald statistic: Punch CA3 vs DG"          
    ## [5] "Wald test p-value: Punch CA3 vs DG"       
    ## [6] "BH adjusted p-values"

    ## for variance stablized gene expression and log transformed data
    rld <- rlog(dds, blind=FALSE)
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    vsd.fast <- vst(dds, blind=FALSE)
    head(assay(rld), 3)

    ##                142C_CA1   142C_DG 143A-CA3-1 143A-DG-1 143B-CA1-1
    ## 0610007P14Rik 4.2677519 4.1863350   4.246738  4.097138   4.597665
    ## 0610009B22Rik 3.0045188 2.5678001   2.532332  2.473324   2.892293
    ## 0610009L18Rik 0.7135633 0.7258563   1.197163  1.350874   2.032248
    ##               143B-DG-1 143C_CA1   143C_DG 143C-CA1-1 143D-CA1-3 143D-DG-3
    ## 0610007P14Rik  4.168940 3.437211 3.2857395   3.823459   4.280930  4.813606
    ## 0610009B22Rik  2.283801 2.415943 2.4041845   2.673166   1.563376  1.578111
    ## 0610009L18Rik  1.768150 1.070580 0.7317905   1.100745   0.833524  1.392884
    ##               144A-CA1-2 144A-CA3-2 144A-DG-2 144B-CA1-1 144B-CA3-1
    ## 0610007P14Rik   4.315784   5.024361 4.3040022  4.3802738   4.541039
    ## 0610009B22Rik   2.781585   3.323114 1.9484477  2.3612882   2.597699
    ## 0610009L18Rik   1.624981   2.569098 0.7387162  0.7513143   1.751894
    ##               144C-CA1-2 144C-CA3-2 144C-DG-2 144D-CA3-2 144D-DG-2
    ## 0610007P14Rik  3.9505098   4.198110  4.143447   4.022004  4.482764
    ## 0610009B22Rik  2.2556415   3.355823  2.543533   2.359550  2.288544
    ## 0610009L18Rik  0.9921251   1.676502  1.565258   1.807928  1.570758
    ##               145A-CA1-2 145A-CA3-2 145A-DG-2 145B-CA1-1 145B-DG-1
    ## 0610007P14Rik   4.361676   3.768841 4.3495000   4.322968 3.6685090
    ## 0610009B22Rik   2.575870   3.269867 2.7493582   2.637826 2.5176543
    ## 0610009L18Rik   1.868978   1.062914 0.8012805   1.301150 0.8019297
    ##               146A-CA1-2 146A-CA3-2 146A-DG-2 146B-CA1-2 146B-CA3-2
    ## 0610007P14Rik   3.837906   4.667209  4.027906  3.5106277   3.714284
    ## 0610009B22Rik   2.739585   2.080128  2.323969  2.4018882   3.364989
    ## 0610009L18Rik   1.895555   1.759342  1.919371  0.8372935   1.119121
    ##               146B-DG-2 146C-CA1-4 146C-CA3-4 146C-DG-4 146D-CA1-3
    ## 0610007P14Rik  3.375583   4.420197   3.941227 4.8579052  4.5119131
    ## 0610009B22Rik  2.326103   2.448272   3.695130 2.9049812  2.8073485
    ## 0610009L18Rik  1.361462   2.155573   1.123693 0.9407276  0.9977117
    ##               146D-CA3-3 146D-DG-3 147-CA1-4 147-CA3-4 147-DG-4 147C-CA1-3
    ## 0610007P14Rik   4.442570  5.140646  3.217367 6.1678903 3.292688   4.098955
    ## 0610009B22Rik   2.631906  2.386915  2.196470 3.7000257 4.230663   2.308821
    ## 0610009L18Rik   1.501784  1.407197  1.263970 0.8986993 1.309172   1.013464
    ##               147C-CA3-3 147C-DG-3 147D-CA3-1 147D-DG-1 148-CA1-2
    ## 0610007P14Rik   4.511103  4.015787   3.910596  4.281375  4.253270
    ## 0610009B22Rik   2.393126  2.764598   2.705518  2.688821  2.057086
    ## 0610009L18Rik   1.450377  1.042930   1.423413  2.046641  2.150567
    ##               148-CA3-2  148-DG-2 148A-CA1-3 148A-CA3-3 148A-DG-3
    ## 0610007P14Rik 4.2261134 3.9673519   4.240433   4.111058  4.321223
    ## 0610009B22Rik 3.3295193 2.5447483   2.856001   2.732990  2.175491
    ## 0610009L18Rik 0.7566737 0.7623257   1.615507   1.572545  1.595821
    ##               148B-CA1-4 148B-CA3-4 148B-DG-4  16-116B  16-116D   16-117D
    ## 0610007P14Rik   3.822119   4.636994 4.0324351 3.923428 4.569736 4.3194252
    ## 0610009B22Rik   1.983420   2.990018 1.6532431 2.746618 2.309903 2.8673326
    ## 0610009L18Rik   1.081673   1.671509 0.8718916 1.499306 1.349279 0.8036222
    ##                16-118B  16-118D  16-119B   16-119D  16-120B   16-120D
    ## 0610007P14Rik 3.231132 3.961983 3.917615 3.6949404 3.432437 3.8641769
    ## 0610009B22Rik 2.177749 2.840042 2.544861 2.4267154 2.360040 1.8160630
    ## 0610009L18Rik 1.031291 1.867782 1.048538 0.7489247 1.150785 0.7620841
    ##                16-122B  16-122D   16-123B  16-123D  16-124D  16-125B
    ## 0610007P14Rik 4.002160 4.070174 3.4883578 3.866692 3.758364 3.730555
    ## 0610009B22Rik 3.246434 2.753782 1.6739283 2.524970 2.723266 1.948988
    ## 0610009L18Rik 1.589150 1.367223 0.8810047 1.403181 1.190391 1.049445
    ##                 16-125D  16-126B
    ## 0610007P14Rik 3.6909519 3.860763
    ## 0610009B22Rik 2.2786624 2.601611
    ## 0610009L18Rik 0.7669431 2.118743

pca plot
--------

    pcaData <- plotPCA(rld, intgroup = c( "Group", "Punch"), returnData=TRUE)
    pcaData

    ##                    PC1           PC2            group      Group Punch
    ## 142C_CA1   -17.6217815  -8.778768944 consistent : CA1 consistent   CA1
    ## 142C_DG     36.7572071  -2.106459059  consistent : DG consistent    DG
    ## 143A-CA3-1  -4.2428789  -0.595082613   conflict : CA3   conflict   CA3
    ## 143A-DG-1   36.4855509  -2.614206001    conflict : DG   conflict    DG
    ## 143B-CA1-1 -15.4635085  -7.387143733    control : CA1    control   CA1
    ## 143B-DG-1   33.8991255  -1.152194586     control : DG    control    DG
    ## 143C_CA1   -17.4766988  -9.774124786 consistent : CA1 consistent   CA1
    ## 143C_DG     31.1351747  -1.183908770  consistent : DG consistent    DG
    ## 143C-CA1-1 -17.1673933  -9.113686485 consistent : CA1 consistent   CA1
    ## 143D-CA1-3 -17.4245996  -6.064502867    control : CA1    control   CA1
    ## 143D-DG-3   34.1204441  -0.048242947     control : DG    control    DG
    ## 144A-CA1-2  -7.1425486  -5.368832267   conflict : CA1   conflict   CA1
    ## 144A-CA3-2  -1.9840174   3.751744059   conflict : CA3   conflict   CA3
    ## 144A-DG-2   37.0901881  -1.742644446    conflict : DG   conflict    DG
    ## 144B-CA1-1 -18.1426634  -6.802403037    control : CA1    control   CA1
    ## 144B-CA3-1   1.8082446   1.466670519    control : CA3    control   CA3
    ## 144C-CA1-2  -8.8119235  -6.174864238 consistent : CA1 consistent   CA1
    ## 144C-CA3-2  -2.2729942   0.441746449 consistent : CA3 consistent   CA3
    ## 144C-DG-2   38.2025940  -1.995787053  consistent : DG consistent    DG
    ## 144D-CA3-2  -2.1217724  -0.183534562    control : CA3    control   CA3
    ## 144D-DG-2   37.1857930  -3.072705100     control : DG    control    DG
    ## 145A-CA1-2 -16.7764103  -8.952248185   conflict : CA1   conflict   CA1
    ## 145A-CA3-2   1.9751354   3.496102810   conflict : CA3   conflict   CA3
    ## 145A-DG-2   34.5523840  -1.661233548    conflict : DG   conflict    DG
    ## 145B-CA1-1 -16.0955292  -6.882171508    control : CA1    control   CA1
    ## 145B-DG-1   35.9221363   0.006470756     control : DG    control    DG
    ## 146A-CA1-2 -15.4035815  -7.637142256   conflict : CA1   conflict   CA1
    ## 146A-CA3-2   0.1484925   1.681960387   conflict : CA3   conflict   CA3
    ## 146A-DG-2   34.7691092  -1.439048199    conflict : DG   conflict    DG
    ## 146B-CA1-2 -13.3613563  -5.295131113    control : CA1    control   CA1
    ## 146B-CA3-2  -4.0099302   0.928017454    control : CA3    control   CA3
    ## 146B-DG-2   13.9568512   4.632470634     control : DG    control    DG
    ## 146C-CA1-4 -15.3245921  -3.457425290 consistent : CA1 consistent   CA1
    ## 146C-CA3-4  -2.9811742  83.347088004 consistent : CA3 consistent   CA3
    ## 146C-DG-4   33.0593342   0.563385690  consistent : DG consistent    DG
    ## 146D-CA1-3 -11.3995090  -0.501401973    control : CA1    control   CA1
    ## 146D-CA3-3  -1.0281664   1.557537566    control : CA3    control   CA3
    ## 146D-DG-3   14.3993347   5.766273516     control : DG    control    DG
    ## 147-CA1-4  -11.2615152   2.068628743   homecage : CA1   homecage   CA1
    ## 147-CA3-4   -7.0457427   3.538462580   homecage : CA3   homecage   CA3
    ## 147-DG-4    16.2203328   5.820075649    homecage : DG   homecage    DG
    ## 147C-CA1-3 -15.0959865  -8.005842598 consistent : CA1 consistent   CA1
    ## 147C-CA3-3  -0.6353310   1.331362524 consistent : CA3 consistent   CA3
    ## 147C-DG-3   37.9838581  -2.347806624  consistent : DG consistent    DG
    ## 147D-CA3-1  -3.2182667   0.322962501    control : CA3    control   CA3
    ## 147D-DG-1   38.7695080  -1.884763487     control : DG    control    DG
    ## 148-CA1-2  -18.3957892  -6.332675357   homecage : CA1   homecage   CA1
    ## 148-CA3-2   -2.5466454   0.368838366   homecage : CA3   homecage   CA3
    ## 148-DG-2    35.1144881   0.612154297    homecage : DG   homecage    DG
    ## 148A-CA1-3 -11.4230932  -7.343152415   conflict : CA1   conflict   CA1
    ## 148A-CA3-3  -3.7926034  -0.292736065   conflict : CA3   conflict   CA3
    ## 148A-DG-3   35.1750697  -1.301996302    conflict : DG   conflict    DG
    ## 148B-CA1-4 -14.2366518   3.750228794    control : CA1    control   CA1
    ## 148B-CA3-4  -2.8141316  -0.184686262    control : CA3    control   CA3
    ## 148B-DG-4   16.3324050   1.580922657     control : DG    control    DG
    ## 16-116B    -18.9453556  -6.185036361    control : CA1    control   CA1
    ## 16-116D    -25.2068349 108.013164644    control : CA1    control   CA1
    ## 16-117D    -17.4188173  -4.183739050    control : CA1    control   CA1
    ## 16-118B    -19.5408802  -8.827323819    control : CA1    control   CA1
    ## 16-118D    -18.1700561  -8.144639400    control : CA1    control   CA1
    ## 16-119B    -19.9230415  -8.857167110    control : CA1    control   CA1
    ## 16-119D    -19.5394293  -8.117712123    control : CA1    control   CA1
    ## 16-120B    -19.2926068  -8.129499180    control : CA1    control   CA1
    ## 16-120D    -18.9786403  -7.154410723    control : CA1    control   CA1
    ## 16-122B    -17.6784812  -7.979936464    control : CA1    control   CA1
    ## 16-122D    -18.2406541  -8.506924677    control : CA1    control   CA1
    ## 16-123B    -15.3116106   0.036736389    control : CA1    control   CA1
    ## 16-123D    -17.7563578  -7.828255918    control : CA1    control   CA1
    ## 16-124D    -20.0920276  -5.840201214    control : CA1    control   CA1
    ## 16-125B    -16.2569691   7.251512972    control : CA1    control   CA1
    ## 16-125D    -18.0891205  -6.769995182    control : CA1    control   CA1
    ## 16-126B    -17.9030922  -8.131124068    control : CA1    control   CA1
    ##                  name
    ## 142C_CA1     142C_CA1
    ## 142C_DG       142C_DG
    ## 143A-CA3-1 143A-CA3-1
    ## 143A-DG-1   143A-DG-1
    ## 143B-CA1-1 143B-CA1-1
    ## 143B-DG-1   143B-DG-1
    ## 143C_CA1     143C_CA1
    ## 143C_DG       143C_DG
    ## 143C-CA1-1 143C-CA1-1
    ## 143D-CA1-3 143D-CA1-3
    ## 143D-DG-3   143D-DG-3
    ## 144A-CA1-2 144A-CA1-2
    ## 144A-CA3-2 144A-CA3-2
    ## 144A-DG-2   144A-DG-2
    ## 144B-CA1-1 144B-CA1-1
    ## 144B-CA3-1 144B-CA3-1
    ## 144C-CA1-2 144C-CA1-2
    ## 144C-CA3-2 144C-CA3-2
    ## 144C-DG-2   144C-DG-2
    ## 144D-CA3-2 144D-CA3-2
    ## 144D-DG-2   144D-DG-2
    ## 145A-CA1-2 145A-CA1-2
    ## 145A-CA3-2 145A-CA3-2
    ## 145A-DG-2   145A-DG-2
    ## 145B-CA1-1 145B-CA1-1
    ## 145B-DG-1   145B-DG-1
    ## 146A-CA1-2 146A-CA1-2
    ## 146A-CA3-2 146A-CA3-2
    ## 146A-DG-2   146A-DG-2
    ## 146B-CA1-2 146B-CA1-2
    ## 146B-CA3-2 146B-CA3-2
    ## 146B-DG-2   146B-DG-2
    ## 146C-CA1-4 146C-CA1-4
    ## 146C-CA3-4 146C-CA3-4
    ## 146C-DG-4   146C-DG-4
    ## 146D-CA1-3 146D-CA1-3
    ## 146D-CA3-3 146D-CA3-3
    ## 146D-DG-3   146D-DG-3
    ## 147-CA1-4   147-CA1-4
    ## 147-CA3-4   147-CA3-4
    ## 147-DG-4     147-DG-4
    ## 147C-CA1-3 147C-CA1-3
    ## 147C-CA3-3 147C-CA3-3
    ## 147C-DG-3   147C-DG-3
    ## 147D-CA3-1 147D-CA3-1
    ## 147D-DG-1   147D-DG-1
    ## 148-CA1-2   148-CA1-2
    ## 148-CA3-2   148-CA3-2
    ## 148-DG-2     148-DG-2
    ## 148A-CA1-3 148A-CA1-3
    ## 148A-CA3-3 148A-CA3-3
    ## 148A-DG-3   148A-DG-3
    ## 148B-CA1-4 148B-CA1-4
    ## 148B-CA3-4 148B-CA3-4
    ## 148B-DG-4   148B-DG-4
    ## 16-116B       16-116B
    ## 16-116D       16-116D
    ## 16-117D       16-117D
    ## 16-118B       16-118B
    ## 16-118D       16-118D
    ## 16-119B       16-119B
    ## 16-119D       16-119D
    ## 16-120B       16-120B
    ## 16-120D       16-120D
    ## 16-122B       16-122B
    ## 16-122D       16-122D
    ## 16-123B       16-123B
    ## 16-123D       16-123D
    ## 16-124D       16-124D
    ## 16-125B       16-125B
    ## 16-125D       16-125D
    ## 16-126B       16-126B

    percentVar <- round(100 * attr(pcaData, "percentVar"))

    ggplot(pcaData, aes(PC1, PC2, color=Group, shape = Punch)) + geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()

![](../results/all/pca-1.png)

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
    ## [22] labeling_0.3         splines_3.3.1        BiocParallel_1.8.1  
    ## [25] geneplotter_1.52.0   stringr_1.1.0        foreign_0.8-67      
    ## [28] RCurl_1.95-4.8       munsell_0.4.3        htmltools_0.3.5     
    ## [31] nnet_7.3-12          tibble_1.2           gridExtra_2.2.1     
    ## [34] htmlTable_1.7        Hmisc_4.0-0          XML_3.98-1.4        
    ## [37] bitops_1.0-6         grid_3.3.1           xtable_1.8-2        
    ## [40] gtable_0.2.0         DBI_0.5-1            scales_0.4.0        
    ## [43] KernSmooth_2.23-15   stringi_1.1.2        XVector_0.14.0      
    ## [46] genefilter_1.56.0    latticeExtra_0.6-28  Formula_1.2-1       
    ## [49] RColorBrewer_1.1-2   tools_3.3.1          survival_2.40-1     
    ## [52] yaml_2.1.14          AnnotationDbi_1.36.0 colorspace_1.2-7    
    ## [55] cluster_2.0.5        caTools_1.17.1       knitr_1.15.1
