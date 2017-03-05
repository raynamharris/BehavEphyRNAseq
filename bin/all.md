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

    ## Warning in if (file.exists(kallistoFiles)) kallistoData =
    ## lapply(kallistoFiles, : the condition has length > 1 and only the first
    ## element will be used

    ##     142C_CA1           142C_DG          143A-CA3-1      
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.0   Median :    0.00  
    ##  Mean   :   82.63   Mean   :   69.9   Mean   :   56.68  
    ##  3rd Qu.:   32.00   3rd Qu.:   26.0   3rd Qu.:   17.00  
    ##  Max.   :46130.00   Max.   :17459.0   Max.   :23989.00  
    ##    143A-DG-1          143B-CA1-1         143B-DG-1       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   89.92   Mean   :   29.29   Mean   :   35.51  
    ##  3rd Qu.:   35.00   3rd Qu.:    9.00   3rd Qu.:   15.00  
    ##  Max.   :76185.00   Max.   :30026.00   Max.   :21691.00  
    ##     143C_CA1          143C_DG           143C-CA1-1        143D-CA1-3      
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    0.0   Median :    0.00   Median :    0.0   Median :    0.00  
    ##  Mean   :  125.3   Mean   :   61.17   Mean   :   37.7   Mean   :   18.59  
    ##  3rd Qu.:   45.0   3rd Qu.:   25.00   3rd Qu.:   13.0   3rd Qu.:    5.00  
    ##  Max.   :44580.0   Max.   :33983.00   Max.   :21143.0   Max.   :14110.00  
    ##    143D-DG-3          144A-CA1-2         144A-CA3-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.000  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.000  
    ##  Median :    0.00   Median :    0.00   Median :    0.000  
    ##  Mean   :   17.78   Mean   :   50.77   Mean   :    7.173  
    ##  3rd Qu.:    7.00   3rd Qu.:   19.00   3rd Qu.:    2.000  
    ##  Max.   :10111.00   Max.   :44270.00   Max.   :12302.000  
    ##    144A-DG-2          144B-CA1-1         144B-CA3-1     
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    0.00   Median :    0.00   Median :    0.0  
    ##  Mean   :   54.67   Mean   :   43.53   Mean   :   17.5  
    ##  3rd Qu.:   21.00   3rd Qu.:   14.00   3rd Qu.:    5.0  
    ##  Max.   :33414.00   Max.   :35177.00   Max.   :12714.0  
    ##    144C-CA1-2         144C-CA3-2        144C-DG-2       
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.0   Median :    0.00  
    ##  Mean   :   56.18   Mean   :   21.1   Mean   :   37.88  
    ##  3rd Qu.:   20.00   3rd Qu.:    6.0   3rd Qu.:   15.00  
    ##  Max.   :35028.00   Max.   :26898.0   Max.   :15607.00  
    ##    144D-CA3-2         144D-DG-2         145A-CA1-2      
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    1.0   Median :    0.00  
    ##  Mean   :   39.57   Mean   :   79.9   Mean   :   79.72  
    ##  3rd Qu.:   13.00   3rd Qu.:   32.0   3rd Qu.:   27.00  
    ##  Max.   :46442.00   Max.   :39287.0   Max.   :73533.00  
    ##    145A-CA3-2          145A-DG-2          145B-CA1-1     
    ##  Min.   :    0.000   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.000   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    0.000   Median :    0.00   Median :    0.0  
    ##  Mean   :    5.886   Mean   :   24.45   Mean   :   34.4  
    ##  3rd Qu.:    1.000   3rd Qu.:    9.00   3rd Qu.:   10.0  
    ##  Max.   :15478.000   Max.   :11567.00   Max.   :39155.0  
    ##    145B-DG-1         146A-CA1-2         146A-CA3-2      
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.0   Median :    0.00   Median :    0.00  
    ##  Mean   :   25.7   Mean   :   29.21   Mean   :   46.94  
    ##  3rd Qu.:   10.0   3rd Qu.:   10.00   3rd Qu.:   14.00  
    ##  Max.   :15903.0   Max.   :28994.00   Max.   :89329.00  
    ##    146A-DG-2          146B-CA1-2         146B-CA3-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   20.46   Mean   :   18.11   Mean   :   36.53  
    ##  3rd Qu.:    8.00   3rd Qu.:    6.00   3rd Qu.:   11.00  
    ##  Max.   :13368.00   Max.   :10320.00   Max.   :23297.00  
    ##    146B-DG-2          146C-CA1-4         146C-CA3-4      
    ##  Min.   :   0.000   Min.   :    0.00   Min.   :   0.000  
    ##  1st Qu.:   0.000   1st Qu.:    0.00   1st Qu.:   0.000  
    ##  Median :   0.000   Median :    0.00   Median :   0.000  
    ##  Mean   :   1.977   Mean   :   23.16   Mean   :   4.391  
    ##  3rd Qu.:   0.000   3rd Qu.:    8.00   3rd Qu.:   0.000  
    ##  Max.   :7802.000   Max.   :15918.00   Max.   :6848.000  
    ##    146C-DG-4          146D-CA1-3          146D-CA3-3      146D-DG-3      
    ##  Min.   :   0.000   Min.   :    0.000   Min.   :    0   Min.   :   0.00  
    ##  1st Qu.:   0.000   1st Qu.:    0.000   1st Qu.:    0   1st Qu.:   0.00  
    ##  Median :   0.000   Median :    0.000   Median :    0   Median :   0.00  
    ##  Mean   :   8.382   Mean   :    6.666   Mean   :   51   Mean   :   1.54  
    ##  3rd Qu.:   3.000   3rd Qu.:    0.000   3rd Qu.:   16   3rd Qu.:   0.00  
    ##  Max.   :2795.000   Max.   :19366.000   Max.   :38276   Max.   :6004.00  
    ##    147-CA1-4          147-CA3-4           147-DG-4       
    ##  Min.   :   0.000   Min.   :    0.00   Min.   :   0.000  
    ##  1st Qu.:   0.000   1st Qu.:    0.00   1st Qu.:   0.000  
    ##  Median :   0.000   Median :    0.00   Median :   0.000  
    ##  Mean   :   2.709   Mean   :   11.74   Mean   :   2.372  
    ##  3rd Qu.:   0.000   3rd Qu.:    0.00   3rd Qu.:   0.000  
    ##  Max.   :8956.000   Max.   :28974.00   Max.   :7543.000  
    ##    147C-CA1-3         147C-CA3-3          147C-DG-3      
    ##  Min.   :    0.00   Min.   :     0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:     0.00   1st Qu.:    0.0  
    ##  Median :    0.00   Median :     0.00   Median :    1.0  
    ##  Mean   :   52.32   Mean   :    98.01   Mean   :   74.1  
    ##  3rd Qu.:   20.00   3rd Qu.:    37.00   3rd Qu.:   34.0  
    ##  Max.   :37687.00   Max.   :150301.00   Max.   :46988.0  
    ##    147D-CA3-1         147D-DG-1         148-CA1-2       
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.0   Median :    0.00  
    ##  Mean   :   78.77   Mean   :  199.3   Mean   :   32.38  
    ##  3rd Qu.:   27.00   3rd Qu.:   82.0   3rd Qu.:    8.00  
    ##  Max.   :95754.00   Max.   :91299.0   Max.   :24841.00  
    ##    148-CA3-2           148-DG-2          148A-CA1-3     
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    0.00   Median :    0.00   Median :    0.0  
    ##  Mean   :   39.91   Mean   :   38.01   Mean   :   89.6  
    ##  3rd Qu.:   13.00   3rd Qu.:   12.00   3rd Qu.:   31.0  
    ##  Max.   :23437.00   Max.   :23190.00   Max.   :52783.0  
    ##    148A-CA3-3         148A-DG-3          148B-CA1-4      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   45.58   Mean   :   68.45   Mean   :    5.74  
    ##  3rd Qu.:   13.00   3rd Qu.:   29.00   3rd Qu.:    0.00  
    ##  Max.   :32891.00   Max.   :31971.00   Max.   :33665.00  
    ##    148B-CA3-4         148B-DG-4          16-116B        
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.0   Median :    0.00  
    ##  Mean   :   59.38   Mean   :   13.6   Mean   :   35.47  
    ##  3rd Qu.:   19.00   3rd Qu.:    5.0   3rd Qu.:   11.00  
    ##  Max.   :37680.00   Max.   :10089.0   Max.   :23015.00  
    ##     16-116D            16-117D            16-118B        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :    4.37   Mean   :   24.49   Mean   :   49.45  
    ##  3rd Qu.:    0.00   3rd Qu.:    8.00   3rd Qu.:   16.00  
    ##  Max.   :51919.00   Max.   :18898.00   Max.   :15698.00  
    ##     16-118D            16-119B            16-119D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   37.33   Mean   :   44.62   Mean   :   44.17  
    ##  3rd Qu.:   13.00   3rd Qu.:   14.00   3rd Qu.:   14.00  
    ##  Max.   :15547.00   Max.   :20043.00   Max.   :13698.00  
    ##     16-120B            16-120D            16-122B        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   48.87   Mean   :   37.38   Mean   :   47.32  
    ##  3rd Qu.:   16.00   3rd Qu.:   11.00   3rd Qu.:   16.00  
    ##  Max.   :21743.00   Max.   :21114.00   Max.   :21913.00  
    ##     16-122D            16-123B            16-123D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   54.55   Mean   :   14.42   Mean   :   43.46  
    ##  3rd Qu.:   18.00   3rd Qu.:    3.00   3rd Qu.:   15.00  
    ##  Max.   :21640.00   Max.   :33189.00   Max.   :16867.00  
    ##     16-124D            16-125B            16-125D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   44.21   Mean   :    8.27   Mean   :   34.99  
    ##  3rd Qu.:   14.00   3rd Qu.:    2.00   3rd Qu.:   12.00  
    ##  Max.   :22218.00   Max.   :90359.00   Max.   :13073.00  
    ##     16-126B        
    ##  Min.   :    0.00  
    ##  1st Qu.:    0.00  
    ##  Median :    0.00  
    ##  Mean   :   45.99  
    ##  3rd Qu.:   16.00  
    ##  Max.   :24322.00

    ##     142C_CA1           142C_DG          143A-CA3-1      
    ##  Min.   :    0.00   Min.   :   0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:   0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :   0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :  17.03   Mean   :   17.03  
    ##  3rd Qu.:    9.00   3rd Qu.:  11.00   3rd Qu.:    7.00  
    ##  Max.   :19994.00   Max.   :9036.00   Max.   :13308.00  
    ##    143A-DG-1          143B-CA1-1         143B-DG-1       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    8.00   3rd Qu.:    6.00   3rd Qu.:    8.00  
    ##  Max.   :25544.00   Max.   :27782.00   Max.   :18912.00  
    ##     143C_CA1           143C_DG           143C-CA1-1      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    9.00   3rd Qu.:   10.00   3rd Qu.:    7.00  
    ##  Max.   :13133.00   Max.   :18685.00   Max.   :18125.00  
    ##    143D-CA1-3         143D-DG-3          144A-CA1-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    5.00   3rd Qu.:    9.00   3rd Qu.:    7.00  
    ##  Max.   :19703.00   Max.   :16880.00   Max.   :24613.00  
    ##    144A-CA3-2         144A-DG-2          144B-CA1-1      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    3.00   3rd Qu.:    9.00   3rd Qu.:    7.00  
    ##  Max.   :38260.00   Max.   :19211.00   Max.   :23010.00  
    ##    144B-CA3-1         144C-CA1-2         144C-CA3-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    5.00   3rd Qu.:    8.00   3rd Qu.:    5.00  
    ##  Max.   :20084.00   Max.   :19363.00   Max.   :30938.00  
    ##    144C-DG-2          144D-CA3-2         144D-DG-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    9.00   3rd Qu.:    6.00   3rd Qu.:    9.00  
    ##  Max.   :13913.00   Max.   :35044.00   Max.   :15967.00  
    ##    145A-CA1-2         145A-CA3-2         145A-DG-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    6.00   3rd Qu.:    1.00   3rd Qu.:    9.00  
    ##  Max.   :27499.00   Max.   :38595.00   Max.   :14783.00  
    ##    145B-CA1-1         145B-DG-1          146A-CA1-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    6.00   3rd Qu.:    9.00   3rd Qu.:    6.00  
    ##  Max.   :31778.00   Max.   :18049.00   Max.   :25334.00  
    ##    146A-CA3-2         146A-DG-2          146B-CA1-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    4.00   3rd Qu.:    9.00   3rd Qu.:    8.00  
    ##  Max.   :43024.00   Max.   :19323.00   Max.   :12770.00  
    ##    146B-CA3-2         146B-DG-2           146C-CA1-4      
    ##  Min.   :    0.00   Min.   :     0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:     0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :     0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :    17.03   Mean   :   17.03  
    ##  3rd Qu.:    6.00   3rd Qu.:     0.00   3rd Qu.:    8.00  
    ##  Max.   :17434.00   Max.   :118921.00   Max.   :20521.00  
    ##    146C-CA3-4         146C-DG-4          146D-CA1-3      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    0.00   3rd Qu.:    6.00   3rd Qu.:    0.00  
    ##  Max.   :26486.00   Max.   :10508.00   Max.   :48182.00  
    ##    146D-CA3-3         146D-DG-3           147-CA1-4       
    ##  Min.   :    0.00   Min.   :     0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:     0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :     0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :    17.03   Mean   :   17.03  
    ##  3rd Qu.:    7.00   3rd Qu.:     0.00   3rd Qu.:    0.00  
    ##  Max.   :21103.00   Max.   :110585.00   Max.   :51833.00  
    ##    147-CA3-4           147-DG-4           147C-CA1-3      
    ##  Min.   :    0.00   Min.   :     0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:     0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :     0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :    17.03   Mean   :   17.03  
    ##  3rd Qu.:    0.00   3rd Qu.:     0.00   3rd Qu.:    8.00  
    ##  Max.   :36128.00   Max.   :129255.00   Max.   :21638.00  
    ##    147C-CA3-3         147C-DG-3          147D-CA3-1      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    7.00   3rd Qu.:   10.00   3rd Qu.:    6.00  
    ##  Max.   :38426.00   Max.   :18962.00   Max.   :30897.00  
    ##    147D-DG-1          148-CA1-2          148-CA3-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    9.00   3rd Qu.:    5.00   3rd Qu.:    6.00  
    ##  Max.   :15216.00   Max.   :22177.00   Max.   :16496.00  
    ##     148-DG-2          148A-CA1-3         148A-CA3-3      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    7.00   3rd Qu.:    7.00   3rd Qu.:    5.00  
    ##  Max.   :12329.00   Max.   :18801.00   Max.   :19835.00  
    ##    148A-DG-3          148B-CA1-4         148B-CA3-4      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:   10.00   3rd Qu.:    0.00   3rd Qu.:    6.00  
    ##  Max.   :15700.00   Max.   :56713.00   Max.   :19770.00  
    ##    148B-DG-4           16-116B            16-116D         
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :     0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:     0.00  
    ##  Median :    0.00   Median :    0.00   Median :     0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :    17.03  
    ##  3rd Qu.:    7.00   3rd Qu.:    8.00   3rd Qu.:     0.00  
    ##  Max.   :21584.00   Max.   :54046.00   Max.   :215104.00  
    ##     16-117D            16-118B            16-118D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    8.00   3rd Qu.:    9.00   3rd Qu.:    8.00  
    ##  Max.   :35992.00   Max.   :22172.00   Max.   :12942.00  
    ##     16-119B            16-119D            16-120B        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    7.00   3rd Qu.:    9.00   3rd Qu.:    8.00  
    ##  Max.   :14617.00   Max.   :21162.00   Max.   :43054.00  
    ##     16-120D            16-122B            16-122D        
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :   17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:    8.00   3rd Qu.:    8.00   3rd Qu.:    8.00  
    ##  Max.   :41812.00   Max.   :15864.00   Max.   :13117.00  
    ##     16-123B             16-123D            16-124D        
    ##  Min.   :     0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:     0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :     0.00   Median :    0.00   Median :    0.00  
    ##  Mean   :    17.03   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:     3.00   3rd Qu.:    9.00   3rd Qu.:    8.00  
    ##  Max.   :169061.00   Max.   :19103.00   Max.   :29286.00  
    ##     16-125B          16-125D            16-126B        
    ##  Min.   :     0   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:     0   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :     0   Median :    0.00   Median :    0.00  
    ##  Mean   :    17   Mean   :   17.03   Mean   :   17.03  
    ##  3rd Qu.:     3   3rd Qu.:    8.00   3rd Qu.:    7.00  
    ##  Max.   :547257   Max.   :27181.00   Max.   :16647.00

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
