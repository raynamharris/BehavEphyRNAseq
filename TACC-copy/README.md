# RNAseq Workflow 

The bioinformatic workflow for the RNAseq portion of this project will be mostly done on TACC. To keep my scripts under version control, I keep a copy on my computer in this "TACC-copy directory". 

## Resources

I constantly referred to these webpages while optimzing this workflow.

- [Explain Shell - a interactive webpage with shell argument explainations](http://explainshell.com/)
- [Stampede User Guide](https://portal.tacc.utexas.edu/user-guides/stampede)
- [wiki page for launcher_creator.py](https://wikis.utexas.edu/display/bioiteam/launcher_creator.py)
- [Dhivya's RNAseq course](https://wikis.utexas.edu/display/bioiteam/Introduction+to+RNA+Seq+Course+2016)
- [Misha's Tag-seq Workflow](https://github.com/z0on/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)
- [Samtools Tutuorial](http://biobits.org/samtools_primer.html)

## RNAseq Jobs

All of the sample information is contained in `../data/sample_info/punches.csv` and `../data/sample_info/punches.csv`. 

These samples have been process in multiple sequencing jobs. 

| RNAseq Job | Data | Description
| :--- | :---: | :--- |
JA16268 | May 24, 2016 | paired end
JA16443 | July 26, 2016 | Tag-seq
JA16444 | October 5, 2016 | paired end

For storing and working with my data, I created a directory on scratch `$SCRATCH/BehavEphyRNAseq`. Inside, there is a subdirectory for each job and then a subdirectory for the raw data.

## Workflow JA16444

I use the command line to create a lot of the commands files and batch scripts, so this readme documents the workflow I use to create these files.  My hope is that this style will help me easily reproduce the workflow (either for reanalysing this project or for applying the workflow to new datasets).

I decided to keep my data and my scripts and my commands files in the same directory, gasp! I tried to keep them separate for a while, but having the commands files in the same directory with the launcher script and the data just makes life easier. I do move the fastqc output files to their own directory after processing because that makes me happy.

So, inside `rawdata` I have the data from job JA16444 as well as the following scripts and commands files:

| Scripts | Description
| :--- | :--- | 
00_gsaf_download | Bash script to download data from amazon cloud. See    https://wikis.utexas.edu/display/GSAF/How+to+download+your+data
01_fastqc | Quality control of raw data 
02_kallistoquant | gene quantification using kallisto with the Mus transcriptome
03_kallistoquantcandidategenes | gene quantification using kallisto with only a subset of candidate genes from the Mus transcriptome

##  00_gsaf_download 

Inside the project directory on scratch `$SCRATCH/BehavEphyRNAseq`, create subdirectory for the job and then a subdirectory for the raw data. Navidate to the directory

~~~ {.bash}
mkdir -p $SCRATCH/BehavEphyRNAseq/JA16444/00_rawdata
cd $SCRATCH/BehavEphyRNAseq/JA16444/00_rawdata
~~~ 

To download, I'll copy a  bash script named `00_gsaf_download.sh` from a one of my folders. This bash script was provide by the GSAF (https://wikis.utexas.edu/display/GSAF/How+to+download+your+data). The script was copied wihtout modification, but I added a few lines at the end to clean up the directory (by making a wgetlog directory for the download stats files) and to make all the sequences files read only (with the `chmod` command).

~~~ {.bash}
# move into directory and copy script from other location
cp $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata/00_gsaf_download.sh .
~~~ 

Then I executed the following command. Note: the key is only active for 72 hours. 

~~~ {.bash}
00_gsaf_download.sh "http://gsaf.s3.amazonaws.com/JA16444.SA16163.html?AWSAccessKeyId=AKIAIVYXWYWNPBNEDIAQ&Expires=1476550876&Signature=Ia7MlL%2FVTW9x0TLTBL0BIGEY%2FT8%3D"
~~~ 


## 01_fastqc 
The first thing to do is check the quality of the reads. 

I used the following for loop, to create a commands file for fastq

~~~ {.bash}
for file in *.fastq.gz
do
     echo $file
     echo "fastqc $file" >> 01_fastqc.cmds
done
~~~

Check to see that the commands file looks like it should

~~~ {.bash}
cat 01_fastqc.cmds
~~~

Then, I ran these two comamands to create and launch a fastqc job

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 01_fastqc.cmds -n fastqc -l 01_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 01_fastqc.slurm
~~~

Then, I moved all the output files to a separate folder, with the date pre-appended (for JA1443) or with command order "01_" pre-appended (for JA1444).

~~~ {.bash}
mkdir ../01_fastqc
mv *.html ../01_fastqc
mv *.zip ../01_fastqc
mv fastqc.* ../01_fastqc
~~~

One must save the data locally in order to view the html files. 

In a new terminal window:

~~~ {.bash}
cd /Users/raynamharris/Github/BehavEphyRNAseq/TACC-copy
mkdir JA16444
cd JA16444
scp rmharris@stampede.tacc.utexas.edu:/scratch/02189/rmharris/BehavEphyRNAseq/JA16444/01_fastqc/*html .
~~~

####################
## Build a kallisto_index ## this should go somewhere else because its only done once and not for every project!!!

Download mouse transcriptome from https://www.gencodegenes.org/mouse_releases/current.html

~~~ {.bash}
mkdir $SCRATCH/BehavEphyRNAseq/refs
cd $SCRATCH/BehavEphyRNAseq/refs
curl -O 
curl -O ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.pc_transcripts.fa.gz
~~~

Then, create the commands file.

~~~ {.bash}
cd $SCRATCH/BehavEphyRNAseq/JA16444/00_rawdata
echo "kallisto index -i gencode.vM11.pc_transcripts_kallisto.idx $SCRATCH/BehavEphyRNAseq/refs/gencode.vM11.pc_transcripts.fa.gz" > 02_kallisto_index.cmds
cat 02_kallisto_index.cmds
~~~

Then create the launcher script. 

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 02_kallisto_index.cmds -n kallistoindex -l 02_kallisto_index.slurm -A NeuroEthoEvoDevo -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 02_kallisto_index.slurm
~~~

Moved the outputs to the refs folder to keep it all the same


##################

### 02_kallistoquant

Quanitify gene expression with kallisto quant.
See https://pachterlab.github.io/kallisto/manual for details.


Create the commands file. 

~~~ {.bash}
mkdir ../02_kallistoquant_largemem
for R1 in *R1_001.fastq.gz; do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    samp=$(basename $R1 _R1_001.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -i $SCRATCH/BehavEphyRNAseq/refs/gencode.vM11.pc_transcripts_kallisto.idx  -o ../02_kallistoquant_largemem/${samp} $R1 $R2" >> 02_kallistoquant.cmds
done
~~~

Create the launcher script and run. **Note**: I actually had to many files to run this on the largemem and stay within the max cores limits, so this was submitted as two jobs. 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 02_kallistoquant.cmds -n kallistoquant -l 02_kallistoquant.slurm -A NeuroEthoEvoDevo -q largemem -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 02_kallistoquant.slurm
~~~

Now, save the data locally

In a new terminal window:

~~~ {.bash}
cd /Users/raynamharris/Github/BehavEphyRNAseq/TACC-copy/JA16444/
scp -r rmharris@stampede.tacc.utexas.edu:/scratch/02189/rmharris/BehavEphyRNAseq/JA16444/02_kallistoquant_largemem .
~~~

Now, remove the uninformative bits of the "sample name" so they match up with the actual sample name. 

~~~ {.bash}
for file in *
do
    sample=${file//_S*/}
    echo $file $sample
    mv $file $sample
done
~~~

Then, replace the `_` with `-`

~~~ {.bash}
for file in *
do
    sample=${file//_/-}
    echo $file $sample
    mv $file $sample
done
~~~

### 03_kallistoquantcandidategenes

Quanitify gene expression with kallisto quant USING the candidate genes only file.
See https://pachterlab.github.io/kallisto/manual for details.

Create the commands file. 

~~~ {.bash}
mkdir ../03_kallistoquantcandidategenes
for R1 in *R1_001.fastq.gz; do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    samp=$(basename $R1 _R1_001.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -i $SCRATCH/BehavEphyRNAseq/refs/gencode.vM7.transcripts_candidategenes_kallisto.idx  -o ../03_kallistoquantcandidategenes/${samp} $R1 $R2" >> 03_kallistoquantcandidategenes.cmds
done
~~~

Create the launcher script and run. 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -q normal -j 03_kallistoquantcandidategenes.cmds -n kallistoquantcandidategenes -l 03_kallistoquantcandidategenes.slurm -A NeuroEthoEvoDevo -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 03_kallistoquantcandidategenes.slurm
~~~

Now, save the data locally

In a new terminal window:

~~~ {.bash}
cd /Users/raynamharris/Github/BehavEphyRNAseq/TACC-copy/JA16444/
scp -r rmharris@stampede.tacc.utexas.edu:/scratch/02189/rmharris/BehavEphyRNAseq/JA16444/03_kallistoquantcandidategenes .
~~~

Now, remove the uninformative bits of the "sample name" so they match up with the actual sample name. 

~~~ {.bash}
for file in *
do
    sample=${file//_S*/}
    echo $file $sample
    mv $file $sample
done
~~~

Then, replace the `_` with `-`

~~~ {.bash}
for file in *
do
    sample=${file//_/-}
    echo $file $sample
    mv $file $sample
done
~~~


### 04 Filter and Trimming Reads (filtrimmedreads)

Okay, so kallisto is awesome because you can do it on unprocessed data. But, let's clean up our reads some. 

First, lets remove adapters with cutadapt.


~~~ {.bash}
cd $SCRATCH/BehavEphyRNAseq/JA16444/00_rawdata
mkdir ../04_filtrimmedreads
for R1 in *R1_001.fastq.gz
do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    R1filtrim=$(basename $R1 fastq.gz)filtrim.fastq.gz
    R2filtrim=$(basename $R2 fastq.gz)filtrim.fastq.gz
    echo $R1 $R2 $R1filtrim $R2filtrim
    echo "cutadapt -q 15,10 -a GATCGGAAGAGCACACGTCTGAACTCCA -A ATCGTCGGACTGTAGAACTCTGAACGTG -m 22 -o $R1filtrim -p $R2filtrim $R1 $R2" >> 04_filtrimmedreads.cmds
done
~~~

~~~ {.bash}
launcher_creator.py -t 4:00:00 -j 04_filtrimmedreads.cmds -n filtrimreads -l 04_filtrimmedreads.slurm -A NeuroEthoEvoDevo -e rayna.harris@utexas.edu 
sbatch 04_trimreads.slurm
~~~

### 05_kallistoquant_largemem

Quanitify gene expression with kallisto quant.
See https://pachterlab.github.io/kallisto/manual for details.

Note: This program needs to be run on the largemem node, which has some compute limitations. I have too many samples to process them all as one job, so I need to split it up in two.  I'll use the lane identifies (L002 and L003) to process the samples in two batches. 

~~~ {.bash}
## first make a directory for the collective output
mkdir ../05_kallistoquant_largemem
~~~

Create the L002 commands file and launcher. 

~~~ {.bash}
for R1 in *L002_R1_001.filtrim.fastq.gz
do
    R2=$(basename $R1 L002_R1_001.filtrim.fastq.gz)L002_R2_001.filtrim.fastq.gz
    samp=$(basename $R1 _L002_R1_001.filtrim.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -i $SCRATCH/BehavEphyRNAseq/refs/gencode.vM11.pc_transcripts_kallisto.idx  -o ../05_kallistoquant_largemem/${samp} $R1 $R2" >> 05_kallistoquant_L002.cmds
done
~~~

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 05_kallistoquant_L002.cmds -n kallistoquant -l 05_kallistoquant_L002.slurm -A NeuroEthoEvoDevo -q largemem -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 05_kallistoquant_L002.slurm
~~~

Create the L003 commands file and launcher. 

~~~ {.bash}
for R1 in *L003_R1_001.filtrim.fastq.gz
do
    R2=$(basename $R1 L003_R1_001.filtrim.fastq.gz)L003_R2_001.filtrim.fastq.gz
    samp=$(basename $R1 _L003_R1_001.filtrim.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -i $SCRATCH/BehavEphyRNAseq/refs/gencode.vM11.pc_transcripts_kallisto.idx  -o ../05_kallistoquant_largemem/${samp} $R1 $R2" >> 05_kallistoquant_L003.cmds
done
~~~

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 05_kallistoquant_L003.cmds -n kallistoquant -l 05_kallistoquant_L003.slurm -A NeuroEthoEvoDevo -q largemem -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 05_kallistoquant_L003.slurm
~~~

Now, save the data locally

In a new terminal window:

~~~ {.bash}
cd /Users/raynamharris/Github/BehavEphyRNAseq/TACC-copy/JA16444/
scp -r rmharris@stampede.tacc.utexas.edu:/scratch/02189/rmharris/BehavEphyRNAseq/JA16444/05_kallistoquant_largemem .
~~~

Now, remove the uninformative bits of the "sample name" so they match up with the actual sample name. 

~~~ {.bash}
for file in *
do
    sample=${file//_S*/}
    echo $file $sample
    mv $file $sample
done
~~~

Then, replace the `_` with `-`

~~~ {.bash}
for file in *
do
    sample=${file//_/-}
    echo $file $sample
    mv $file $sample
done
~~~


### 06_kallistoquantcandidategenes

Quanitify gene expression with kallisto quant using the "candidate genes transcriptome".

Create the directory for storing the output. 

~~~ {.bash}
mkdir ../06_kallistoquantcandidategenes
~~~

Create the commands file. 

~~~ {.bash}
for R1 in *R1_001.filtrim.fastq.gz
do
    R2=$(basename $R1 R1_001.filtrim.fastq.gz)R2_001.filtrim.fastq.gz
    samp=$(basename $R1 _R1_001.filtrim.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -i $SCRATCH/BehavEphyRNAseq/refs/gencode.vM7.transcripts_candidategenes_kallisto.idx  -o ../06_kallistoquantcandidategenes/${samp} $R1 $R2" >> 06_kallistoquantcandidategenes.cmds
done
~~~

Create the launcher script and run. 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -q normal -j 06_kallistoquantcandidategenes.cmds -n kallistoquantcandidategenes -l 06_kallistoquantcandidategenes.slurm -A NeuroEthoEvoDevo -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 06_kallistoquantcandidategenes.slurm
~~~

Now, save the data locally

In a new terminal window:

~~~ {.bash}
cd /Users/raynamharris/Github/BehavEphyRNAseq/TACC-copy/JA16444/
scp -r rmharris@stampede.tacc.utexas.edu:/scratch/02189/rmharris/BehavEphyRNAseq/JA16444/06_kallistoquantcandidategenes .
cd 06_kallistoquantcandidategenes
~~~

Now, remove the uninformative bits of the "sample name" so they match up with the actual sample name. 

~~~ {.bash}
for file in *
do
    sample=${file//_S*/}
    echo $file $sample
    mv $file $sample
done
~~~

To match the folder name with the sample name, we must replace the `_` with `-`.

~~~ {.bash}
for file in *
do
    sample=${file//_/-}
    echo $file $sample
    mv $file $sample
done
~~~










## Workflow JA16443

### JA16443
So, inside `2016-07-26-rawdata` I have the data from job JA16443 as well as the following scripts and commands files:

| Scripts | Description
| :--- | :--- | 
00_gsaf_download | Bash script to download data from amazon cloud. See    https://wikis.utexas.edu/display/GSAF/How+to+download+your+data
01_fastqc | Quality control of raw data 
02_gunzip | Gunzip, keep original files
03_clean | 3 step trimming and filtering using Fastx toolkit, adapted from Misha https://github.com/raynamharris/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt
04_fastqc | Quality control of clean data

##  00_gsaf_download 

For storing and working with my data, I created a directory on scratch `$SCRATCH/BehavEphyRNAseq`. Inside a created inside for each job and a subdirectory for the raw data.

~~~ {.bash}
mkdir -p $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata
cd $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata
~~~ 
    
To download, I'll copy a  bash script named `00_gsaf_download.sh` from a one of my folders. This bash script was provide by the GSAF (https://wikis.utexas.edu/display/GSAF/How+to+download+your+data). The script was copied wihtout modification, but I added a few lines at the end to clean up the directory (by making a wgetlog directory for the download stats files) and to make all the sequences files read only (with the `chmod` command).

~~~ {.bash}
# move into directory and copy script from other location
cp $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata/00_gsaf_download.sh .
~~~ 

Then I executed the following command. Note: the key is only active for 72 hours. 

~~~ {.bash}
00_gsaf_download.sh "http://gsaf.s3.amazonaws.com/JA16444.SA16163.html?AWSAccessKeyId=AKIAIVYXWYWNPBNEDIAQ&Expires=1476550876&Signature=Ia7MlL%2FVTW9x0TLTBL0BIGEY%2FT8%3D"
~~~ 


### 01_fastqc 
The first thing to do is check the quality of the reads. 

I used the following for loop, to create a commands file for fastq

~~~ {.bash}
for file in *.fastq.gz
do
     echo $file
     echo "fastqc $file" >> 01_fastqc.cmds
done
~~~

Check to see that the commands file looks like it should

~~~ {.bash}
cat 01_fastqc.cmds
~~~

Then, I ran these two comamands to create and launch a fastqc job

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 01_fastqc.cmds -n fastqc -l 01_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 01_fastqc.slurm
~~~

Then, I moved all the output files to a separate folder, with the date pre-appended (for JA1443) or with command order "01_" pre-appended (for JA1444).

~~~ {.bash}
### for JA1443
mkdir ../<date>-fastqc
mv *.html ../<date>-fastqc
mv *.zip ../<date>-fastqc
mv fastqc.* ../<date>-fastqc
    
    
## 02_ gunzip

I'll use this for loop to create a list of gunzip commands that will have a new output with the `.fastq` extension AND it will keep the orginal `.fastq.gz` file.

~~~ {.bash}
for file in *.fastq.gz
do
	newfile=${file//.fastq.gz/.fastq}
	echo $file $newfile
	echo "gunzip -c $file > $newfile" >> 02_gunzip.cmds 
done
~~~

Check to see that the commands file looks like it should

~~~ {.bash}
cat 02_gunzip.cmds
~~~

Then, I used this launcher command to unzip the files. I removed the job stderr and stdout files after the job was complete with `rm gunzip.*`. 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 02_gunzip.cmds -n gunzip -l 02_gunzip.slurm -A NeuroEthoEvoDevo
sbatch 02_gunzip.slurm
~~~

Then, I like to rename these to `.fq` files to help keep things visually distinct

~~~ {.bash}
for file in *.fastq
do
	newfile=${file//.fastq/.fq}
	echo $file $newfile
	mv $file $newfile
done
~~~

Now, we can look at the raw data with head commands.

~~~ {.bash}
for file in *.fq
do
	echo $file
	head -32 $file | grep -E '^[NACGT]+$'
done
~~~

There are a lot of strings of As and Ts and some low quality stuff, so let's clean up the reads.

### 03_clean

Here, I'm creating a 3 step cleaning workshop using the fastx_clipper and fastq_quality_filter functions from the fastx toolkit. 

~~~ {.bash}
for file in *.fq
do
	newfile=${file//.fq/.trim.fq}
	echo $file $newfile
	echo "fastx_clipper -v -a AAAAAAAA -l 20 -Q33 -i $file | fastx_clipper -v  -a TTTTTTTT -l 20 -Q33 | fastq_quality_filter -v  -Q33 -q 20 -p 90 -o $newfile"  >> 03_clean.cmds
done
~~~ 

Then, I created and launched a job with this command.

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 03_clean.cmds -n clean -l 03_clean.slurm -A NeuroEthoEvoDevo -m 'module load fastx_toolkit/0.0.13.2'
sbatch 03_clean.slurm
~~~ 

Now, I can look at a few lines of each to see if indeed things were trimmed (and compare it to the old sequences. 

~~~ {.bash}
for file in *.fq
do
	echo $file
	head -32 $file | grep -E '^[NACGT]+$'
done
~~~ 


## 02_cutadapt  ######crap! don't like this cause its not paired end trimming
A colleque recommmends cutadapt for trimming adapters of paired end reads.

~~~ {.bash}
for R1 in *R1_001.fastq.gz 
do
     	trimmedR1="$(basename $R1 .fastq.gz).trimmed.fastq"
     	R2="$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz"
     	trimmedR2="$(basename $R2 .fastq.gz).trimmed.fastq"
     	echo $R1 $R2 $trimmedR1 $trimmedR2
     	echo "cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC -A GTCGGGTGGGTCACCGACCTTTTGCACGCCGCGAGAGGAAGTATTGGCGT -m 22 -o  $trimmedR1 -p $trimmedR2 $R1 $R2" >> 02_cutadapt.cmds
done
cat 02_cutadapt.cmds
~~~

Then, I used this launcher command to process the files. 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 02_cutadapt.cmds -n cutadapt -l 02_gunzip.slurm -A NeuroEthoEvoDevo
sbatch 02_gunzip.slurm
~~~

Clean up the folder

~~~ {.bash}
mkdir ../02_cutadapt
mv *.trim.fq ../02_cutadapt
mv cutadapt.* ../02_cutadapt
~~~

### 03_fastqc JA16444

Now, I check to see how well the cleaning worked.

Create the commands file: 

~~~ {.bash}
for file in *.trim.fq
do
     echo $file
     echo "fastqc $file" >> 03_fastqc.cmds
done
~~~ 

Create and launch the slurm file.

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 03_fastqc.cmds -n fastqc -l 03_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 03_fastqc.slurm
~~~

Then, I moved all the output files to the folder with all the other fastqc files.

~~~ {.bash}
mv *.html ../03_fastqc
mv *.zip ../03_fastqc
mv fastqc.* ../03_fastqc
cd ../03_fastqc
~~~    
    
       
### 04_fastqc  

Now, I check to see how well the cleaning worked.

Create the commands file: 

~~~ {.bash}
for file in *.trim.fq
do
     echo $file
     echo "fastqc $file" >> 04_fastqc.cmds
done
~~~ 

Create and launch the slurm file.

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 04_fastqc.cmds -n fastqc -l 04_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 04_fastqc.slurm
~~~

Then, I moved all the output files to the folder with all the other fastqc files.

~~~ {.bash}
mv *.html ../2016-08-22-fastqc
mv *.zip ../2016-08-22-fastqc
mv fastqc.* ../2016-08-22-fastqc
cd ../2016-08-22-fastqc
~~~

Used filezilla to transfer the files to my local computer and then check them online.


### Reference Genomes and Indexes

Get the mm9 reference genome. Note: I did this from the root directory, but I probably should have submitted a job!

See http://support.illumina.com/sequencing/sequencing_software/igenome.html for lists of available genomes.

~~~ {.bash}
cd $SCRATCH/BehavEphyRNAseq/
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz
tar -zxvf Mus_musculus_UCSC_mm9.tar.gz
~~~

To ways to find the names of chromosomes in genome file (.fa) and in the annotation file (.gtf)

~~~ {.bash}
grep '^>' $SCRATCH/BehavEphyRNAseq/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
cut -f 1 $SCRATCH/BehavEphyRNAseq/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf | sort | uniq -c
~~~ 


### 06_BWA

https://wikis.utexas.edu/pages/viewpage.action?pageId=67797479#Mappingtutorial(bowtie2,bwa)(GVA14)-MappingwithBWA

Start with the single end Tag-seq data. Make an output directory for the output. Create a soft link to the indexes for easier commands.

~~~ {.bash}
cd $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata
mkdir ../2016-08-24-bwa
ln -s $SCRATCH/BehavEphyRNAseq/Mus_musculus/UCSC/mm9/Sequence/BWAIndex
~~~ 

Command for  Illumina single-end reads shorter than ~70bp:
bwa aln ref.fa reads.fq > reads.sai; bwa samse ref.fa reads.sai reads.fq > aln-se.sam

~~~ {.bash}
for file in *trim.fq
do
	saifile=${file//.fq/.sai}
	samfile=${file//.fq/.sam}
	echo $file $saifile $samfile
	echo "bwa aln BWAIndex/genome.fa $file > $saifile; bwa samse BWAIndex/genome.fa $saifile $file > $samfile"  >> 06_bwa.cmds
done
cat 06_bwa.cmds
~~~ 

Create and launch job

~~~ {.bash}
launcher_creator.py -t 01:30:00 -j 06_bwa.cmds -n bwa -l 06_bwa.slurm -A NeuroEthoEvoDevo -m 'module load bwa/0.7.7'
sbatch 06_bwa.slurm
~~~ 


Did it work?

~~~ {.bash}
for file in *.sam
do
echo $file
samtools flagstat $file
done
~~~ 

Here is a slimmed version of the output. Note: I ran it twice and got different results for a few samples. Yikes!!!

~~~
142C_CA1_S_S19_L003_R1_001.trim.sam
6815744 + 0 in total (QC-passed reads + QC-failed reads)
3092451 + 0 mapped (45.37%:-nan%)

142C_CA1_S_S19_L003_R2_001.trim.sam
7864320 + 0 in total (QC-passed reads + QC-failed reads)
3783480 + 0 mapped (48.11%:-nan%)

0 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 mapped (-nan%:-nan%)

142C_DG_S_S21_L003_R1_001.trim.sam
13369344 + 0 in total (QC-passed reads + QC-failed reads)
3651980 + 0 mapped (27.32%:-nan%)

142C_DG_S_S21_L003_R2_001.trim.sam
10892560 + 0 in total (QC-passed reads + QC-failed reads)
3510025 + 0 mapped (32.22%:-nan%)

143C_CA1_S_S20_L003_R1_001.trim.sam
5767168 + 0 in total (QC-passed reads + QC-failed reads)
3609402 + 0 mapped (62.59%:-nan%)

143C_CA1_S_S20_L003_R2_001.trim.sam
0 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 mapped (-nan%:-nan%)

18690052 + 0 in total (QC-passed reads + QC-failed reads)
3265799 + 0 mapped (17.47%:-nan%)

143C_DG_S_S22_L003_R1_001.trim.sam
18690052 + 0 in total (QC-passed reads + QC-failed reads)
3265799 + 0 mapped (17.47%:-nan%)

143C_DG_S_S22_L003_R2_001.trim.sam
0 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 mapped (-nan%:-nan%)

143C_DG_S_S22_L003_R2_001.trim.sam
14276719 + 0 in total (QC-passed reads + QC-failed reads)
2860619 + 0 mapped (20.04%:-nan%)


ALL Tag-seq data

0 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 mapped (-nan%:-nan%)
~~~



No, not for the tag-seq data.


### 07_Samtools

http://biobits.org/samtools_primer.html
http://www.htslib.org/doc/samtools-1.1.html

The commands
- samtools view -b -S -o alignments/sim_reads_aligned.bam alignments/sim_reads_aligned.sam
- samtools sort -T /tmp/aln.sorted -o aln.sorted.bam aln.bam
- samtools index alignments/sim_reads_aligned.sorted.bam



#### 07a View

~~~ {.bash}
for file in *.sam
do
	bamfile=${file//.sam/.bam}
	echo $file $bamfile
	echo "samtools view -b -S $file > $bamfile" >> 07a_samtools.cmds
done
cat 07a_samtools.cmds
~~~ 

~~~ {.bash}
launcher_creator.py -t 01:00:00 -j 07a_samtools.cmds -n samtoolsa -l 07a_samtools.slurm -A NeuroEthoEvoDevo -m 'module load samtools/1.3'
sbatch 07a_samtools.slurm
~~~

#### 07b sort

Usage: samtools sort -T /tmp/aln.sorted -o aln.sorted.bam aln.bam

~~~ {.bash}
mkdir tmp
for file in *.bam
do
	sortedtmp=${file//.bam/.temp.sorted}
	sortedbam=${file//.bam/_sorted.bam}
	echo $file $sortedfile $sortedtmp
	echo "samtools sort -T temp/$sortedtmp -o $sortedbam $file" >> 07b_samtools.cmds
done
cat 07b_samtools.cmds
~~~ 

~~~ {.bash}
launcher_creator.py -t 01:00:00 -j 07b_samtools.cmds -n samtoolsb -l 07b_samtools.slurm -A NeuroEthoEvoDevo -m 'module load samtools/1.3'
sbatch 07b_samtools.slurm
~~~

#### 07c index

- samtools index alignments/sim_reads_aligned.sorted.bam


~~~ {.bash}
for file in *.bam
do
	sortedfile=${file//.bam/_sorted.bam}
	echo $file $sortedfile
	echo "samtools sort $file $sortedfile" >> 07b_samtools.cmds
done
cat 07b_samtools.cmds
~~~ 

~~~ {.bash}
launcher_creator.py -t 01:00:00 -j 07b_samtools.cmds -n samtoolsb -l 07b_samtools.slurm -A NeuroEthoEvoDevo -m 'module load samtools/1.3'
sbatch 07b_samtools.slurm
~~~


#### not working samptools pipe

But, I will do all these together with a pipe

samtools view -b -S $file | samtools sort | samtools index - $bamfile


~~~ {.bash}
for file in *.sam
do
	bamfile=${file//.sam/.bam}
	echo $file $bamfile
	echo "samtools view -b -S -o $file | samtools sort | samtools index - $bamfile"  >> 07_samtools.cmds
done
cat 07_samtools.cmds
~~~ 

~~~ {.bash}
launcher_creator.py -t 01:00:00 -j 07_samtools.cmds -n samtools -l 07_samtools.slurm -A NeuroEthoEvoDevo -m 'module load samtools/1.3'
sbatch 07_samtools.slurm
~~~


Did it work?

~~~ {.bash}
module load samtools
for file in *.bam 
do
echo $file
samtools flagstat $file
done


### 05_bowtie2 Summary Statistics

~~~ {.bash}
module load samtools
samtools flagstat JA16268.bam
~~~

~~~
10027976 + 0 in total (QC-passed reads + QC-failed reads)
5114817 + 0 mapped (51.01%:-nan%)
~~~

Repeat for JA16443

~~~ {.bash}
module load samtools
samtools flagstat JA16443.bam
~~~

~~~
14866 + 0 in total (QC-passed reads + QC-failed reads)
5877 + 0 mapped (39.53%:-nan%)
~~~

**Interpretation:** The tag-seq data produced ~5K mapped reads where the regular RNA-seq data produced ~5M mapped reads.



