# RNAseq Workflow

The bioinformatic workflow for the RNAseq portion of this project will be mostly done on TACC. To keep my scripts under version control, I keep a copy on my computer in this "TACC-copy directory". 

## RNAseq Jobs

All of the sample information is contained in `../data/sample_info/punches.csv` and `../data/sample_info/punches.csv`. 

These samples have been process in multiple sequencing jobs. 

| RNAseq Job | Data | Description
| :--- | :---: | :--- |
JA16268 | May 24, 2016 | paired end
JA16443 | July 26, 2016 | Tag-seq
JA16444 | TBA | Tag-seq

## Workflow

I use the command line to create a lot of the commands files and batch scripts, so this readme documents the workflow I use to create these files.  My hope is that this style will help me easily reproduce the workflow (either for reanalysing this project or for applying the workflow to new datasets).

I decided to keep my data and my scripts and my commands files in the same directory, gasp! I tried to keep them separate for a while, but having the commands files in the same directory with the launcher script and the data just makes life easier. I do move the fastqc output files to their own directory after processing because that makes me happy.

So, inside `2016-07-26-rawdata` I have the data from job JA16443 as well as the following scripts and commands files:

| Scripts | Description
| :--- | :--- | 
00_gsaf_download | Bash script to download data from amazon cloud. See    https://wikis.utexas.edu/display/GSAF/How+to+download+your+data
01_fastqc | Quality control of raw data 
02_gunzip | Gunzip, keep original files
03_clean | 3 step trimming and filtering using Fastx toolkit, adapted from Misha https://github.com/raynamharris/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt
04_fastqc | Quality control of clean data


## Resources

I constantly referred to these webpages while optimzing this workflow.

- [Explain Shell - a interactive webpage with shell argument explainations](http://explainshell.com/)
- [Stampede User Guide](https://portal.tacc.utexas.edu/user-guides/stampede)
- [wiki page for launcher_creator.py](https://wikis.utexas.edu/display/bioiteam/launcher_creator.py)
- [Dhivya's RNAseq course](https://wikis.utexas.edu/display/bioiteam/Introduction+to+RNA+Seq+Course+2016)
- [Misha's Tag-seq Workflow](https://github.com/z0on/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)
- [Samtools Tutuorial](http://biobits.org/samtools_primer.html)

###  00_gsaf_download

For storing and working with my data, I created a directory on scratch `$SCRATCH/BehavEphyRNAseq`. Inside a created inside for each job and a subdirectory for the raw data.

~~~ {.bash}
mkdir $SCRATCH/BehavEphyRNAseq
mkdir -p $SCRATCH/BehavEphyRNAseq/JA16268/2016-05-24-rawdata
mkdir -p $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata
~~~ 

To download the raw data, I navidated to the appropriate directory use nano to create a bash script named `00_gsaf_download.sh` with the bash script written by the GSAF stored here: https://wikis.utexas.edu/display/GSAF/How+to+download+your+data

~~~ {.bash}
cd $SCRATCH/BehavEphyRNAseq/JA16443/2016-07-26-rawdata
nano 00_gsaf_download.sh
~~~ 

Then I executed the following command. Note: the key is only active for 72 hours. 

~~~ {.bash}
00_gsaf_download.sh "http://gsaf.s3.amazonaws.com/JA16443.SA16131.html?AWSAccessKeyId=AKIAIVYXWYWNPBNEDIAQ&Expires=1469736784&Signature=EBDWx1Ke55tIRekbuN0WPrt4d6s%3D"
~~~ 

Then, I cleaned up my directory a bit with the following code:

~~~ {.bash}
mkdir wgetlog
mv *.wget.log wgetlog/
mv files.html wgetlog/
~~~

Finally, I made the files *read-only* with this command.

~~~ {.bash}
chmod u-w *.fastq.gz 
~~~ 

### 01_fastqc 

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

Then, I moved all the output files to a separate folder, with the date pre-appended.

~~~ {.bash}
mkdir ../<date>-fastqc
mv *.html ../<date>-fastqc
mv *.zip ../<date>-fastqc
mv fastqc.* ../<date>-fastqc
~~~

### 02_ gunzip

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



