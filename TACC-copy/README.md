# TaqSeq Workflow

Small RNAseq project designed to see if TaqSeq will for for hippocampal punches

## $SCRATCH/JA16443/

I decided to keep my data and my scripts and my commands files in the same directory, gasp! I tried to keep them separate for a while, but having the commands files in the same directory with the launcher script and the data just makes life easier. I do move the fastqc output files to their own directory after processing because that makes me happy.

So, inside `2016-07-26-rawdata` I have the data from job JA16443 as well as the following scripts and commands files:

| Scripts | Description | Notes
| :--- | :--- | :--- | 
00_gsaf_download | Bash script to download data from amazon cloud | See    https://wikis.utexas.edu/display/GSAF/How+to+download+your+data
01_fastqc | Quality control of raw data | Output moved to ../2016-08-22-fastqc
02_gunzip | Gunzip, keep original | deleted .e and .o files
03_clean | 3 step trimming and filtering | adapted from Misha's 2 step method https://github.com/raynamharris/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt
04_fastqc | Quality control of clean data | Output moved to ../2016-08-22-fastqc

## Workflow

More of less, this is what had to be executed at the command line. 

~~~ {.bash}
# start by creating a space on scratch where I'll download and process the raw data
mkdir $SCRATCH/JA16443/2016-07-26-rawdata
cd $SCRATCH/JA16443/2016-07-26-rawdata
~~~  

###  00_gsaf_download

Use nano to create a bash script named `00_gsaf_download.sh` with the bash script written by the GSAF stored here: https://wikis.utexas.edu/display/GSAF/How+to+download+your+data

~~~ {.bash}
nano 00_gsaf_download.sh
~~~ 

Then I executed the following command, replacing <"linktokey"> with the real webaddress inside quotes

~~~ {.bash}
00_gsaf_download.sh <"linktokey">
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
     cat >> 01_fastqc.cmds <<EOF
fastqc $file
EOF
done
~~~

Then, I ran these two comamands to create and launch a fastqc job

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 01_fastqc.cmds -n fastqc -l 01_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 01_fastqc.slurm
~~~

Then, I moved all the output files to a separate folder

~~~ {.bash}
mkdir ../2016-08-22-fastqc
mv *.html ../2016-08-22-fastqc
mv *.zip ../2016-08-22-fastqc
mv fastqc.* ../2016-08-22-fastqc
cd ../2016-08-22-fastqc
~~~

### 02_ gunzip

I'll use this for loop to create a list of gunzip commands that will have a new output with the `.fq` extension AND it will keep the orginal `.fastq.gz` file.

~~~ {.bash}
for file in *.fastq.gz
do
	newfile=${file//.fastq.gz/.fq}
	echo $file $newfile
	cat >> 02_gunzip.cmds <<EOF
	gunzip -c $file $newfile
EOF
done
~~~

Then, I used this launcher command to unzip the files. I removed the job stderr and stdout files after the job was complete with `rm gunzip.*`. 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 02_gunzip.cmds -n gunzip -l 02_gunzip.slurm -A NeuroEthoEvoDevo
sbatch 02_gunzip.slurm
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
	cat >> 03_clean.cmds <<EOF
	fastx_clipper -v -a AAAAAAAA -l 20 -Q33 -i $file | fastx_clipper -v  -a TTTTTTTT -l 20 -Q33 | fastq_quality_filter -v  -Q33 -q 20 -p 90 -o $newfile 
EOF
done
~~~ 

Then, I created and launched a job with this command.

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 03_clean.cmds -n clean -l 03_clean.slurm -A NeuroEthoEvoDevo -m 'module load fastx_toolkit/0.0.13.2'
sbatch 03_clean.slurm
~~~ {.bash}

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
     cat >> 04_fastqc.cmds <<EOF
fastqc $file
EOF
done
~~~ {.bash}

Create and launch the slurm file.

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 04_fastqc.cmds -n fastqc -l 04_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 04_fastqc.slurm
~~~

Used filezilla to transfer the files to my local computer and then check them online.


