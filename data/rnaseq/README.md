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

## RNAseq Project and Jobs

All of the sample information is contained in `../data/sample_info/punches.csv` and `../data/sample_info/punches.csv`. 

These samples have been process in multiple sequencing jobs. 

| RNAseq Job | Data | Description
| :--- | :---: | :--- |
JA16268 | May 24, 2016 | paired end
JA16443 | July 26, 2016 | Tag-seq
JA16444 | October 5, 2016 | paired end
JA17009 | TBD | paired end

Rather than set hard paths to project-specific directories, I will set enviornment varaiable that specific the folders for a given project. This way, I just need to update the variable to rerun the workflow, rather than continually modify the path name. 

~~~ {.bash}
RNAseqProject=<insernameofmainproject>  ## replace with project name (i.e. BehavEphysRNAseq) 
RNAseqJob=<insernaseqjobname>  ## replace with job number (i.e. JA16268)

# now we can use these variablles to call the path
mkdir -p $SCRATCH/$RNAseqProject/$RNAseqJob
~~~ 

For this project, $RNAseqProject will always be "BehavEphysRNAseq" but $RNAseqJob will be either JA16268 or JA16444

## Workflow 

I use the command line to create a lot of the commands files and batch scripts, so this readme documents the workflow I use to create these files.  My hope is that this style will help me easily reproduce the workflow (either for reanalysing this project or for applying the workflow to new datasets).

I decided to keep my data and my scripts and my commands files in the same directory, gasp! I tried to keep them separate for a while, but having the commands files in the same directory with the launcher script and the data just makes life easier. I do move the fastqc output files to their own directory after processing because that makes me happy.

##  00_rawdata

Inside the project directory on scratch `$SCRATCH/BehavEphyRNAseq`, create subdirectory for the job and then a subdirectory for the raw data. Navidate to the directory

~~~ {.bash}
mkdir -p $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
cd $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
~~~ 

To download data from the GSAF, download and run the script provided on their website.
To download, I'll copy a  bash script named `00_gsaf_download.sh` from a one of my folders. This bash script was provide by the GSAF (https://wikis.utexas.edu/display/GSAF/How+to+download+your+data). 

~~~ 
To Download the data:
1. Visit https://wikis.utexas.edu/display/GSAF/How+to+download+your+data
2. Copy the bash script and save it to a file called `gsaf_download.sh` as suggested.
3. Copy the amazon key that corresponds to your data. It will look something like `http://gsaf.s3.amazonaws.com/...`
4. Execute the bash script with the following synatx: `gsaf_download.sh "amazon key"`
~~~ 

Once you've downloaded the data, let's write-protect all the files so we don't accidentally modify them.

~~~ {.bash}
chmod u-w *.fastq.gz
~~~ 

Then, I like to move all the log files to a subdirectory so they are out of the way. This isn't necessary, but it is a way to keep the files that might come in handy in the future. 

~~~ {.bash}
mkdir wgetlog
mv *.wget.log wgetlog/
mv files.html wgetlog/
mv md5.txt wgetlog
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

Then, I use the python function `launchercreator.py` to prepare and submit a batch job to process all the files at once. For details on the arguments, see: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j 01_fastqc.cmds -n 01_fastqc -l 01_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 01_fastqc.slurm
~~~

Once the job completes, I moved all the output files to a new directory on the same levels as the raw data.

~~~ {.bash}
mkdir ../01_fastqc
mv *.html ../01_fastqc  # this is the html file for easy viewing
mv *.zip ../01_fastqc   # this allow command line viewing of results
mv 01_fastqc.e* ../01_fastqc  #TACC standard error file
mv 01_fastqc.o* ../01_fastqc	#TACC standard output file
~~~

One must save the data locally in order to view the html files. 

In a new terminal window:

~~~ 
1. Create a directory on your local machine. 
2. Navigate to this directory
3. Use the scp command to copy the files `scp -r <username>@stampede.tacc.utexas.edu:/scratch/02189/rmharris/BehavEphyRNAseq/JA16444/01_fastqc/ .` # be sure to include the period
~~~

### 02 Filter and Trimming Reads (filtrimmedreads)

Next we use cutadapt to quality filter reads in trim adaptors. The cutadapt program is not maintained by TACC, so we need to use the program stored in a colleagues working directory. To access the program, set the path.

~~~ {.bash}
PATH=/work/01184/daras/bin/cutadapt-1.3/bin:$PATH 
~~~

Again, we will use a for loop to create the commands file. For details on Cutatdapt, visit: http://cutadapt.readthedocs.io/en/stable/guide.html.

~~~ {.bash}
cd $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
for R1 in *R1_001.fastq.gz
do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    R1filtrim=$(basename $R1 fastq.gz)filtrim.fastq.gz
    R2filtrim=$(basename $R2 fastq.gz)filtrim.fastq.gz
    echo $R1 $R2 $R1filtrim $R2filtrim
    echo "cutadapt -q 15,10 -a GATCGGAAGAGCACACGTCTGAACTCCA -A ATCGTCGGACTGTAGAACTCTGAACGTG -m 22 -o $R1filtrim -p $R2filtrim $R1 $R2" >> 02_filtrimmedreads.cmds
done
~~~

Then, launch a batch script. 

~~~ {.bash}
launcher_creator.py -t 4:00:00 -j 02_filtrimmedreads.cmds -n 02_filtrimreads -l 02_filtrimmedreads.slurm -A NeuroEthoEvoDevo -e rayna.harris@utexas.edu -q 'normal'
sbatch 02_trimreads.slurm
~~~

Once that's done, let's clean up a bit and move our files to a new directory. 

~~~ {.bash}
# moving just the data, leaving the output files where they are.
mkdir ../02_filtrimmedreads
mv *filtrim.fastq.gz ../02_filtrimmedreads
~~~


##  Building a Kallisto index of the mouse transcriptome
An Kallisto index must be created for mapping and counting reads with Kallisto. This only have to be done one. See https://pachterlab.github.io/kallisto/manual for details on using Kallisto.

For the mouse, the reference transciptome can be downloaded https://www.gencodegenes.org/mouse_releases/current.html. For this project, I downloaded version M11.

~~~ {.bash}
mkdir $SCRATCH/$RNAseqProject/refs
cd $SCRATCH/$RNAseqProject/refs
curl -O ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.pc_transcripts.fa.gz
~~~

Then, create the commands file to build a kallisto-specific index.

~~~ {.bash}
echo "kallisto index -i gencode.vM11.pc_transcripts_kallisto.idx $SCRATCH/$RNAseqProject/refs/gencode.vM11.pc_transcripts.fa.gz" > kallisto_index.cmds
cat kallisto_index.cmds
~~~

Then create the launcher script. Again, Kallsito isn't a standard module on TACC, so we must specific the path to Kallisto stored in a colleagues work directory.

~~~ {.bash}
launcher_creator.py -t 0:30:00 -j kallisto_index.cmds -n kallisto_index -l kallisto_index.slurm -A NeuroEthoEvoDevo -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch kallisto_index.slurm
~~~


### 03_kallistoquant

Okay, now we are ready to use our index for mapping and quantification. See https://pachterlab.github.io/kallisto/manual for details on using Kallisto.

~~~ {.bash}
# Let's navigate back to where our trimmed reads are
cd $SCRATCH/$RNAseqProject/$RNAseqJob/02_filtrimmedreads
~~~ 

Create the commands file. We also need to create the directory where these files will go so they don't override each other. 

~~~ {.bash}
mkdir $SCRATCH/$RNAseqProject/$RNAseqJob/03_kallistoquant
for R1 in *R1_001.filtrim.fastq.gz
do
    R2=$(basename $R1 R1_001.filtrim.fastq.gz)R2_001.filtrim.fastq.gz
    samp=$(basename $R1 _R1_001.filtrim.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -i $SCRATCH/$RNAseqProject/refs/gencode.vM11.pc_transcripts_kallisto.idx  -o ../03_kallistoquant/${samp} $R1 $R2" >> 03_kallistoquant.cmds
done
~~~

Create the launcher script and run. **Note**: We need to use the largemem cluster because this program requires lots of working memory! The largemem cluster has a max core limit that is easy to reach. If you have too many samples, split the job into pieces. Note: This cluster has long wait times! Check first 

~~~ {.bash}
launcher_creator.py -t 1:00:00 -j 03_kallistoquant.cmds -n 03_kallistoquant -l 03_kallistoquant.slurm -A NeuroEthoEvoDevo -q largemem -m 'module use -a /work/03439/wallen/public/modulefiles; module load gcc/4.9.1; module load hdf5/1.8.15; module load zlib/1.2.8; module load kallisto/0.42.3'
sbatch 03_kallistoquant.slurm
~~~

Now, save the data locally for analysis. 

~~~ 
1. Make a local directory for storing the files. 
2. Navigate to that directory.
3. Copy the files with `scp -r rmharris@stampede.tacc.utexas.edu:/$SCRATCH/$RNAseqProject/JA16444/02_kallistoquant .
~~~

Now, I need to do some processing so that the file name corresponds to the sample name for downstream data analysis. Let's remove the uninformative bits. The GSAF has added `_S` followed by some numbers to each sample, so let's remove that. Note: Maybe check that this will indeed do the correct renaming before actually renaming the file.  

~~~ {.bash}
for file in *
do
    sample=${file//_S*/}
    echo $file $sample
    mv $file $sample
done
~~~

Also, all my sample names have dashes, but the file names have underscores, so let's replace the `_` with `-`

~~~ {.bash}
for file in *
do
    sample=${file//_/-}
    echo $file $sample
    mv $file $sample
done
~~~

### Now, we are ready for R!!!




