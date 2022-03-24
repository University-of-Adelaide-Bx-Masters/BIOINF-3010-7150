# Week 5 Practical Part 1: Short and Long Read Alignment
## By Nathan Watson-Haigh, updated/adapted by Dave Adelson
{:.no_toc}

* TOC
{:toc}

As with previous week's practicals, you will be using RStudio to interact with your VM.

See [week 1's practical](../Bash_Practicals/1_IntroBash.md#rstudio) to remind yourself how to connect to your VM.

# Setup for today

## Working Directory

First we will set up a directory for today's practical.
In general it is very worthwhile to keep all your project-specific code and data organised into a consistent location and structure.
This are not essential, but is very useful and is in general good practice.
If you don't follow this step, you will be making your life immeasurably harder for the duration of this practical.

To make and enter the directory that you will be working in, run the following commands in the terminal pane.

```bash
# Setup project working directory
mkdir --parents ~/Project_4/data/
cd ~/Project_4/

# load the required software environment
conda activate assembly
```

## Get the Data

```bash
# Make the directory for the reference genome
mkdir --parents ~/Project_4/data/reference/

# Make subdirectories for the various data sets
mkdir --parents ~/Project_4/data/{illumina_pe,pacbio}/

# Get the data
#####
# RefSeq E. coli K-12 substr. MG1655
cp ~/data/S_LR_Alignment/NC_000913.3.fasta.gz ~/Project_4/data/reference/
# Illumina PE
cp --link ~/data/S_LR_Alignment/36_ACGCACCT-GGTGAAGG_L002_R?_001_40x.fastq.gz ~/Project_4/data/illumina_pe/
# PacBio
cp --link ~/data/S_LR_Alignment/lima.bc1106--bc1106_40x.subreadset.fastq.gz ~/Project_4/data/pacbio/
# Nanopore
# TODO
```

### Questions

 - *What does the `--link` argument to `cp` do? Hint: Use the `man` page and google to work it out.*
 - *Why might using `--link` be useful with genomics data files?*

# QC

Using FastQC, check the Illumina data and determine if you need to perform read trimming.

```bash
fastqc \
  --threads 2 \
  data/illumina_pe/*.fastq.gz
```

## Questions

 - *What is the minimum read length present in the Illumina data? Hint: Look at the FastQC reports.*
 - *What do you think of the "Per base sequence content" plot in the FastQC reports?*

## Read Trimming and Filtering

If you deem it necessary to quality trim, adapter trim or length filter your raw reads then look back at last week's code and use either Trimmomatic or fastp to perform the trimming/filtering.

# Illumina Read Mapping

Once you're happy your reads are good to align to the reference genome

```bash
# Index the reference genome for use with BWA
bwa index data/reference/NC_000913.3.fasta.gz

# Create an output directory for read mappings
mkdir --parents mappings/

# Align reads
time bwa mem \
  -t 2 \
  -T 30 \
  data/reference/NC_000913.3.fasta.gz \
  data/illumina_pe/36_ACGCACCT-GGTGAAGG_L002_R1_001_40x.fastq.gz \
  data/illumina_pe/36_ACGCACCT-GGTGAAGG_L002_R2_001_40x.fastq.gz \
| samtools view -F 4 -u \
| samtools sort \
  --threads 2 -l 7 \
  -o mappings/36_ACGCACCT-GGTGAAGG_L002_40x_T30.bam \
  /dev/stdin
```

## Questions

 - *What does the `-T` argument to `bwa mem` do and what default value does `bwa mem` use?*
 - *What does the `-F 4` argument to `samtools view` do?*
 - *What does the `-u` argument to `samtools view` do and why is this beneficial in the middle of a pipeline?*
 - *What does the `-l 7` argument to `samtools sort` do?*

# PacBio Read Mapping

```bash
# Index the reference genome for minimap2
minimap2 \
  -d data/reference/NC_000913.3.fasta.gz.mmi \
  data/reference/NC_000913.3.fasta.gz

# Align the reads
time minimap2 \
  -ax map-pb \
  -t 2 \
  data/reference/NC_000913.3.fasta.gz \
  data/pacbio/lima.bc1106--bc1106_40x.subreadset.fastq.gz \
| samtools view -F 4 -u \
| samtools sort \
  --threads 2 -l 7 \
  -o mappings/lima.bc1106--bc1106_40x.bam
```

## Questions

 - *Have a look at the contents of the read file. What do you make of the quality strings?*

# IGV

As with the last practical we will use [IGV-web](https://igv.org/app/) to visualise the `.bam` files. 

Decompress the E. coli K-12 genome and index it with `samtools faidx`
Index the `.bam` files with `samtools index`

Download the reference, the `.bam` files and the index files to your computer.

load the E. coli K-12 reference and index and the `.bam` files and index files as in the previous practical on [IGV-web](https://igv.org/app/).

What do you make of these regions:

 * `NC_000913.3:276291-294244`
 * `NC_000913.3:4295777-4296810`

**Optional exercise**

If you want to try to download and run IGV on your local computer, head to the [IGV download page](https://software.broadinstitute.org/software/igv/download) and grab the version for your operating system. This will only work if you are installing on your personal computer, not one of the computer suite workstations. 
If you don't have Java or you don't want to install the Java version that comes with it, or have troubles running it, just use [IGV-web](https://igv.org/app/) instead.

**NOTE: If you're using IGV on your own machine, you will not need to decompress the genome sequence FASTA file or load `.bai` files.**

# BWA Alignments

In the above `bwa mem` alignment step, we used the default value (`30`) for `-T`.
I asked you what the option does.
Now I want you to explore the effect of increasing this value on the number of mismatches observed in the read alignments when you load the resulting BAM files into IGV.
