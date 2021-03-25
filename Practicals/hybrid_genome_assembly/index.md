# Week 4 Practical Part 2: Hybrid Genome Assembly
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
mkdir --parents ~/Project_5/data/
cd ~/Project_5/

# load the required software environment
conda activate assembly
```

## Get the Data

```bash
# Make subdirectories for the various data sets
mkdir --parents data/{reference,scripts,illumina_pe,pacbio}/

# Get the data
#####
# RefSeq E. coli K-12 substr. MG1655
cp ~/data/NC_000913.3.fasta.gz data/reference/
# Illumina PE
cp --link ~/data/36_ACGCACCT-GGTGAAGG_L002_R?_001_*x.fastq.gz data/illumina_pe/
# PacBio
cp --link ~/data/lima.bc1106--bc1106_*x.subreadset.fastq.gz data/pacbio/
# R script
cp --link ~/data/plot_delta.R data/scripts/
```

# Hybrid Genome Assembly

Using what you have learned previously, regarding read quality control, trimming and assembly tools, process the available Illumina PE and PacBio reads to assemble the E. coli K-12 substr. MG1655 genome.

Use the tool SPAdes (`spades.py`) to perform a _de novo_ assembly from **BOTH** Illumina and PacBio reads.
You will need to view the SPAdes [documentation/manual](https://github.com/ablab/spades/blob/spades_3.15.2/README.md#pacbio) to figure out how to specify both Illumina PE and PacBio reads in the same assembly.

As before, there are several subsamplings (`2`, `4`, `5`, `8`, `10`, `20` and `40`) of the original reads for both the Illumina and PacBio data.
When performing the assembly, consider how much of each data type you might need/use to generate a good assembly in reasonable time.

**HINT: time how long each `spades.py` command takes to complete.**

## Questions

 - *Increasing the number of input reads results in additional computational time requirements. If you double the input coverage, how much more time might be required to perform an assembly?*

## Comparing Assemblies

You should use what you learned about comparing an assembly against a reference genome to ascertain if an assembly is reasonable or whether more/less data might be required.

Most of these hybrid assemblies will take more than 10mins each, so you won't be able to do many within the time limit of this practical.

### Questions

 - *How many contigs and scaffolds do you get for your assemblies?*

Given the limited amount of time in this practical:

 - *What sort of coverage would be wise to start with high or low?*
 - *To explore assemblies with different mixes of Illumina/PacBio coverage, how as a group, might you be able to explore more assemblies with different levels of input coverage?*

## Visualisation

Visualise your assembly(ies) in IGV or IGV-web and look to see if you can see any inconsistencies between the aligned Illumina and PacBio data.

<!--

Sample command line for hybrid assembly.

```bash
ILLUMINA_COV=2
PACBIO_COV=2

mkdir --parents de_novo_hybrid

time spades.py \
  --threads 2 \
  -o de_novo_illumina/36_ACGCACCT-GGTGAAGG_L002_${ILLUMINA_COV}x_${PACBIO_COV}x \
  -1 data/illumina_pe/36_ACGCACCT-GGTGAAGG_L002_R1_001_${ILLUMINA_COV}x.fastq.gz \
  -2 data/illumina_pe/36_ACGCACCT-GGTGAAGG_L002_R2_001_${ILLUMINA_COV}x.fastq.gz \
  --pacbio data/pacbio/lima.bc1106--bc1106_${PACBIO_COV}x.subreadset.fastq.gz \
| tee de_novo_hybrid/36_ACGCACCT-GGTGAAGG_L002_${ILLUMINA_COV}x_${PACBIO_COV}x.log
```
-->
