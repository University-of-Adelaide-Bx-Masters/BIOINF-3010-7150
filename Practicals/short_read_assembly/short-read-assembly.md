# Week 3 Practical Part 2: SARS-CoV-2 Short Read Assembly
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
mkdir --parents ~/Project_3/data/
cd ~/Project_3/

# load the required software environment
conda activate assembly
```

## Get the Data

```bash
# Make the directory for the reference genome
mkdir --parents ~/Project_3/data/reference/

# Get the SARS-CoV-2 reference genome
cp ~/data/COVID-19.fasta.gz ~/Project_3/data/reference/

# Make the directory for the Illumina PE reads
mkdir --parents ~/Project_3/data/illumina_pe/

# Get the subsampled Illumina PE data
cp ~/data/SRR111407{44,46,48,50}_?_*x.fastq.gz ~/Project_3/data/illumina_pe/

# Make the directory for scripts
mkdir --parents ~/Project_3/data/scripts/

# Get the scripts
cp ~/data/plot_delta.R ~/Project_3/data/scripts/
```

Lets see what files and directories we have under our current working directory:

```bash
tree
```

### Questions

 - *What did `{44,46,48,50}` mean/do in the above `cp` command while grabbing the Illumina PE data?*

## Data

We're going to be using some very topical data, SARS-CoV-2 data, the virus causing COVID-19.
We will be looking at this data in different ways, providing insights into how bioinformaticians analyses these types of data.

### RefSeq Genome Assembly

We will use the SARS-CoV-2 genome assembly [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).

This assembly was generated from a sample collected in Dec 2019 and submitted to NCBI on 13th Jan 2020 as [NC_045512.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.1).
Subsequently, [NC_045512.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.1) was replaced by [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) on the 17th Jan 2020.

### Public Sequencing Data

Public sequence data is being released via NCBI's [SARS-CoV-2](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) page.

We will be looking at some data released on [21st Feb 2020](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP250294) for a clinical swab obtained from a confirmed case in Madison, WI, USA.
For further information see: [https://openresearch.labkey.com/wiki/ZEST/Ncov/page.view?name=SARS-CoV-2%20Deep%20Sequencing](https://openresearch.labkey.com/wiki/ZEST/Ncov/page.view?name=SARS-CoV-2%20Deep%20Sequencing)

We will be working with the `SRR11140748` sample (Illumina data) but you are welcome to also look at any of the other 3 samples if you have time.
Here is a table of information linking to the orginal source of the data:

| Description  | Accession/URL                                                                                | Coverage |
|:-------------|:--------------------------------------------------------------------------------------------:|---------:|
| RefSeq       | [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)                 |          |
| **Illumina** |                                                                                              |          |
| veroSTAT-1KO | [SRR11140744](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11140744)                    |    7,586 |
| veroE6       | [SRR11140746](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11140746)                    |    5,328 |
| vero76       | [SRR11140748](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11140748)                    |    6,353 |
| swab         | [SRR11140750](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11140750)                    |      252 |

# Short Read Assembly

## Initial Goals

 1. Assemble a genome using a de Bruijn graph assembler (SPAdes)
 2. Visualise assembly graphs
 3. Compare assemblies against an existing reference (MUMmer)


# De Novo Genome Assembly

Lets perform some *de novo* genome assemblies of Illumina reads.
I encourage you to explore the effects of using differing amounts of coverage of the input data on the assembly and explore the assembly of the other 3 samples.
We have made data files available for `10x`, `20x`, `40x`, `80x` and `100x` coverage of all 4 Illumina samples.

```bash
time spades.py \
  --threads 2 \
  -o ~/Project_3/de_novo_illumina/SRR11140748_10x_PE \
  -1 ~/Project_3/data/illumina_pe/SRR11140748_1_10x.fastq.gz \
  -2 ~/Project_3/data/illumina_pe/SRR11140748_2_10x.fastq.gz \
| tee ~/Project_3/de_novo_illumina/SRR11140748_10x_PE.log
```
## Questions

 - *We used `tee` in the above pipeline; what did it do?*
 - *Look at the log file generated, what k-mer sizes were used for the assembly?*

## Comparisons to Reference Genome

Compare the assembly to the SARS-CoV-2 RefSeq assembly using [MUMmer](http://mummer.sourceforge.net/manual/#nucmer)'s `nucmer` tool:

```bash
# Install MUMmer
conda install -c bioconda mummer

# Decompress the reference genome for MUMmer
pigz -dcp2 \
  < ~/Project_3/data/reference/COVID-19.fasta.gz \
  > ~/Project_3/data/reference/COVID-19.fasta

# Run nucmer from MUMmer package
nucmer \
  -maxmatch \
  -minmatch 100 \
  -mincluster 500 \
  -prefix ~/Project_3/de_novo_illumina/SRR11140748_10x_PE \
  ~/Project_3/data/reference/COVID-19.fasta \
  ~/Project_3/de_novo_illumina/SRR11140748_10x_PE/contigs.fasta
```

This will generate a `.delta` file which describes the alignments between your assembled contig(s) and the reference sequence(s).
While the file is plain text, it's not meant for human consumption.
Instead, we need to use some other tools to visualise the information it contains.

### Visualisation of Delta Files

The online tool [Assemblytics](http://assemblytics.com/) can be used to generate useful plots of the information contained within a `.delta` file.

Alternatively, your can use some R code to plot the information.
We have provided an R script for this purpose, which can be executed as follows:

```bash
Rscript --vanilla \
  ~/Project_3/data/scripts/plot_delta.R \
  ~/Project_3/de_novo_illumina/SRR11140748_10x_PE.delta
```

This will create a PDF file called `~/Project_3/de_novo_illumina/SRR11140748_10x_PE.delta.pdf`

### Questions

 - *How many contigs are in your assembly?*
 - *Do those contigs cover all the bases of the reference genome?*
 - *Are there any variants detected between this sample and the reference genome? NOTE: you may need to modify Assemblytics settings to see small variants.*

# Explore k-mer Length Effects

SPAdes uses several k-mer lengths during the assembly process to produce a de Bruijn graph for each k-mer length and then combines the information to produce a final assembly.
Using what you learned above about exploring assemblies, lets look at how an assembly is affected by our choice of k-mer length.

First, lets look at an assembly produced with a single short k-mer of length of 11:

```bash
k=11

time spades.py \
  --threads 2 \
  -k ${k} \
  -o ~/Project_3/de_novo_illumina/SRR11140748_10x_PE-k${k} \
  -1 ~/Project_3/data/illumina_pe/SRR11140748_1_10x.fastq.gz \
  -2 ~/Project_3/data/illumina_pe/SRR11140748_2_10x.fastq.gz \
| tee ~/Project_3/de_novo_illumina/SRR11140748_10x_PE-k${k}.log
```

Next, lets look at an assembly produced with a long k-mer of length of 127:

```bash
k=127

time spades.py \
  --threads 2 \
  -k ${k} \
  -o ~/Project_3/de_novo_illumina/SRR11140748_10x_PE-k${k} \
  -1 ~/Project_3/data/illumina_pe/SRR11140748_1_10x.fastq.gz \
  -2 ~/Project_3/data/illumina_pe/SRR11140748_2_10x.fastq.gz \
| tee ~/Project_3/de_novo_illumina/SRR11140748_10x_PE-k${k}.log
```

Now modify the above commands to produce another 2 assemblies but using the `100x` coverage data.

## Questions

 - *How many contigs were produced with the above 4 assemblies?*
 - *Compare each of the assemblies with the reference genome using MUMmer and Assemblytics and the R script provided.*

