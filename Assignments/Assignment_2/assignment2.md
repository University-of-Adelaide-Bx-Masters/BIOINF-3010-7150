# BIOINF3010/7150 Genome Assembly - Assignment 2 (36 marks which will be weighted to count for 20% of your final mark)

## Introduction

In this assignment, you aim is to de novo assemble one of the SMALLEST eukaryotic genomes, Encephalitozoon intestinalis. E. intestinalis belongs to Microsporidia, and it's a parasite (microbial fungi), which causes microsporidiosis (an oppotunistic intestinal infection that causes diarrhea and wasting in immunocompromised individuals, such as HIV). If you want to understand more about E. intestinalis, please hava a read at [wikipedia](https://en.wikipedia.org/wiki/Encephalitozoon_intestinalis). 

Although the genome of E. intestinalis is very small (~2.5 Mb), it has 11 chromsomes. If you want to know more about the genome statistics of E. intestinalis, please have a look at [here](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_024399295.1/).

## Data

You are provided with 8 fastq files containing sequencing reads using different sequencing platforms. These fastq files can be found in `2_assignment/raw_data` folder, and are including:

| File(s)                                    | Platform | Coverage | Description                                            |
|--------------------------------------------|----------|----------|--------------------------------------------------------|
| illumina_SR_30x_1.fq, illumina_SR_30x_2.fq | Illumina | ~30x     | Paried-end short reads from Illumina Miniseq           |
| nanopore_LR_15x.fq                         | Nanopore | ~15x     | Long reads from Nanopore MinION                        |
| nanopore_LR_15x_filt.fq                    | Nanopore | ~15x     | Long reads (reads >= 10 kb) from Nanopore MinION       |
| nanopore_LR_30x.fq                         | Nanopore | ~30x     | Long reads from Nanopore MinION                        |
| pacbio_LR_15x.fq                           | PacBio   | ~15x     | Long reads from PacBio_SMRT Sequel II                  |
| pacbio_LR_15x_filt.fq                      | PacBio   | ~15x     | Long reads (reads >= 10 kb) from PacBio_SMRT Sequel II |
| pacbio_LR_30x.fq                           | PacBio   | ~30x     | Long reads from PacBio_SMRT Sequel II                  |

The original dataset can be found [here](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=594722).

The illumina dataset will be used to do genome survey analysis (genome size estimation), and then you will be generating an assembly from each of six long reads datasets using Flye (v2.8.1) and then comparing the quality of these assemblies.

You are also provided with the E. intestinalis reference (taken from [here](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_024399295.1/) if you want to have a look). 

It is in the `/data/assignment2/DB` directory and is called `GCA_024399295.1_ASM2439929v1_genomic.fna`. 

In addition to the sequencing files, you will be also given two scripts which will be used to do sequence/genome statistics and genome survey analysis. These two scritps are stored in the folder of "2_assignment/bin", and are including:

| Script       | Description                                 | Link                                               |
|----------------|---------------------------------------------|----------------------------------------------------|
| assembly-stats | A light tool to do basic genome statistics  | https://github.com/sanger-pathogens/assembly-stats |
| genomescope.R  | A R script to do genome survey analysis     | https://github.com/schatzlab/genomescope           |

### Part 1, QC

**1.1** Run `assembly-stats` on each LR (Long Reads) dataset and find out the following info for each dataset:

* total number of bases *- 1 marks*
* number of reads *- 1 marks*
* average read length*- 1 marks*
* largest read length*- 1 marks*

**1.2** Predict which datasets will produce the best and worst assemblies. 
Don't worry if your predictions don't match up with your results later. 
Just try to justify your predictions based on the information you've collected and your current knowledge. *- 2 marks*

### Part 2

**2.1** Use `jellyfish` and `genomescope.R` to perform genome survey analysis. What is the estimated genome size? Provide the generated figure showing the fitted model for k-mer distribution (Hint: plot.png in your genomescope output foloder). *- 4 marks*

### Part 3

**3.1** Write a bash script to assemble all 6 LR datasets using Flye. [Hint: Assembly all 6 LR datasets will take ~50 mins in total, so be patient if you see the Flye is running for a while.] *- 6 marks*

**3.2** Run `assembly-stats` on each assembled genome and find out the following info for the assembled genomes:

* draft genome size *- 1 marks*
* number of contigs *- 1 marks*
* largest contig length*- 1 marks*
* N50 *- 1 marks*

### Part 4

**4.1** Compare these assemblies with each other using QUAST and provide the report (Hint: report.pdf in your QUAST output folder). *- 2 marks*

**4.2** Comment on the contig length distribution. Is this what you expected? *- 2 marks*

**4.3** Explain why contiguity isn't a good measure of assembly accuracy but is still relevant to the overall assessment of assembly quality. *- 2 marks*

### Part 5

**5.1** Run BUSCO on your 6 assemblies using an appropriate lineage and create a comparison image using the `generate_plot.py` script. Note that this script will only work when the `BUSCO` conda environment is activated. Provide this image [Hint: the "busco_figure.png" file in short summaries folder]. *- 4 marks*

**5.2** Comment on the BUSCO results. Which assembly appears to be the best and which is the worst? *- 2 marks*

### Part 6

**6.1** Considering your findings so far, justify which assembly you think is the best and which is the worst.
If your findings don't match your predictions from Part 1, try to explain why this might be. *- 4 marks*

## Check-list for your assignment2 submission

* A document including answers to all questions
* A image/figure showing the k-mer distribution model fitting on genome survey analysis in Part 2.1
* A bash script to assemble all 6 LR datasets using Flye in Part 3.1
* A report from the QUAST assessment in Part 4.1
* A comparison figure showing the BUSCO scores for all 6 draft assemblies in Part 5.1.
