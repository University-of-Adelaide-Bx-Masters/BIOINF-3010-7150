# BIOINF3010/7150 Genome Assembly - Assignment 2 (worth 20% of your final mark)

## Introduction

You are provided with four fastq files containing PacBio Continuous Long Reads which were subset from Whole Genome Sequencing of Lactobacillus acidophilus. 
The original dataset can be found [here](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18189666).
Lactobacillus acidophilus is a common probiotic that is an ingredient in yoghurt and other fermented foods.

You will be generating an assembly from each dataset and then comparing the quality of these assemblies. 

## Data

Your datasets are as follows:

* 20x.fq - random subset of original dataset to ~20x coverage 
* 20x_filt.fq - random subset of original dataset to ~20x coverage of reads greater than 10,000 bases in length 
* 25x.fq - random subset of original dataset to ~25x coverage
* 35x.fq - random subset of original dataset to ~35x coverage

These files can be found in the `/data/assignment2` directory

You are also provided with the Lactobacillus acidophilus reference (taken from [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_020883435.1) if you want to have a look). 

It is also in the `/data/assignment2` directory and is called `reference.fasta`. 

### Part 1

Canu requires a `genomeSize` parameter in order to run. 

**1.1** What genome size will you use? 
Justify this choice with a source.

### Part 2

**2.1** Run FastQC on each dataset and inspect the sequence length distributions. 
In addition, find out the following for each dataset: 

* total number of bases
* number of reads
* approximate genome coverage calculated using your chosen genome size

**2.2** Predict which datasets will produce the best and worst assemblies. 
Don't worry if your predictions don't match up with your results later. 
Just try to justify them based on the information you've collected and your current knowledge. 

### Part 3

**3.1** Assemble all four datasets using Canu. Document the parameters you used. 

### Part 4

**4.1** Compare these assemblies with each other using QUAST and provide an image of the contig length distribution. 

**4.2** Comment on the contig length distribution. Is this what you expected? 

**4.3** Explain why contiguity isn't a good measure of assembly accuracy but is still relevant to the overall assessment of assembly quality. 

### Part 5

**5.1** Run BUSCO on your four assemblies using an appropriate lineage and create a comparison image using the `generate_plot.py` script. Note that this script will only work when the `quality` conda environment is activated.
Provide this image. 

**5.2** Which lineage did you choose?

**5.3** Comment on the BUSCO results. Which assembly appears to be the best and which is the worst? 

**5.4** What does the presence of duplicated genes likely indicate in the assembly of a haploid genome?

### Part 6

**6.1** Compare your four Canu assemblies with the reference using QUAST and have a look through the report. 
Comment on the "# indels per 100 kbp" metric remembering that an indel is a small insertion or deletion in the assembly that isn't present in the reference. 
How does coverage impact this metric? 

**6.2** Would you expect that an assembly created from 200x read coverage would have a "# indels per 100kbp" value of 0? Explain why or why not.

### Part 7

**7.1** Considering your findings so far, justify which assembly you think is the best and which is the worst.
If your findings don't match your predictions from Part 2, try to explain why this might be.

