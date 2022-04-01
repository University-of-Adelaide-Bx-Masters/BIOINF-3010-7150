# Week 6 Practical Part 1: **Genome Assembly** 
## By Chelsea Matthews 
{:.no_toc}

* TOC
{:toc}

Today we will be looking at *de novo* assembly with more of a focus on larger, more complex genomes.  

# Get Started!

As with previous weeks' practicals, you will be using RStudio to interact with your VM.

While the focus of this section of the coursework is on the assembly of large complex genomes, it's not feasible for us to actually assemble a large complex genome because, frankly, our VM's are too small and sitting around for hours waiting for an assembly to finish is boring. 
So, we'll be assembling a small genome (Campylobacter hepaticus) from PacBio reads.  
Even though the Campylobacter hepaticus genome is fairly small, it will still take ~10-15 minutes for the assembly to run so we'll get it going straight away. 

## 1. Set up

Planning your directory structure from the beginning makes your work easier to follow for both yourself and others. 
The following directory structure is recommended for this practical. 
Note that I've used a numbered prefix for the `0_data` directory. 
I like to number directories within the project directory in the order in which they are created. 
This is a personal preference but I find that it makes it easier to navigate as I can remember the order I did things in but not always the specific directory name. 
For this practical we will use a numbering system so that you can see if you like it. 

```bash
mkdir --parents ~/Project_6/0_data/
cd ~/Project_6/
```

We will be using PacBio CLR reads for our assembly today. Copy the appropriate data from the `~/data` directory to your local folder.

```
cp ~/data/project_6/c_hepaticus.fastq ./0_data/
```

## 2. Create conda environments

We will also need a new conda environment.
To create our conda environment, we will first install mamba and then use it to create an environment using a configuration file. 
Using a configuration file means that we can easily track which versions of tools we are using. 

First, copy over the config file from `~/data/`.

```
cp -r ~/data/project_6/genomeassembly.yaml .
```
 
Have a quick look at the contents of the .yaml file to see which tools we will be installing in the environment. 

```
# install mamba 
conda install -c conda-forge mamba

# create the conda environment
mamba env create -f genomeassembly.yaml
```

This hopefully won't take too long. 

Now, when you run `conda env list`, a new environment should appear called `genomeassembly`.

We will now activate the `genomeassembly` environment. 

```
conda activate genomeassembly
```

## 3. Quick look at the input data

Today we are assembling a Campylobacter hepaticus genome from PacBio CLR reads using Canu.
Campylobacter is a bacteria and is one of the most common causes of baterial gastroenteritis in humans.
It was selected for the practical because of its small size (~1.5Mbp) to reduce computation time.

We will use seqtk to get a quick summary of the fastq data. 
We could alternatively use FastQC for this but the first three lines of the seqtk report are all we really want as PacBio reads don't have quality scores.

```
seqtk fqchk 0_data/c_hepaticus.fastq > raw_fastq_summary.txt
head -n 3 raw_fastq_summary.txt
``` 

* *How many reads are there in total in the `c_hepaticus.fastq` file?*
* *What are the minimum and maximum read lengths?*
* *Assuming that our genome is approximately 1.5Mbp long, what coverage do these reads give us?*
* *Do you think this is sufficient to generate a good assembly? Why or why not?*

## 4. Assemble!

We will be using the assembler Canu today. Canu was installed as part of the `genomeassembly` conda environment. 
Test that it works by running `canu`.
The Canu usage instructions should be printed to the screen.

Now let's run the Canu assembly. 

```bash
# run canu assembly
canu -p c_hepaticus -d 1_canu_assembly genomeSize=1.5m corThreads=2 -pacbio 0_data/c_hepaticus.fastq
```

It should take between 10 and 15 minutes for the assembly to complete.   

## 5. While we wait...

### Resources for genome assembly

Larger, more complex genomes contain high percentages of repeats. Assembling these repeats accurately with short reads is not effective. Long reads that are able to span these repeat regions result in much better assemblies. There are many long read assembly tools available including the following: 

* Canu (which we are using today)
* Flye
* NextDenovo
* Wtdbg2 (pronounced "redbean" - I don't understand this either)
* Raven
* Falcon
* Shasta
* miniasm

These tools all do more or less the same job (ie. long read assembly) but with slight differences. 

When we are assembling a small genome, memory (RAM) and CPU (compute threads/cores) requirements are relatively low and so we can choose our assembler based almost entirely on the quality of the resulting assembly. 
While we can still do this to some degree for larger genomes, resource allocations and limitations in High Performance Computing (HPC) environments mean that some tools may not be suitable for our particular dataset (we may not be able to run the assembly to completion due to time or resource limitations) and we need to more carefully consider the parameters we use for the assembly so that we don't need to re-run it unnecessarily. 
Even without considering the cost of wasted resources (and this can be measured in dollars where we are using a service like Amazon Web Service (AWS)), re-running these assemblies can take a long time. 

A number of tools have been designed to enable us to assemble larger genomes where resources are limited and Canu is one of these tools.
Read the [Canu Quick Start](https://canu.readthedocs.io/en/latest/quick-start.html#quickstart) summary to see what kind of features Canu has. 

One of the difficulties with assembling larger genomes is working out what resources a tool will require and deciding whether it will be suited to assembling your genome of interest within the resources available to you. 
Have a look at the [Canu FAQ](https://canu.readthedocs.io/en/latest/faq.html) and see what it says about the resources required (specifically CPU hours and memory) for bacterial and mammalian assemblies (it's the first question in the FAQ) and answer the following questions:

* *How many CPU hours (approximately) are required to assemble a bacteria?*
* *How many CPU hours (approximately) are required to assemble a mammal with high raw-read coverage?*
* *How many CPU hours (approximately) are required to assemble a mammal from HiFi reads?*
* *Why is the number of CPU hours and memory requirement lower for a HiFi assembly compared with a a non-HiFi read assembly?*
* *If you had access to a single node with 48 compute threads and sufficient memory, how long (in hours) would it take to assemble a "well-behaved large genome"?* 
* *How long would this same assembly take on your local VM (assuming memory was not a limitation) with two compute threads?*

You can see now why we aren't assembling a larger genome. 

### Canu

As it says above, Canu has been designed to run where resources are somewhat limited. 
It does this by breaking the assembly process down into many steps (hundreds when the genome is very large) and submitting each small job separately to the scheduler. 
While it may submit hundreds of small separate jobs, Canu effectively runs in three steps.

They are: 

1. Correction - Increases the accuracy of bases in all reads
2. Trimming - Remove low quality bases and discard suspicious sequences (for example, the PacBio Adapter if PacBio reads are used)
3. Assembly - Reads are ordered based on their overlaps to create contigs and produce the final consensus sequences - this step implements the OLC algorithm

Look through the [Canu documentation](https://canu.readthedocs.io/en/latest/tutorial.html) to answer the following questions:

* *What do the following options used in our Canu assembly command mean?*

-p 

-d 

genomeSize

corThreads

-pacbio

* *What is the "correctedErrorRate" default parameter for PacBio reads?*
* *Why is this parameter higher for Nanopore reads than for PacBio reads?*
* *Given that we have only approximately 16x coverage of reads for the C. hepaticus genome, is this correctedErrorRate still appropriate? If not, what should it be adjusted to and why?* 
* *Canu requires a genomeSize parameter. What does it use this parameter for?*
* *What is the default minimum read length (minReadLength)? Hint - you may need to look in a different section of the Canu documentation.*

## 6. Looking at the assembly report

Canu should have finished running by now so let's see what it's done by having a look at the summary report it has prepared.  

```
cd 1_canu_assembly

# have a look at the contents of the directory created by Canu 
# c_hepaticus.contigs.fasta is the final assembly
ls

# page through the assembly report
less c_hepaticus.report
```

The report should begin with a heading that says [CORRECTION/READS]. Remember that read correction is the first step in the Canu workflow. Take a look through this first section and answer the following questions. 

* *How many reads did Canu find in the raw data?*
* *What depth of coverage does that equate to?*
* *Why is the number of reads different to the number of reads in the input .fastq file?*

After correction, Canu trims the data and then the final step - unitigging - starts. 

Scroll down to the UNITIGGING section. Unitigging is the assembly stage and a unitig can be thought of as a high-confidence contig. 

* *How many reads make it through to the unitigging step?*
* *How many contigs does the final assembly have?*
* *How long is the final assembly and how does this compare with the estimated size of the genome?*

When you are finished, exit out of the report with `q`. 

We will have a closer look at the assembly quality in the next practical. 

## 7. Another Canu assembly

Re-do your Canu assembly with the same options as before but add in an adjustment of the correctedErrorRate parameter to the value that you suggested in Part 5. 
The output directory should be 1_canu_adjusted_assembly. 
Leave the other options the same. 

We will use this assembly in the next practical. 

## 8. Diploid assembly

C. hepaticus is a small haploid genome which means that it is fairly straightforward nowadays to assemble. 
What is more challenging are larger diploid and polyploid genomes.

* *Can you think of a few reasons why diploid or polyploid genomes are more difficult to assemble than haploid genomes? Relate these reasons back to the OLC algorithm.*
* *If we were assembling a diploid genome with a genome size of 2Gbp, under what circumstances might our resulting assembly be approximately 2Gbp in length?*
* *With the same assumptions, under what circumstances might our resulting assembly be greater than 2Gbp in length?* 

Canu is capable of assembling both haplotypes of a diploid genome - termed phasing. 
Have a look at the paper below (mainly the first figure) to understand how this works.

[Haplotype-Resolved Cattle Genomes Provide Insights Into Structural Variation and Adaptation](https://www.biorxiv.org/content/10.1101/720797v3.full)

* *Briefly summarise in your own words how trio binning works in Canu.*


