# Week 6 Practical Part 2: Assembly Quality Assessment
## By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

This practical will continue on from Week 6 Practical 1. 

## What are we doing today?

In the first practical we assembled the C. hepaticus genome twice using Canu.
Today we will look at measuring *de novo* assembly quality with a focus on reference free quality metrics. 
In particular we will look at contig length distributions and BUSCO. 

## 1. Setting up

Copy over the pre-made assembly from `~/data` and load your conda environment for todays practical.

```
mkdir Project_6/1_flye_assembly
cp ~/data/project_6/assembly.fasta ./Project_6/1_flye_assembly/

# deactivate the current conda environment if you have one loaded and then activate the `genomeassembly` conda environment
conda deactivate 
conda activate genomeassembly

cd Project_6
```

Your `Project_6` directory should now contain three assembly directories, each containing a different assembly. 

```
1_canu_assembly
1_canu_adjusted_assembly
1_flye_assembly
```


## 2. Contig length distribution

You've previously learnt about N50 as a measure of assembly contiguity ([refresh your memory here](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50)) but N50 on its own only gives us a snapshot into assembly contiguity.
A better way to look at assembly contiguity is to inspect a cumulative contig length plot.
This is particularly helpful when we want to compare multiple different assemblies of the same genome or species. 

The tool QUAST (QUality ASsesment Tool) (documentation here)[http://quast.sourceforge.net/docs/manual.html] can be run on one or more assemblies at once and produces some handy comparison statistics and graphics including a cumulative contig length plot.
It can also be run with a reference genome to which it will compare the other provided assemblies and produce additional statistics.
Here we will run it without a reference. 

The basic command is below.
Make sure the names you provide to the `--labels` parameter are descriptive and that they match the order of the assembly.fa files.

```
quast -o 2_quast -t 2 --labels "name1, name2, name3" path/to/assembly1.fa path/to/assembly2.fa path/to/assembly3.fa
```

This won't take very long.

Open the `icarus.html` file located in the `2_quast` directory in a web browser. 
Have a look at the information provided at the two links on the web page and answer the following questions.

* *Spend some time interpreting the cumulative length plot. Under what circumstances would a cumulative length plot be more informative than just comparing N50s?*
* *What would an assembly with only one contig look like on this plot?*
* *Looking at the Contig size viewer, do you think it is more appropriate to compare N50 or NG50? Why? Find out what NG50 is if you're not sure.*
* *Based on the QUAST report, which assembly appears to be the best?* 

The image below is a slightly different version of a cumulative contig length plot. 
The assemblies represented are publicly available cannabis assemblies. 
The dashed horizontal lines are the approximate sizes of the male and female genomes with male being the longer.
Inspect the plot and answer the questions below. 

![Cumulative contig length plot](images/cum_contig_length_dist.png)

* *What benefit is there to having contig length on the x axis instead of the contig number in size order as in the QUAST plot?*
* *One of these assemblies was generated from short illumina reads. Which do you think it is?*
* *Can you tell which are contig level assemblies and which are scaffolded?*

Next, let's run BUSCO. 

## 3. BUSCO (Benchmarking Universal Single Copy Orthologs) 

The idea behind BUSCO is that we should be able to find almost all single copy orthologs in a complete genome assembly because these chosen genes should only have one copy per genome and are only very occasionally duplicated within a genome.
As these genes are well conserved across evolution for the particular taxa, missing or fragmented genes indicate that there is a problem with the assembly.

Run BUSCO on each of your three assemblies. 
Give each assembly a different descriptive name (eg. canu, canu_adjust, and flye) but use the same --out_path  of `3_busco`.
This will put all of our results into the same directory and keep things neat and tidy. 

```
busco -i 1_canu_assembly/c_hepaticus.contigs.fasta -o canu -l bacteria_odb10 -c 2 --out_path 3_busco --mode geno
```

You should get output printing to the terminal that includes something similar to the following:

```
         ***** Results: ****                                                                          
         C:89.0%[S:84.3%,D:0.0%],F:2.1%,M:8.8%,n:124                   
         4359    Complete BUSCOs (C)                                     
         4125    Complete and single-copy BUSCOs (S)                     
         234     Complete and duplicated BUSCOs (D)                      
         104     Fragmented BUSCOs (F)                                   
         433     Missing BUSCOs (M)                                      
         4896    Total BUSCO groups searched 

```

These results are also stored in the `short_summary*.txt` file in `3_busco/assembly_name/` directory. 

* *Which assembly is the best according to the BUSCO metrics? Why?*
* *C. hepaticus has a haploid genome so it is unsurprising that there are no duplicated BUSCOs in our assemblies. What would the presence of a large number of duplicated BUSCOs suggest if we were looking at a diploid assembly? What about if a diploid assembly had no duplicated BUSCOs?*

## 4. QUAST with a reference

From looking at the QUAST and BUSCO results above, it probably appears as though the Flye assembly is better than the Canu assemblies. 
However, let's look at the command I used to generate the Flye assembly. 

```
flye --pacbio-corr 0_data/c_hepaticus.fastq --genome-size 1500000 --threads 2
```

The `--pacbio-corr` option means that the input reads are corrected PacBio reads but this isn't actually the case. 
The input reads used for this assembly are the same as for the Canu assemblies we generated in the first part of the practical - raw PacBio CLRs.

So what does that mean? 

Is the Flye assembly actually the best of our three assemblies or are our methods for measuring assembly quality insufficient?

Let's compare our three assemblies with a reference assembly of C. hepaticus and see if this helps us to decide. 
Note that this reference is not for the same sample as our raw PacBio reads came from. 

```
quast -r 0_data/reference.fasta -o 4_quast_ref -t 2 --labels "canu,canu_adjusted,flye" assembly1.fa assembly2.fa assembly3.fa
```

Again examine the `icarus.html` file as well as the QUAST documentation to answer the following:

* *QUAST has compared all three of the provided assemblies with the reference genome and gives us some new information. Do you still agree with your earlier assessment of which assembly is the best? Why or why not?*
* *Keeping in mind that we don't always have a reference genome available, what do these results say about the importance of knowing how a dataset or assembly was generated?* 
* *If we were assessing the quality of a single de novo assembly, why might it not make sense to do a QUAST comparison with the species reference?*

## 5. Merqury

In addition to BUSCO and QUAST, there are other methods and tools for measuring assembly quality.

Merqury is one such tool. 

* *Do some research to find out how Merqury evaluates assembly quality.*
