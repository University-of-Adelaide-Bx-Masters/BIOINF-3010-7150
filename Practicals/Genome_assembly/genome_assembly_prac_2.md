# Week 8 Practical session 1: Assembly Quality Assessment
## By Chelsea Matthews

This practical will continue on from Week 6 Practical 1. 

## What are we doing today?

In the first practical we assembled the C. hepaticus genome twice using Canu.
Today we will look at measuring *de novo* assembly quality with a focus on reference free quality metrics. 
In particular we will look at contig length distributions and BUSCO. 

## 1. Setting up

Along with the two Canu assemblies you generated on Tuesday, you will also need to copy over three other assemblies from the `/data` directory. 
These are:

* Flye assembly with 16x PacBio CLRs - the same data as we used in Canu yesterday
* Canu assembly with ~60x PacBio CLRs
* a C. hepaticus reference genome

```
mkdir Project_6/1_flye_assembly Project_6/1_canu60x_assembly 
cp data/project_6/assembly.fasta Project_6/1_flye_assembly
cp data/project_6/canu60x.fasta Project_6/1_canu60x_assembly
cp data/project_6/reference.fasta Project_6/0_data/
```

We will now make a new conda environment that will include the tools BUSCO and QUAST. 

```
conda create -n quality -c conda-forge -c bioconda busco=5.3.1 quast=5.0.2
conda activate quality
```

Your `Project_6` directory should now contain four assembly directories, each containing a different assembly. 
These directories should be named as follows. 

```
1_canu60x_assembly
1_canu_assembly
1_canu_adjusted_assembly
1_flye_assembly
```

If they're not, re-name them to match using `mv oldname newname` command. 

Now let's compare these assemblies. 

## 2. QUAST


The tool QUAST (QUality ASsesment Tool) (documentation here)[http://quast.sourceforge.net/docs/manual.html] can be run on one or more assemblies at once and produces some handy comparison statistics and graphics.
It can be run with or without a reference genome. 
If we provide a reference genome, each of the assemblies will be compared to this reference and QUAST will produce some additional statistics. 

Here we will run it first without a reference. 

```
quast -o 2_quast -t 2 --labels "canu, canu_adjusted, canu60x, flye" \
	1_canu_assembly/c_hepaticus.contigs.fasta \
	1_canu_adjusted_assembly/c_hepaticus.contigs.fasta \
	1_canu60x_assembly/canu60x.fasta \
	1_flye_assembly/assembly.fasta
```

This won't take very long.

Open the `icarus.html` file located in the `2_quast` directory in a web browser. 

Click on the QUAST Report link and have a look at the metrics produced by QUAST for each of the four assemblies. 

* *Based on the N50 and N75 scores, which is the best assembly?*

You've previously learnt about N50 as a measure of assembly contiguity ([refresh your memory here](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50)) but N50 on its own only gives us a snapshot into assembly contiguity.
A better way to look at assembly contiguity is to inspect a cumulative contig length plot.
This is particularly helpful when we want to compare multiple different assemblies of the same genome or species. 

Scroll down and inspect the Cumulative Length plot. Look at the Nx plot and GC Content plots too while you're there. 

* *Given that the C. hepaticus genome should be about 1.5Mbp long, which assembly do you think looks the best based on the contig length distributions?*
* *Looking just at the canu and canu_adjusted assemblies, do you think it was worth adjusting the correctedErrorRate?*
* *Spend some time interpreting the cumulative length plot. Under what circumstances would a cumulative length plot be more informative than just comparing N50s?*
* *What would an assembly with only one contig look like on this plot?*

The image below is a slightly different version of a cumulative contig length plot. 
The assemblies represented are publicly available cannabis assemblies. 
The dashed horizontal lines are the approximate sizes of the male and female genomes with the male being the longer.
Inspect the plot and answer the questions below. 

![Cumulative contig length plot](images/cum_contig_length_dist.png)

* *What benefit is there to having contig length on the x axis instead of the contig number in size order as in the QUAST plot?*
* *One of these assemblies was generated from short illumina reads. Which do you think it is?*
* *Can you tell which are contig level assemblies and which are scaffolded?*

Now let's have a look at the Icarus contig size viewer. 

Click on the "View in Icarus contig browser" link at the top of the page. 

This isn't really new information but it can be nice to visually see the lengths of our contigs. We can see that the Flye assembly has many small contigs compared with the other three assemblies. 
It's also a nice way to see the N50 and N75 values distributed on a visual representation of the assembly contigs. 

* *Based on what you've seen in the QUAST report, which assembly do you think looks the best?*

## 3. BUSCO (Benchmarking Universal Single Copy Orthologs)

BUSCO is a tool that helps us measure how well we have assembled the gene space of an assembly. 
It does this by searching our assembly for a list of genes that should be present in our assembly in a single copy. 
These genes are known as "Universal Single Copy Orthologs" and are genes that are present in at least 90% of the species/clade members and are only present in a single copy within at least 90% of the species/clade.
It is essentially a list of genes that we are almost certain should be present within our genome once.

Obviously this list of genes is different for different species/clades and so BUSCO has a list of datasets that we can choose from. 

Have a look at this list using the command below. 

```
busco --list-datasets
```

We have assembled Campylobacter hepaticus and so will be using the `campylobacterales_odb10` database. 
We could use the bacteria database, because C.hepaticus is a bacteria, but the campylobacterales database will contain single copy orthologs that are more specific to the Campylobacter genus. 

* *Do you think that the campylobacterales database will have more or less single copy orthologs than the generic bacteria database?*

You can see the number of universal single copy orthologs in each BUSCO database at [this](https://busco.ezlab.org/list_of_lineages.html) website. 

Check to see if you were right. 

Let's run BUSCO on our assemblies.

```
# run BUSCO on the first Canu assembly
busco -i 1_canu_assembly/c_hepaticus.contigs.fasta -o canu \
	-l campylobacterales_odb10 -c 2 --out_path 3_busco --mode geno
```

You should get output printing to the terminal that includes something similar to the following. 

Don't worry if the numbers are different! The results below are for a cannabis assembly. 

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

These results are also stored in the `short_summary*.txt` file in `3_busco/canu/` directory. 

Run BUSCO on your remaining assemblies and the c_hepaticus reference genome as follows:

```
busco -i 1_canu_adjusted_assembly/c_hepaticus.contigs.fasta -o canu_adjust \
	-l campylobacterales_odb10 -c 2 \
	--out_path 3_busco --mode geno
busco -i 1_flye_assembly/assembly.fasta -o flye \
	-l campylobacterales_odb10 -c 2 \
	--out_path 3_busco --mode geno
busco -i 0_data/reference.fasta -o reference \
	-l campylobacterales_odb10 -c 2 \
	--out_path 3_busco --mode geno
busco -i 1_canu60x_assembly/canu60x.fasta -o canu60x \
	-l campylobacterales_odb10 -c 2 \
	--out_path 3_busco --mode geno
```

Now let's visualise the results. 

We will use the `generate_plot.py` script that is usually installed with BUSCO. 
We haven't got it because we installed BUSCO using conda. 
It is stored in a gitlab repository [here](https://gitlab.com/ezlab/busco/-/blob/master/scripts/generate_plot.py) but you can get the file by running the command below.

```
mkdir scripts
cd scripts
wget https://gitlab.com/ezlab/busco/-/raw/master/scripts/generate_plot.py?inline=false -O generate_plot.py
cd ..
```

This script requires that all of the BUSCO short summaries be placed into a single folder. 

```bash
cd 3_busco
mkdir all_summaries

cp canu/short*.txt all_summaries
cp flye/short*.txt all_summaries
cp canu_adjust/short*.txt all_summaries
cp canu60x/short*.txt all_summaries
cp reference/short*.txt all_summaries
```

Do a quick check in the `all_summaries` folders to check that there are five files named like below:

`short_summary.specific.campylobacterales_odb10.assemblyName.txt`

We will also need to install the package `ggplot2` in R to create the plot. 

Move from the Terminal to the Console window in RStudio and run the commands below. 

```R
install.packages("rlang")
install.packages("ggplot2")
install.packages("grid")
```

Now, go back to the Terminal window and run the python script. The command below should do the trick. 

Note, you should do this from the `Project_6` directory. 

```
python3 scripts/generate_plot.py -wd 3_busco/all_summaries
``` 

If all goes well, this should produce a `busco_figure.png` in the `all_summaries` directory. 

Navigate to it and open the png. 

* *Which assembly is the best according to the BUSCO metrics? Why?*
* *C. hepaticus has a haploid genome but there are duplicated genes present in the canu60x assembly. What do you think this means?* 

* *Based on the BUSCO results, and excluding the reference assembly, which assembly do you think looks the best?*

## 4. QUAST with a reference

We don't always have a reference genome available for a species but in this case, we do. 

We will now run QUAST on our four assemblies but will provide the reference genome this time. 

```
quast -r 0_data/reference.fasta -o 4_quast_ref -t 2 \
	--labels "canu, canu_adjusted, canu60x, flye" \
	1_canu_assembly/c_hepaticus.contigs.fasta \
	1_canu_adjusted_assembly/c_hepaticus.contigs.fasta \
	1_canu60x_assembly/canu60x.fasta \
	1_flye_assembly/assembly.fasta
```

Again examine the `icarus.html` file in the `4_quast_ref` directory. 
Don't forget to click the "Extended report" link in the QUAST report. 
Have a good look through the new information and try to understand it. 

* *When interpreting the misassemblies and mismatches sections in the QUAST report, what do you think we might need to keep in mind?*
* *Looking at the Contig size viewer, do you think it is more appropriate to compare N50 or NG50? Why? Find out what NG50 is if you're not sure.*
* *If we were assessing the quality of a single de novo assembly, why might it not make sense to do a QUAST comparison with the species reference?*
* *QUAST has compared all four of the provided assemblies with the reference genome and gives us some new information. Do you still agree with your earlier assessment of which assembly is the best? Why or why not?*

None of these assemblies are particularly good. 16x of PacBio CLRs is not sufficient to make a good assembly.
However, even the Canu 60x assembly isn't anywhere near as good as the reference genome. 

* *What would be some strategies to improve one of these assemblies or, to make a better assembly if you were starting from scratch?*

