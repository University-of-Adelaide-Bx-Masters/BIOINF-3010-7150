# Week 7 Practical: **Genome Assembly** 
*By Zhipeng Qu and Chelsea Matthews* 

Today we will be looking at *de novo* assembly with more of a focus on larger, more complex genomes.  

# Get Started!

As with previous weeks' practicals, you will be using RStudio to interact with your VM.

While the focus of this section of the coursework is on the assembly of large complex genomes, it's not feasible for us to actually assemble a large complex genome (i.e. human genome) because, frankly, our VM's are too small and sitting around for hours waiting for an assembly to finish is boring (You can check the benchmark for the computational cost using Flye from [here](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-benchmarks)). 
So, we'll be assembling one of the small eukaryotic genomes, fission yeast (*Schizosaccharomyces pombe*), using Long Reads (LR) from nanopore and pacbio sequencing.   
Even though the fission yeast genome is fairly small (~15 Mb with 3 chromosomes), but it's still more complicated than the prokaryotic organisms. It will still take ~20 minutes for the assembly to run even we just use 5x genome coverage. 

## Dataset

Due to limitation of computating resources in our VM, we will use subsets of raw sequencing reads. The original sequencing dataset is very big, you can access it from this [link](https://www.ncbi.nlm.nih.gov/sra?term=SRP352919). Here are some useful genome informations about fission yeast:

- Reference genome: ASM294v2
- Number of chromosomes: 3 nucleus chromosomes
- Genome size (reference): 12,591,251 bp
- Ploidy: Haploid

## Tools/packages

Tools/packages used in this Prac:

| Too/Package    | Version      | URL                                                        |
|----------------|--------------|------------------------------------------------------------|
| fastQC         | v0.11.9      | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| assembly-stats | v1.0.1       | https://github.com/sanger-pathogens/assembly-stats         |
| jellyfish      | v2.2.10      | https://github.com/gmarcais/Jellyfish                      |
| genomescope    | v1           | https://github.com/schatzlab/genomescope                   |
| flye           | v2.8.1-b1676 | https://github.com/fenderglass/Flye                        |
| QUAST          | V5.2.0       | https://github.com/ablab/quast                             |
| BUSCO          | v5.4.4       | https://busco.ezlab.org/                                   |

## Running time estimation (based on teaching VM)

| Step              | Tool/Package        | Estimated run time  |
|-------------------|---------------------|---------------------|
| QC                | fastqc              | < 1 min             |
| QC                | assembly-stat       | < 1 min             |
| genome survey     | jellyfish           | < 5 mins            |
| genome survey     | genomescope         | < 1 min             |
| genome assembly   | flye + nanopore_5x  | ~20 mins            |
| genome assembly   | flye + pacbio_5x    | ~20 mins (optional) |
| genome assembly   | flye + nanopore_10x | ~30 mins (optional) |
| genome assembly   | flye + pacbio_10x   | ~30 mins (optional) |
| genome assessment | QUAST               | < 5 min             |
| genome assessment | BUSCO               | ~10-20 mins (each)  |

The following table shows the estimated run time on VM for the different processes:

## What you will learn in this Prac

- Practice bash commands you have learned
- Practice QC for illumina data you have learned
- Learn how to do de novo genome assembly using Flye with long reads
- Learn how to assess your genome assembly quality using QUAST and BUSCO

# Let's start!

There are 4 different major steps (5 parts actually) in this Prac.  Please follow the instructions __in order__, because some commands will rely on the results from previous commands. Feel free to talk to tutors/instructors if you have a problem/question. 

## Part 1. Set up and project preparation

### 1.1 Prepare folder structure

Planning your directory structure from the beginning makes your work easier to follow for both yourself and others. 

I put initial input data into a `data` folder, and output files from different processing stages into separate folders in a `results` folder. If there are databases involved, I also create a `DB` folder. I store all scripts in a separate `scripts` folder (we won't use this folder in this Prac). I also use a numbered prefix such as `01_raw_data` to label different folders. All of these naming rules are just personal preference, and feel free to build your own project folder structure rules, and keep it consistent for your different projects in future.

The following is the folder structure for this genome assembly project:

```
./prac_genome_assembly/
├── 01_bin
├── 02_DB
├── 03_raw_data
├── 04_results
│   ├── 01_QC
│   ├── 02_genome_survey
│   ├── 03_genome_assembly
│   └── 04_genome_assessment
└── 05_scripts
```

We can use following commands to build this folder structure:

```bash
cd ~/
mkdir prac_genome_assembly
cd prac_genome_assembly
mkdir 01_bin 02_DB 03_raw_data 04_results 05_scripts
cd 04_results
mkdir 01_QC 02_genome_survey 03_genome_assembly 04_genome_assessment
```

If you want to check your folder structure:
```bash
cd ~/
tree ./prac_genome_assembly
```

### 1.2 Input data and additional tools
We will be using raw sequencing data from different sequencing platforms in this Prac. These fastq files can be found in `~/data/prac_genome_assembly/01_raw_data` folder, and are including:

| File(s)                                    | Platform | Coverage | Description                                                |
|--------------------------------------------|----------|----------|------------------------------------------------------------|
| illumina_SR_20x_1.fq, illumina_SR_20x_2.fq | Illumina | ~20x     | Paried-end (PE150) short reads from Illumina MGISEQ-2000RS |
| nanopore_LR_5x.fq                          | Nanopore | ~5x      | Long reads from Nanopore PromethION                        |
| nanopore_LR_10x.fq                         | Nanopore | ~10x     | Long reads from Nanopore PromethION                        |
| pacbio_LR_5x.fq                            | PacBio   | ~5x      | Long reads from PacBio_SMRT Sequel                         |
| pacbio_LR_10x.fq                           | PacBio   | ~10x     | Long reads from PacBio_SMRT Sequel                         |

The original dataset can be found [here](https://www.ncbi.nlm.nih.gov/sra?term=SRP352919).

We will use the paired-end illumina short reads to do genome survey analysis (estimate the genome size), and use the other four Long Reads (LR) files to do genome assembly separately.

The reference genome is also provided, and you can find it in folder `~/data/prac_genome_assembly/02_DB`.

Next we can copy the corresponding input files into our corresponding project folders:

```bash
cd ~/prac_genome_assembly/02_DB
cp ~/data/prac_genome_assembly/02_DB/* ./
cd ~/prac_genome_assembly/03_raw_data
cp ~/data/prac_genome_assembly/03_raw_data/*.fq ./
```
Because `assembly-stats` and `genomescope` are not globally installed in VM, we need to maually put them somewhere we know. You can find these two scripts in folder `~/data/prac_genome_assembly/01_bin/`, and we will put them in the `01_bin` folder of our project.

```bash
cd ~/prac_genome_assembly/01_bin
cp ~/data/prac_genome_assembly/01_bin/* ./
```

Now all the setup work is done. Let's move to part 2.

## Part 2, QC on input data

In this part, we will be using `fastQC` to check the sequencing quality of illumina reads, and using `assembly-stats` to acquire some statistics about the long reads. 

### 2.1 QC for illumina reads

Let's have a look at the short reads (illumina reads) using `fastQC` first

Before we start the actual analysis, we need to activate the `conda` environment so that we can use packages/tools that are only installed in specific environments (It's always a good idea to check what `conda` environment you are in before you do following analyses). To activate `bioinf` environment, you can run:

```bash
source activate bioinf
```

After you `activate` this `bioinf` environment, your terminal prompt should be changed to `(bioinf) axxxxxxx@ip-xx-xxx-x-xx`, which indicates that now you are in the `bioinf` environment.

Then we can run `fastqc` to check the quality of the short reads:

```bash
cd ~/prac_genome_assembly/04_results/01_QC
fastqc ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_1.fq ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_2.fq -o ./ -t 2
```

* *How many sequences are there in the dataset?*
* *How good the illumina reads are?* 

### 2.2 Long reads
We will be using `assembly-stats` to get some statistics for the long reads.

```bash
cd ~/prac_genome_assembly/04_results/01_QC
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/nanopore_LR_5x.fq
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/pacbio_LR_5x.fq
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/nanopore_LR_10x.fq
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/pacbio_LR_10x.fq
```

* *Which dataset has the largest average read length?*
* *Which dataset has the longest individual read?*
* *Assuming that our genome is approximately 15 Mb long, what coverage do these reads give us? Remember that Coverage = (total number of bases in reads)/genome size*
* *Do you think this is sufficient to generate a good assembly? Why or why not?*

## Part 3, genome survery analysis (genome size estimation)
When we start a whole genome sequencing project for a new species, we normally need to collect some genomics information before we do the de novo genome assembly, such as we need to know how big the genome is. In the lecture, we had learned that we can use lab-based flow cytometry to estimate the genome size, and we can also use short reads (illumina reads) to do this computationally. In this part, we will learn how to estimate the genome size using short reads.

### 3.1 get k-mer distribution
The first step of genome size estimation is to get the k-mer distribution using the available short reads. We can use `jellyfish` to do this.

```bash
cd ~/prac_genome_assembly/04_results/02_genome_survey

jellyfish count -C -m 21 -s 4G -o illumina_SR_20x.21mer_out \
    ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_1.fq \
    ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_2.fq

jellyfish histo -o illumina_SR_20x.21mer_out.histo illumina_SR_20x.21mer_out

```

`jellyfish` will break short reads into fixed length short sequences, which we call them [k-mers](https://en.wikipedia.org/wiki/K-mer) (we use 21-mer in this project). The size of k-mers should be large enough allowing the k-mer to map uniquely to the genome (a concept used in designing primer/oligo length for PCR). However, too large k-mers leads to overuse of computational resources. `21` is normally a good start. In the `jellyfish count` command, `-C` means we count k-mers at both strands, `-s 4G` is used to control the memory usage, and `-m 21` means we will count 21-mers. After we count k-mers, we use `jellyfish histo` to get the frequency of k-mers with certain copy numbers.

### 3.2 genome survey analysis using genomescope
After we get the `histo` file from the `jellyfish` run, we can use that to do genome survey with `genomescope`. `genomescope` is a R script, we can run it with following command:

```bash
cd ~/prac_genome_assembly/04_results/02_genome_survey

Rscript ~/prac_genome_assembly/01_bin/genomescope.R illumina_SR_20x.21mer_out.histo 21 150 illumina_SR_20x.21mer
```

In the command, `21` means we are using 21-mer, `150` is the short reads length, and `illumina_SR_20x.21mer` will be the output folder. We can check the k-mer distribution by checking the file `plot.png`, which is normally located in the output folder `illumina_SR_20x.21mer`


## Part 4, genome assembly 

### 4.1 assembly using Flye
Okey! Now we can do the de novo genome assembly. There are tons of assemblers to do de novo genome assembly. Some popular ones are `canu`, `Flye`, `smartdenovo`, `wtdbg`, etc. In this prac, we will be using `Flye`.

Do a test run to see if it works using the following (make sure you are under the `bioinf` conda environment):

```
flye
```

The Flye usage instructions should be printed to the screen. 

Now let's using the flye to do assembly on the `nanopore_LR_5x` dataset. 

```bash
cd ~/prac_genome_assembly/04_results/03_genome_assembly

flye --nano-raw ~/prac_genome_assembly/03_raw_data/nanopore_LR_5x.fq --out-dir nanopore_LR_5x --threads 2
```

It should take ~20 minutes for the assembly to complete.   

### 4.2. While we wait ...  Resources for genome assembly

Larger, more complex genomes contain high percentages of repeats. Assembling these repeats accurately with short reads is not effective. Long reads that are able to span these repeat regions result in much better assemblies. There are many long read assembly tools available including the following: 

* Canu 
* Flye
* NextDenovo
* Wtdbg2 
* Raven
* Falcon
* Shasta
* miniasm

These tools all do more or less the same job (ie. long read assembly) but with slight differences. 

When we are assembling a small genome, memory (RAM) and CPU (compute threads/cores) requirements are relatively low and so we can choose our assembler based almost entirely on the quality of the resulting assembly. 
While we can still do this to some degree for larger genomes, resource allocations and limitations in High Performance Computing (HPC) environments mean that some tools may not be suitable for our particular dataset (we may not be able to run the assembly to completion due to time or resource limitations) and we need to more carefully consider the parameters we use for the assembly so that we don't need to re-run it unnecessarily. 
Even without considering the cost of wasted resources (and this can be measured in dollars where we are using a service like Amazon Web Service (AWS)), re-running these assemblies can take a long time. 

One of the difficulties with assembling larger genomes is working out what resources a tool will require and deciding whether it will be suited to assembling your genome of interest within the resources available to you. 
The authors of Flye did some benchmarks on computational resources required when assembling genomes from different species with different input data. You can check this from this [link](https://github.com/fenderglass/Flye#flye-benchmarks). Check following questions after you have a look at the table.

* *How many CPU hours (approximately) are required to assemble a bacteria?*
* *How many CPU hours (approximately) are required to assemble a mammal with high raw-read coverage?*
* *How many CPU hours (approximately) are required to assemble a mammal from HiFi reads?*
* *Why is the number of CPU hours and memory requirement lower for a HiFi assembly compared with a non-HiFi read assembly?*
* *If you had access to a single node with 48 compute threads and sufficient memory, how long (in hours) would it take to assemble a "well-behaved large genome"?* 
* *How long would this same assembly take on your local VM (assuming memory was not a limitation) with two compute threads?*

You can see now why we aren't assembling a larger genome (i.e. human genome) using our VM. 

### 4.3 Looking at the assembly report

Flye will generate 5 folders, which store output files from 5 stages of Flye, and some important individual files. These includes:

- params.json: the parameters that Flye used for this run
- assembly_graph.{gv or gfa}: Final repeat graph. The edges of repeat graph represent genomic sequence, and nodes define the junctions. You can get more info about the repeat graph from [here](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-repeat-graph) 
- **assembly.fasta**: This is the final assembly. Contains contigs and possibly scaffolds.
- assembly_info.txt: Extra information about contigs.
- flye.log: Log report showing all running info and summary of final assembly.

You can check the log report using following commands:

```
cd prac_genome_assembly/04_results/03_genome_assembly/nanopore_LR_5x
ls 

less flye.log
```
Type `G` to go to the end of the file, and you will see a brief summary about the assembly.

* *How many contigs does the final assembly have?*
* *How long is the final assembly and how does this compare with the estimated size of the genome?*

When you are finished, exit out of the report with `q`. 

### 4.4 assembly using other datasets

You can do the genome assembly for the other three LR datasets:

```bash
cd ~/prac_genome_assembly/04_results/03_genome_assembly

flye --pacbio-raw ~/prac_genome_assembly/03_raw_data/pacbio_LR_5x.fq --out-dir pacbio_LR_5x --threads 2
flye --nano-raw ~/prac_genome_assembly/03_raw_data/nanopore_LR_10x.fq --out-dir nanopore_LR_10x --threads 2
flye --pacbio-raw ~/prac_genome_assembly/03_raw_data/pacbio_LR_10x.fq --out-dir pacbio_LR_10x --threads 2

```

Each of these will take 20-30 mins. I encourage you to run these jobs before the next session, and go through the reports (flye.log) to get ideas about the different assemblies, and we will look into them in more details in the next session.

### 4.5 Additional information - diploid assembly

The fission yeast is a small haploid genome which means that it is fairly straightforward nowadays to assemble. 
What is more challenging are larger diploid and polyploid genomes.

* *Can you think of a few reasons why diploid or polyploid genomes are more difficult to assemble than haploid genomes? Relate these reasons back to the OLC algorithm if you can.*
* *If we were assembling a diploid genome with a genome size of 2Gbp, under what circumstances might our resulting assembly be greater than 2Gbp in length?*

Canu is capable of assembling both haplotypes of a diploid genome - termed phasing. 
Have a look at the paper below (mainly the first figure) to understand how this works.

[Haplotype-Resolved Cattle Genomes Provide Insights Into Structural Variation and Adaptation](https://www.biorxiv.org/content/10.1101/720797v3.full)

