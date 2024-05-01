# BIOINF 3010 & 7150

## Graph Genomes Practical: Part 1

##### _By Chelsea Matthews - Adapted from Yassine Souilmi_

# 1. Background

Many bioinformatics analyses begin by aligning reads to a linear reference genome and you've done this in earlier practicals. 
Reads are aligned to the reference based on the similarity of their sequence to the reference sequence.
This means that if a newly sequenced genome has a genomic sequence that is very different to the reference genome for that species, reads that originate from these different regions may not align to the reference genome or might align to the wrong region of the reference genome.
When a read doesn't align at all, it is often discarded and so we lose some information about the newly sequenced sample genome.
Similarly, when a read incorrectly aligns, we falsely interpret what the newly sequenced genome looks like.
This bias in an analysis is called _reference bias_ and it can have a significant impact on our findings. 

Pan-genomes, and in this case, pan-genome graphs help to reduce reference bias. 
Instead of representing one possible genome for that species (eg. the E. coli reference you used in an earlier practical) genome graphs represent multiple possibilities for that species genome.
This means that when we align reads to the genome graph, reads are more likely to be similar to a sequence represented within the graph and so alignment rates and alignment accuracy are better.

!["Alignment to a reference vs alignment to a graph"](../../images/ref_vs_graph_figure.png)

**Note:**
Bioinformatics methods where we align reads to a linear reference genome (like the E. coli reference you used in an earlier practical) are cheaper and simpler than building and using a genome graph and these linear reference-based methods are very well documented and tested. 
The methods and software for building and using genome graphs are still under development.
For these reasons, linear reference genomes are currently far more commonly used than graph-based approaches but this may change in the near future.

## Part 1: Variant Graphs

First you will need to load the required tools.

```
source activate bioinf
```

This includes all of the tools we'll be using today.

- vg - for constructing genome graphs

- graphviz - graph visualisation software

- Bandage - graph visualisation software

The following commands should each display some usage instructions for each tool

```
vg
Bandage --help
dot --help #dot is part of the graphviz package
```

### Get the data

Make directories and copy the data over:

```
cd ~
mkdir -p GraphGenomes/{variant,msa,mygraph}
cp /shared/data/Graph_Genomes/prac1/tiny* ~/GraphGenomes/variant/.
cp /shared/data/Graph_Genomes/prac1/HLA_haplotypes.fa ~/GraphGenomes/msa/.
cp /shared/data/Graph_Genomes/prac1/{mystery.fq,cannabis.fasta} ~/GraphGenomes/mygraph/.
cd ~/GraphGenomes/variant
```

You should have the following two files in `~/GraphGenomes/variant`.

- tiny.fa
- tiny.vcf.gz

### Have a look at the data

We have two files. The first, `tiny.fa`, is the reference sequence. 
Take a look at it. 

```
cat tiny.fa
```

- What is the name of the sequence? 
- How long is it? 

The second file is a VCF. 

VCF stands for "Variant Call Format".
VCF files describe genomic variants including:
- the variant location within the reference genome
- the genomic sequence of the reference allele at that location
- the genomic sequence of the variant

Small variants like SNPs, insertions, and deletions are commonly documented in VCFs but larger structural variants like duplications and inversions may also be included. 
These variants can come from one sample or from many members of a population.

Let's have a look.

```
gunzip tiny.vcf.gz
cat tiny.vcf
```

- How many variants are in this VCF? What do you think this means for our graph?

- What do the first 5 columns of the VCF contain?

That's all we need to build our pan-genome! 

### Build a graph with just the reference sequence

To build the graph, use the following command:

```
vg construct -r tiny.fa -m 32 > ref.vg
```

- What does the `-m` option do in `vg construct`? 

We can have a look at the text representation of the graph by inspecting the GFA output produced by vg view. 
GFA stands for Graphical Fragment Assembly. 

It has a header line that begins with `H`.
Lines beginning with `S` are segments (or **nodes**) of the graph.
Lines beginning with `L` are links (or **edges**) and lines beginning with `P` are **paths** through the graph.
A path through a graph represents an observed biological sequence or haplotype and a single graph can have many paths. 

```
vg view ref.vg
```

#### Questions:

- How many nodes does this graph have?
- How many edges does this graph have?
- How many paths does this graph have?

While GFA format is human readable, it's more intuitive to look at graphs with nodes and edges displayed visually.

```
vg view -dp ref.vg | dot -Tpdf -o ref_dot.pdf
```

It's a graph, but not a very exciting or useful one. 
Before we continue looking at visualisation options, let's build a graph with some variation in it. 

### Build a graph that's actually useful

Build a variation graph using the command below and then have a look at it in GFA format.

```
vg construct -r tiny.fa -v tiny.vcf -m 32 > tiny.vg
vg view tiny.vg > tiny.gfa
```

Have a look at the `tiny.gfa` file. 

- How many nodes, edges, and paths does this graph have?

## Part 2: Visualisation

We have a number of options when it comes to visualisation. 

#### Option 1 - dot from the Graphviz package

!["Simple graph"](../../images/tiny_ref_dot.png)

```
vg view -dp tiny.vg | dot -Tpdf -o tiny_dot.pdf 
```

- What do the `-d`, `-p`, `-S`, and `-n` options do in the `vg view` command? 

### Option 2 - Bandage

!["Bandage graph"](../../images/tiny_bandage.png)

For larger graphs, `Bandage` is sometimes a better solution.
It was designed for visualising assembly graphs (these are the graph structures that are made by assembly tools during *de novo* genome assembly).
Bandage can't display path information but it does allow us to visualise the structure of graphs that are many Megabases in size.

Bandage requires the graph to be in .gfa format.

```
vg view tiny.vg > tiny.gfa
Bandage image tiny.gfa tiny_bandage.png
```

#### Option 3 - vg viz

!["tiny.svg"](../../images/tiny.svg.png)

`vg viz`￼￼ produces a linear layout of a graph that scales well between small and larger graphs without losing path information.
Unfortunately, it doesn't seem possible to view .svg files from within our R-studio machines.
To get around this, you will need to download the .svg fle to your local machine.

`vg viz` requires the .xg index of the graph you want to visualise. 

```
vg index -x tiny.xg tiny.vg

vg viz -x tiny.xg -o tiny_viz.svg
```

## Part 3: Align reads to the graph

Now we've built a graph and looked at it but we haven't done anything with it yet. 

It would be useful if we could align reads to the graphs as this is often a first step in a bioinformatics analysis.

Mapping reads to a graph is done in two stages: first, seed hits are identified and then a sequence-to-graph alignment is performed for each individual read.
Seed finding allows vg to spot candidate regions in the graph to which a given read can potentially map to.

To do this, vg requires two indexes, an XG index and a GCSA index.
We already have the xg index so we'll make the gcsa index next. 

```
vg index -g tiny.gcsa -k 16 tiny.vg
```

We don't have any reads for our data so we'll use vg to simulate some using ￼￼`vg sim￼￼` as below. 

```
cd ~/GraphGenomes/variant
mkdir reads

# generate three reads and then have a look at them
vg sim -l 20 -n 3 -e 0.05 -i 0.02 -x tiny.xg tiny.vg -a | vg view -X - > reads/sim_reads.fq
```

- What do the parameters in the `vg sim` command do? 

Now we will map these three reads to our graph. 
When it comes to mapping in `vg`, we have three options. 

`vg giraffe` is designed to be fast for highly accurate short reads, against graphs with haplotype information.

`vg map` is a general-purpose read mapper.

`vg mpmap` does "multi-path" mapping, to allow describing local alignment uncertainty. This is very useful for transcriptomics.

We will be using `vg map`.


```
vg map --fastq reads/sim_reads.fq -x tiny.xg -g tiny.gcsa > tiny.gam
```

Let's re-visualise to see the alignments. 
Notice that the ￼￼`vg view`￼￼ command below uses the ￼￼`-S￼￼` option to simplify the dot output.

```
vg view -dS tiny.vg -A tiny.gam | dot -Tpdf -o alignment.pdf
```

- Why does your graph look different to mine and/or your neighbours?

- What do you think the different colours of the nodes indicate?

Feel free to try it without the ￼￼`-S￼￼` option as well as any other parameter combos you like to see the difference. 

Looking at the reads of a graph can be quite overwhelming when you've got a lot of alignments and so it can sometimes help your visualisation to add these alignments back into the graph as paths, in the same way as we have the reference path.

```
vg augment --label-paths tiny.vg tiny.gam > augmented.vg
```

Let's take a look at the graph with the dot output format as well as with `vg viz`

```
vg view -dp augmented.vg | dot -Tpdf -o augmented_dot.pdf 

vg index -x augmented.xg augmented.vg
vg viz -x augmented.xg -o augmented_viz.svg
```

Note that this path information doesn't include the quality of the mapping, it just tells you where the best alignment was. 

Next practical we'll look at doing more with our read alignments (calling variants, embedding variants back into the graph). 


## Part 4:  Multiple Sequence Alignment graph

As well as the method we used above (a reference genome and a vcf), we can also build graphs using Multiple Sequence Alignment. 
This method is often used when we have a number of high quality haplotype-resolved assemblies available.
We align them with each other to see where each of the sequences is similar and where they differ. 

```
cd ~/GraphGenomes/msa
```

Take a look at the `HLA_haplotypes.fa` file. 

- What does it contain? 

- How many sequences are there? 

Now build the graph with `vg msga` and index it. 
It's going to print out a heap of warnings because the `vg msga` module is not strictly the recommended way to produce multiple sequence alignments. 
We don't have the software for other methods and this works for our purposes. 

```
vg msga -f HLA_haplotypes.fa -t 2 -k 16 | vg mod -U 10 - | vg mod -c -X 256 - > hla.vg
vg index hla.vg -x hla.xg -g hla.gcsa
```

- What do the parameters `-U` and `-c` from `vg mod` do?

And visualise again! 

```
vg view hla.vg > hla.gfa

Bandage image hla.gfa hla.png
vg view -dp hla.vg | dot -Tpdf -o hla.pdf
vg viz -x hla.xg -o hla.svg
```

Are all of these visualisation methods a good/convenient way to look at this graph? 

## Your Turn! 

Build a multiple sequence alignment graph from the cannabis.fasta dataset in the `~/GraphGenomes/mygraph` directory. 
Consider the first sequence in the fasta file the reference sequence.

Try to do/answer the following: 
 
- Visualise your graph however you want.
- What is the major difference between each of the three samples compared with the `pink_pepper` reference. There should be one sequence with a duplication, one with a deletion, and one with an insertion. There are also other small differences between the sequences but these are not the focus. 
- Change the colour of your Bandage plot (this is entirely for fun, it's not essential)
- Align the supplied reads to your graph and try to work out which sample they most likely came from. 




