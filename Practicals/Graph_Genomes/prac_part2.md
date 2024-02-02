
## Graph Genomes Practical: Part 2

## Chelsea Matthews - modified from Yassine Souilmi

### Learning Outcomes

1. See how alignment parameters (in this case minimum read length match) can alter alignment rates.

2. See how read alignment rates to a linear reference compare with alignment rates to a graph incorporating variants found within that population.

3. Have a basic understanding of genotyping.

### Introduction

One of the main goals for graph pan-genomes is to be able to use them as an alternative to a linear reference genome.
In order to do this, we need to be able to align reads to them.
We saw how this was done with a set of three reads in the first part of the practical.
Today, we'll be aligning reads to a much larger graph and comparing the read alignment rate of these reads with the read alignment rate of the same reads to a linear reference genome.⋅
This will hopefully show how read alignment rates are improved by aligning to a graph with some natural variation from the population built in.
This improvement is due to the way that different paths through the graph represent different possible genomic sequences and so a graph is more likely to contain a genomic sequence that is more similar to the newly sequenced sample than a linear reference genome.⋅

Note: The :muscle: emoji indicates that there's a task that is essential for the practical that I have not provided the code for.
You will need to use the help pages and the first practical to work out these snippets but it shouldn't be too difficult.
We will also work through these questions together in the practical session.

### Getting (re)started

Let's load the tools we'll need today. 

```
source activate bioinf
```

The tools we'll need today include:

```
vg
bcftools
jq -h
dot
```

You're already familiar with vg and dot from the first practical. 
You may have used `bcftools` before but if not, it's a tool that is designed for use in variant calling pipelines.
It allows us to manipulate VCF (variant call format) files which is what we'll use it for today. 
And finally, `jq`.
This is a tool for processing and filtering JSON files.
We will be using it, as well as a few other commands, to summarise read alignment rates against different graphs.

We'll create a new directory for todays prac so we don't get too confused by all the files. 

Copy all of the files from `data/graph_genomes/prac2/` into a new directory for todays practical.
`cd` into that directory.

```
mkdir graph_2
cp data/graph_genomes/prac2/* graph_2/.
cd graph_2
```

Check that you have the following files:

- z.fa
- z.fa.fai
- z.vcf.gz
- z.vcf.gz.tbi
- README.md
- yeti.fa
- yeti.vcf
- nylamo.reads

The README.md just gives a brief summary of the `z` data.

### Working with larger graphs

The z files contain 1Mbp of 1000 Genomes data for chr20:1000000-2000000.
As for the tiny example, let's' build one linear graph that only contains the reference sequence and one graph that additionally encodes the known sequence variation.
The reference sequence is contained in `z.fa`, and the variation is contained in `z.vcf.gz`.

Using the examples from the first prac and the `vg construct` help menu, do the following: 

:muscle: Build a reference-only graph named `ref.vg`.

:muscle: Build the same graph but with variation included named `z.vg`.

Default parameters are fine.
Look at the previous examples to figure out the command or use the `vg construct` help menu.

You might be tempted to visualize these graphs (and of course you are welcome to try), but they are sufficiently big already that your machine might run out of memory and crash.

In a nutshell, mapping reads to a graph is done in two stages: first, seed hits are identified and then a sequence-to-graph alignment is performed for each individual read.
Seed finding hence allows vg to spot candidate regions in the graph to which a given read can map potentially map to.
To this end, we need an index.
In fact, vg needs two different representations of a graph for read mapping XG (a succinct representation of the graph) and GCSA (a k-mer based index).
To create these representations, we use `vg index` as follows.

```
vg index -x z.xg z.vg
vg index -g z.gcsa -k 16 z.vg
```

Passing option `-k 16` tells vg to use a k-mer size of *16k*.
The best choice of *k* will depend on your graph and will lead to different trade-offs of sensitivity and runtime during read mapping.

As mentioned above, the whole graph is unwieldy to visualize.
But thanks to the XG representation, we can now quickly **find** individual pieces of the graph.
Let's extract the vicinity of the node with ID 2401 and create a PDF.

```
vg find -n 2401 -x z.xg -c 10 | vg view -dp - | dot -Tpdf -o 2401c10.pdf
```

2401c10.pdf![image](https://user-images.githubusercontent.com/1767457/116373399-74790080-a84c-11eb-8f77-95d57aa2beb2.png)

The option `-c 10` tells `vg find` to include a context of 10 nodes in either direction around node 2401.

Next, we want to map some reads to the graph, just like we did in the first practical. 
Again, we will use vg to simulate these reads.

```
vg sim -x z.xg -l 100 -n 1000 -e 0.01 -i 0.005 -a > z.sim
```

This generates 1000 (`-n`) reads of length (`-l`) with a substitution error rate of 1% (`-e`) and an indel error rate of 0.5% (`-i`).
Adding `-a` instructs `vg sim` to output the true alignment paths in GAM format rather than just the plain sequences.
Map can work on raw sequences (`-s` for a single sequence or `-r` for a text file with each sequence on a new line), FASTQ (`-f`), or FASTA (`-f` for two-line format and `-F` for a reference sequence where each sequence is over multiple lines).

We are now ready to map the simulated read to the graph.

For evaluation purposes, vg has the capability to compare the newly created read alignments to true paths of each reads used during simulation.

```
vg map -x z.xg -g z.gcsa -G z.sim --compare -j
```

This outputs the comparison between mapped and and true locations in JSON format but it's way too much for us to go through manually.

:question: have a look at just the first entry using the `head -n 1` command. Can you understand any of it?

We can process ths information with `jq`, `awk`, and `sed` to get a summary of alignment correctness.
We can use this quickly check if our alignment process is doing what we expect on the variation graph we're working on.
For instance, we could set alignment parameters that cause problems for our alignment and then observe this using the `--compare` feature of the mapper.

Let's summarise the alignment accuracy for the mapping we just did.

```
vg map -x z.xg -g z.gcsa -G z.sim --compare -j | jq .correct | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
```

:question: Can you work out what the `jq .correct` part is doing? 

:question: What is the `sed` part doing?

:question: What is the `awk` part doing? 

In contrast, if we were to set a very high minimum match length we would throw away a lot of the information we need to make good mappings, resulting in a low correctness metric:

```
vg map -k 51 -x z.xg -g z.gcsa -G z.sim --compare -j | jq .correct | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
```

It is essential to understand that our alignment process works against the graph which we have constructed.
This pattern allows us to quickly understand if the particular graph and configuration of the mapper produce sensible results at least given a simulated alignment set.
Note that the alignment comparison will break down if we simulate from different graphs, as it depends on the coordinate system of the given graph.

### Exploring the benefits of graphs for read mapping

Let's now use the reads that we simulated above to compare alignment rates between a linear reference and a graph.

:question: What do you think the 1000 reads we simulated above actually represent?

We've already built the graph that represents the linear reference (`ref.vg`).
We now need to build a graph with a subset of variants from the z.vcf to replicate what it would be like to align reads from a number of new samples to a pan-genome graph.

To do this, we will filter the `z.vcf.gz` file by allele frequency (a common way to select which variation we will include in a pan-genome). 
We will 
We will be building the graph containing only variants that are present in more than 20% of the population.
Variants that appear in less than 20% of the population will be excluded.

You can make a VCF with a minimum allele fequency with the command below noting that AF stands for allele frequency:

```
bcftools filter -i 'AF > 0.2' z.vcf.gz > 0.2_af_filtered.vcf
```

:question: How many variants were present in the original z.vcf.gz?

:question: How many variants are there in the 0.2_af_filtered.vcf?

:muscle: Build a new graph using `vg construct` from the reference sequence and the filtered vcf called filtered.vg.

:muscle: Create .xg and .gcsa indexes for `filtered.vg` and `ref.vg`.

We will now map the simulated reads to both of our graphs and compare the mapping identity.
Earlier when we compared mapping alignment accuracy, we were mapping reads to the graph they were generated from and we knew were they were supposed to align (thanks to the -a option).

We don't have a "truth set" for the alignment of these reads to either the ref.vg or filtered.vg graphs because they didn't originage from either of these graphs.
Therefore, when we filter using `jq`, we will use `.identity` instead of `.correct`. 

```
vg map -x filtered.xg -g filtered.gcsa -G z.sim --compare -j | jq .identity | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'

vg map -x ref.xg -g ref.gcsa -G z.sim --compare -j | jq .identity | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
```

:question: Which graph has higher read mapping identity?

:question: Do you think that the difference in mapping identity rate is worth the effort of building a graph?

:muscle: Align these reads to the z.vg graph (the graph they were generated from) using the `jq .identity` parameter and compare this number with the numbers above.

:question: Why isn't the alignment identity you just generated equal to 1? How could we make it equal 1?

Let's also look at the sizes of the graphs and indexes.

```
ls -sh *.vg
ls -sh *.gcsa*
```

:question: How do these file sizes seem to scale with the minimum allele frequency cutoff compared with the full variant set?

:question: Given that the .gcsa index is a k-mer index that is used for read alignment, how might the addition of a single variant to the graph alter the number of k-mers in the graph?

:question: If we were looking at a population with a large amount of sequence variation, would it make sense to build a graph that contained as much of this variation as possible? What are the possible outcomes of this approach? What would our index look like?

### Genotyping

We will now explore how we can use read alignment against a graph from a newly sequenced sample for genotyping.
Genotyping is where we identify the specific genomic makeup of a sample.
In this case, we will use a graph constructed from known variants and we will identify which of these variants a new sample has. 
This process is used for a range of applications. 
For example, genetic screening. 
If we know of the specific variants that cause certain diseases, we can use genotyping to determine whether a person/organism will have those diseases or not.

You have been provided with another reference sequence (yeti.fa) and a vcf (yeti.vcf).
There are also reads in a file called `nylamo.reads`.
Nylamo is a yeti sample that we will be genotyping. 

:muscle: build a graph called yeti.vg from the two `yeti` files but include the -a option in your `vg construct` command. 

:question: What does the -a option do?

Index this graph as below:

```
vg index -x yeti.xg -L yeti.vg
vg index -g yeti.gcsa -k 16 yeti.vg
```

Align the reads in nylamo.reads to the graph.

```
vg map -T nylamo.reads -x yeti.xg -g yeti.gcsa > nylamo.gam
```

Now, we will use `vg pack` to determine read support for the graph

```
vg pack -x yeti.xg -g nylamo.gam -o nylamo.pack
```

And finally, we will genotype the sample. 

```
vg call yeti.xg -k nylamo.pack -v yeti.vcf > nylamo_genotype.vcf
```

Open the output file so that you can have a look at it.

The genotype of the sample is found at the start of the final column.
For example, 0/0 means that the sample was homozygous for the reference allele.
1/1 indicates that the sample is homozygous for the alternative allele.
1/0 would mean that the sample was heterozygous with one copy of the reference allele and one copy of the alternate allele.
We won't hve any of these entries though because the yeti genome is a haploid.

:question: Which variants are present in the nylamo sample? Which variants are not?

:question: Can you find the quality scores and read depth for each variant call?

Thanks!
