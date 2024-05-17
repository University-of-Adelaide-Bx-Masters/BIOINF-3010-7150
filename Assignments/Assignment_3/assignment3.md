# Assignment 3: Structural Variants and Genome Graphs

### By Paul Wang and Chelsea Matthews

---

Data for this assignment is located in the `/data/assignment3/` directory

## Question 1: Copy number variations  (15 marks)

### Part 1a (3 marks)

![Copy Number](../images/Q1_fig1.tumour_CN.png)

The above figure shows the copy number estimates for various genomic regions inside a hypothetical human tumour sample.

List all the copy number variations that you can see from the figure.

### Part 1b (3 marks)

![VAF](../images/Q1_fig2.tumour_VAF.png)

Similarly, the above figure shows the variant allele frequency of variants inside another hypothetical tumour sample. 

List all the CNVs that you can see from this figure.


### Part 1c (3 marks)

Now interpret the two figure together as data from the same sample.

List the CNVs that you can see by integrating both data sets.


### Part 1d (3 marks)

The above figures are from an hypothetical "pure" tumour sample. But in practice, we often get tumour samples which are mixed with some normal/non-tumourous cells.

For example, for the same tumour sample, now there is a percentage of normal cells included in the data, and the CN estimates and VAF graphs now look like this:

![impure CN](../images/Q1_fig3.impure_CN.png)
![impure VAF](../images/Q1_fig4.impure_VAF.png)

Estimate tumour purity of the sample from these figures.

### Part 1e (3 marks)

1. (1 mark) What do the grey bars in the figures represent? Why do they not contain any data points?

2. (2 marks) If we suspect that p-arm of chromosome 8 has fused with the q-arm of chromosome 17. How will you be able to determine whether a fusion event has occured? Will short-read NGS data (WGS, WES or RNA-seq) be able to detect such a fusion event?

## Question 2 - (5 marks)

---

In `/data/assignment3/Q2` you should find these files:

```
mappable_region.fasta
compound_mono_allele.bam
compound_mono_allele.bam.bai
```

Download these files. Launch IGV. 

1. Load `mappable_region.fasta` as the genome file ("Genomes" -> "Load Genome from File...")
2. Load `compound_mono_allele.bam` as the alignment data ("File" -> "Load from File...")


### Question 2.1

Locate and list any breakpoints you can see in the data.

[2 marks]

### Question 2.2

Identify all SV events and their associated breakpoints. 
Show the steps and reasonings for your answers. Include diagrams if you think it helps. If you want to use hand-drawn diagram, just take and submit a photo of your drawing, but make it's clearly legible.

[3 marks]

## Question 3 - MSA Pan-genome Graphs

---

You have been provided with a file called `svs.vg` in the `~/data/assignment3/graphs` folder that contains a graph constructed using the multiple sequence alignment technique from a reference sequence and 4 other well assembled haplotypes.
Therefore, there are 5 paths embedded within the graph.

### Question 3.1 

Visualise the graph in any way you like and provide a screenshot/image of the most interesting looking region.
Also provide the command you used to generate the visualisation. 

[2 marks]

### Question 3.2

Identify and describe the structural variants you see in each of the four haplotypes embedded in the graph.  

[3 marks]

### Question 3.3 - Conceptual

If we were to align reads generated from sample 1 and sample 2 to the linear reference genome (as represented by the reference path through the graph above), which do you think would have a higher mapping rate? Why?

[2 marks]

## Question 4 - Pan-genome Variant Graphs

You have been provided with a reference sequence `hippogryph.fasta` as well as a .vcf containing known variants found in the hippogryph population in `hippogryph.vcf` in the `~/data/assignment3/graphs/` diretory.

A hippogryph is a mythical creature with the body of a horse and the head and wings of an eagle. 

Your task is to genotype a newly sequenced sample (Frayfeather) which will help us to understand more about the phenotypic characteristics of frayfeather. 

### Question 4.1

Construct a hippogryph pan-genome graph using the `hippogryph.fa` reference sequence and the .vcf.

When constructing this graph, keep in mind that we will be using it for genotyping.

Visualise this graph in a format of your choice and screenshot the first part of the graph up to and including the first variant.
Any visual form is acceptable but you should also include the code that you used to create this visualisation.

[3 marks]

### Question 4.2

 Align the provided reads in `frayfeather.fastq` to the graph. 
Report your alignment command and the average read mapping identity using `jq .identity`.

Would you expect this value to be higher or lower if the graph contained only the **least** common variants found in the hippogryph population. 

[2 marks]

### Question 4.3

Genotype the frayfeather sample using your read alignments to produce a vcf.

**IMPORTANT HIPPOGRYPH INORMATION**

- Hippogryphs with the alternate allele for variant hip748 have particularly strong wings and can fly very long distances. 
- Hippogryphs with the reference allele for variant hip932 love to eat fish and seafood while the alternate allele gives them a seafood intolerance. 

Inspect the .vcf you produced for Frayfeather. 
Did he have the reference allele or alternative allele for variants hip748 and hip932? 

Choose the option from below that best describes Frayfeather:

- **A.** Frayfeather is an excellent choice for long flights over the sea!
- **B.** Frayfeather is great for long flights over the sea as long as you make sure to pack his lunch.
- **C.** Frayfeather would prefer to stay at home in the mountains and eat goats.

[2 marks]

### Question 4.4 - Conceptual 

What are the potential impacts (positive and negative) of including all of the known genomic variants for a population in a single pan-genome graph?

Do you think that this is advisable? Why or why not?
 
[3 marks]

### Question 4.5 - Conceptual

What is reference bias and how can using a pan-genome graph reduce reference bias?

[2 marks]

### Question 4.6 - Conceptual

In transcriptomics experiments, we align reads generated from RNA to a reference genome and determine the read depth for each gene in the reference.
Greater read depth for a gene is equated with higher levels of expression than a gene with lower read depth. 
In this way, we can make an estimate of relative gene expression. 
This system relies on the assumption that reads are aligning equally well to all genes. 
How could using a pan-genome graph potentially increase the accuracy of this type of analysis?

[1 mark]

Total: 40 Marks

