
# BLAST Practical: Part 2 - Dave Adelson
{:.no_toc}

* TOC
{:toc}

Annotation can be tricky. In this practical you will use BLAST to identify the protein encoded by a gene (sequence supplied) and to determine the principal repeat sequence in the interval and then estimate the number of those repeats and the number of base pairs they contribute to the human genome.  

## 1. Background

To find protein coding regions (exons) in the human gene sequence we have provided, you will use BLASTX to find local alignments of your sequence with proteins in Uniprot-swissprot, a curated protein database. 

To find repeats in the human gene sequence we have provided, you will use BLASTN to find local alignments of your sequence with human repeat consensus sequences from RepBase, a curated repeat database. 

To see how difficult it can be to deal with the large numbers of repeats in the genome you will extract one repeat interval identified in your repeat BLAST output using Samtools and you will use that repeat sequence to search the human genome for alignments using BLASTN. 

### Questions you should answer after following the instructions below:

- Identify the protein coding gene present in your query sequence. 

- What is the number of repeat intervals found in your query sequence, and the most dominant family of repeats in your output?

- How many BLAST hits did you obtain in the human genome using your repeat interval?

- What difficulties and problems have you faced in this course of this practical?

- How these are likely to affect the process of genome annotation?

## 2. Practical instructions.

### 2.1  Prepare the data

You will then need to index/format the `humanrReps.fa.gz` consensus sequence so that BLASTN can search it. Rebember that `makeblastdb` will not accept `.gz` compressed files. 

#### Activate NCBI BLAST environment

Open a terminal and at the bash command line prompt type the following to activate the `conda` environment for `blast`. 

```
source activate bioinf
```

Now you are ready to use BLAST.

### 2.2 Alignments to Swissprot proteins

Reminder: general syntax for BLAST searches is as follows:  

```bash
blastn -query [file.fasta] -task [blastn] -db [database file]  -outfmt [0 through 17] -out [outputfile]
```
For the next step you will use BLASTX, not BLASTN so that you can query the protein database using a nucleotide (gene) query. The gene query is the full sequence of a gene from the human genome. 

- I suggest using outfmt 7 and 17, 7 gives you a tab delimited file, 17 gives you a .sam file. 

```bash
blastx -query ~/Project_4/queries/hg38_gene_query.fasta -db ~/Project_4/dbs/sprot \
-num_threads 2 -out ~/Project_4/results/HumGene_blastx_sprot.txt \
-outfmt "7 delim= qaccver qlen sallgi sallacc stitle slen pident length \
mismatch gapopen qstart qend sstart send evalue bitscore"
```
Call your output file whatever you like, as long as it makes sense to you. 

Once BLASTX has completed you can look at your output using "head", "less", "more" or "cat" or open it with a text editor. 

There will be quite a few hits. You can reduce them to a manageable level by parsing the output to find only the alignments with human proteins.

```bash
grep "Homo sapiens" ~/Project_4/results/HumGene_blastx_sprot.txt | less
```

### 2.3 Alignments to human repeat consensus sequences

You will need to make a `BLAST` index for the `~/dbs/humanReps.fa` file before carrying out the next search. Look at the instructions for the previous practical to figure out how to do this. 

```bash
blastn -query ~/Project_4/queries/hg38_gene_query.fasta -task blastn -db ~/Project_4/dbs/humrep -out ~/Project_4/results/gene_blastn_humrep.txt -outfmt 7
```

This will identify all the repeat sequence intervals in your gene sequence. 

To determine the most abundant repeat type in your output you can try:
- just scroll through the output and eyeball it
- use `cut` , `sort` and `uniq` to list all the repeat types
- use `grep -c` to count some of the repeat types to get an objective assessment of how many insertions there are for every repeat type. *Hint: when using `grep` use the shortest search pattern you can, to group repeats of the same type into the count*.

In order to obtain a human repeat sub-sequence for the most abundant repeat type from `hg38_gene_query.fasta` you will need to use [samtools-faidx](https://www.htslib.org/doc/samtools-faidx.html). 

- You will need to identify the coordinates of the repeat interval that you will use to retrieve the sequence. 

- Identify the coordinates of the repeat interval you want to retrieve by inspecting the text output file from the above `blastn` search and select an interval from a robust (*longest or almost longest `alignment length`  with high `bitscore` and low `evalue`*) alignment **for the most abundant type of repeat** in your output. 

**I have used arbitrary coordinate values 12045-12345 in the example below, you will need to use the coordinates you selected.**

`samtools` should be already installed on your VMs and you should be able use it to extract the subsequence you want from the genome. Remember we covered `sam`/`bam` file formats and the use of `samtools` more extensively in an earlier practical. For now we are using this just retrieve a sequence from the database that matched our query. 

```bash

samtools faidx ~/Project_4/queries/hg38_gene_query.fasta hg38:12045-12345 > ~/Project_4/queries/hg38_12045-12345.fasta
```

### 2.4 Alignment of your human repeat sub-sequence to the human genome

This may take a while to run, so be patient.  

```bash

blastn -query ~/Project_4/queries/hg38_12045_12345.fasta -task blastn -db ~/Project_4/dbs/hg38 -num_threads 2 -out ~/Project_4/results/hg38_repeats.txt -outfmt 7
```

The output will be very large, so do not open with the text editor. You can see how many hits you have by using:

```bash
head -n5 ~/Project_4/results/hg38_repeats.txt

```
This gives you an estimate of the number of insertions of that repeat/transposon type in the genome.

You can get the sum of the values for `alignment length` to determine the exact number of base pairs contributed by this repeat to the genome. *Hint: you can use [awk](https://www.gnu.org/software/gawk/manual/gawk.html) with this syntax:* `awk` is a very useful program that allows you to select particular records in a file and perform operations on them. 

```bash
awk '{sum+=$n} END {print sum}' [input file]
```
Where `n` is the column number in your input file you want to sum.

Now you should be able to answer all the questions posed at the beginning of section 2. 

**In Assignment 1, question 13 requires you to understand and be able to identify a repeat sequence.** Make sure you understand all the steps you have carried out in sections 2.3 and 2.4. 
