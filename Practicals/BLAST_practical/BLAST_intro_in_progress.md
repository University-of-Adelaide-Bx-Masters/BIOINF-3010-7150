# Introduction to BLAST 

## 1. Background

NCBI BLAST is a suite of programs that will find local alignments of query sequences to database entries. We will use the command line version of BLAST to identify matching nucleotide or protein coding regions in a human genomic sequence and in the SwissProt database. You can find out about all things BLAST for command line BLAST [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/) or in the [BLAST Book](https://www.oreilly.com/library/view/blast/0596002998/).  

BLAST can align a variety of sequence types:
- nucleotide sequences to nucleotide sequences (BLASTN)
- protein sequences to protein sequences (BLASTP)
- translated nucleotide sequences to protein sequences (BLASTX)
- protein sequences to translated nucleotide sequences (TBLASTN)
- translated nucleotide sequences to translated nucleotide sequences (TBLASTX)

There are a number of parameters that can be varied to increase sensitivity or speed up a search at the expense of sensitivity. There also numerous [parameters](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a) to modify BLAST behaviour and customise output.

## 2. The BLAST algorithm

BLAST algorithm steps:  
- Generate a word list from the query sequence (W letters)  *Smaller values of W give a more sensitive search as the expense of speed*
- For each word, scan database for high scoring word matches (>T)  *T parameter only for protein alignments, smaller values of T give a more sensitive search at the expense of speed*  

[Effect of T on word matches - Fig 5.3 BLAST Book](https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/blastp_T_param.png)

- For each word, find high scoring neighbours and make word clusters  

[Word clusters - Fig 5.4 BLAST book](https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/blast_word_size.png)

- Extend words/clusters to form High Scoring Pairs (HSP)  
- Evaluate significance for each HSP  
- Get local alignment for each HSP  
- Report matches with low `e-value`

## 3. Having a BLAST

This practical aims to familiarise you with the use of NCBI BLAST as a tool for annotation. You will use your VMs for this. Your first task will be to install BLAST and download the databases you will require. 

### 3.1 Install NCBI BLAST

Open a terminal and at the bash command line prompt type the following to create a new `conda` environment for `blast`. 

```bash
conda create --name blast
conda activate blast
conda install -c bioconda blast
```
This will install NCBI BLAST and all its dependencies.

Once you have installed NCBI BLAST you will need to download the data you will need for the practical.

### 3.2 The data

The files you need are in `~/data/`.

See the `README_BLAST.txt` file for descriptions of the data files.

Make sure you create a working directory, along  with useful sub-directories for the practical. 

```bash
mkdir -p ~/BLAST_practical/dbs
mkdir -p ~/BLAST_practical/queries
mkdir -p ~/BLAST_practical/results
```

### 3.3 Prepare the BLAST databases

BLAST searches a special database of nucleotide sequences that have been broken into `kmers` for faster searching. It finds matching words in the database and then extends the matches to create a local alignment. You will need to format the two BLAST databases that you will use.

You will need to decompress the hg 38.fa uniprot_sprot.fasta files:

```bash
unpigz -p 2 uniprot_sprot.fasta.gz
unpigz -p 2 hg38.fa.gz
```
Once you have done this, you will need to index/format these files so that BLAST can search them.

The syntax for `makeblastdb` is as follows:

`makeblastdb -in <reference.fa> -dbtype nucl -parse_seqids -out <database_name> -title "Database title"`
 
In this example your input file is a nucleotide sequence file `reference.fa` so you need to specify the database type `-dbtype nucl`. You should use the `-parse-seqids` flag as it preserves the sequence identifiers in the input file. The `-out` flag specifies the name for db files and the `-title` flag is optional. 

so you will

```bash
makeblastdb -in uniprot_sprot.fasta -dbtype 'prot' -parse_sequids -out sprot
```

This will generate three files that BLASTX uses.

You will then need to index/format the human hg38 chromosome sequences so that BLASTN can search them.

```bash
makeblastdb -in hg38.fa -dbtype 'nucl' -parse_seqids -out hg38
```
Now you are ready to use BLAST.

### 3.4 Have a BLAST

For quick BLASTN help you can type:

```bash
blastn -help
```

For BLASTX:

```bash
blastx -help
```

Try these to see what the allowed syntax, flags and parameters are.  You will probably want to pipe this to `less`.

For detailed documentation for all things BLAST see [here](https://www.ncbi.nlm.nih.gov/books/NBK1762/) and for detailed command line options and flags see [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a)

A typical command line might look something like this:

```bash
blastn -query [file.fasta] -task [megablast] -db [database file]  -outfmt [0 through 17] -out [outputfile]
```

- I suggest using outfmt 7 and/or 17, 7 gives you a tab delimited file, 17 gives you a .sam file. 

Output formatting options shown below.

> -outfmt `<number>`  
   alignment view options:  
     0 = Pairwise,  
     1 = Query-anchored showing identities,  
     2 = Query-anchored no identities,  
     3 = Flat query-anchored showing identities,  
     4 = Flat query-anchored no identities,  
     5 = BLAST XML,  
     6 = Tabular,  
     7 = Tabular with comment lines,  
     8 = Seqalign (Text ASN.1),  
     9 = Seqalign (Binary ASN.1),  
    10 = Comma-separated values,  
    11 = BLAST archive (ASN.1),  
    12 = Seqalign (JSON),  
    13 = Multiple-file BLAST JSON,  
    14 = Multiple-file BLAST XML2,  
    15 = Single-file BLAST JSON,  
    16 = Single-file BLAST XML2,  
    17 = Sequence Alignment/Map (SAM),  
    18 = Organism Report  

Note that when using output formats >4 some options for choosing the number and description of hits do not work or are incompatible (see below). 

> -num_descriptions `<Integer, >=0>`  
   Number of database sequences to show one-line descriptions for  
   Not applicable for outfmt > 4  
   Default = `500'  
    * Incompatible with:  max_target_seqs  
>    
> -num_alignments `<Integer, >=0>`  
   Number of database sequences to show alignments for  
   Default = `250'  
    * Incompatible with:  max_target_seqs  
>    
> -max_target_seqs `<Integer, >=1>`  
   Maximum number of aligned sequences to keep   
   (value of 5 or more is recommended)  
   Default = `500'  
    * Incompatible with:  num_descriptions, num_alignments  


#### 3.4.1 Alignments to Swissprot proteins

You can try the following 

> A "vanilla" BLAST command line that give the basic text output from BLAST

```bash
blastp -query ~/BLAST_practical/queries/protein_query.fa -task blastp -db ~/BLAST_practical/dbs/sprot
```

> A command line that specifies a set of custom columns to report

```bash
blastp -query ~/BLAST_practical/queries/protein_query.fa -task blastp -db ~/BLAST_practical/dbs/sprot -outfmt "7 delim=  qaccver qlen sallgi sallacc slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

These will of course just dump everything to `stdout` but you know how to cope with that.

> A command line that writes output to a file and specifies the number of threads to run and limits output to 50 targets. 

```bash
blastx -query ~/BLAST_practical/queries/hg38_gene_query.fasta -task blastx -db ~/BLAST_practical/dbs/sprot -num_threads 2 -max_target_seqs 50 -out ~/BLAST_practical/results/humgene_blastx_sprot.txt -outfmt 7 
```
Call your output file whatever you like, as long as it makes sense to you. 

Once BLASTX has completed you can look at your output using "head", "less", "more" or "cat" or open it with a text editor. 

There will be quite a few hits. You can reduce them to a manageable level by parsing the output to find only the alignments with human proteins.

```bash
grep HUMAN H15_blastx_sprot.txt | less
```