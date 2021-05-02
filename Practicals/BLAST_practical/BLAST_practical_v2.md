
# BLAST Practical

## 1. Background

To find protein coding regions (exons) in the human sequence we have provided you will use BLASTX to find local alignments of your sequence with proteins in Uniprot-swissprot, a curated protein database. 

To find repeats in the human sequence we have provided you will use BLASTN to find local alignments of your sequence with human repeat consensus sequences from RepBase, a curated repeat database. 

To see how difficult it can be to deal with the large numbers of repeats in the genome you will extract one repeat interval identified in your repeat BLAST output using Samtools and you will use that repeat sequence to search the human genome for alignments using BLASTN. 

### Questions:

- Identify the protein coding gene present in your query sequence. 

- What is the number of repeat intervals found in your query sequence, and the most dominant family of repeats in your output?

- How many BLAST hits did you obtain in the human genome using your repeat interval?

- What difficulties and problems have you faced in this course of this practical?

- How these are likely to affect the process of genome annotation?

## 2. Practical instructions.

This practical aims to familiarise you with the use of NCBI BLAST as a tool for annotation. You will use your VMs for this. Your first task will be to install BLAST and download the databases you will require. 

### 2.1 Install NCBI BLAST

Open a terminal and at the bash command line prompt type the following:


```bash
conda activate assembly
```

```bash
conda install -c conda-forge -c bioconda blast

```
Once you have installed NCBI BLAST you will need to download the data you will need to carry out the practical.

### 2.2 Download the data


Some of the files you need are in a zip archive called BLAST_prac_data.zip in your home directory.

Make sure you create a directory for the practical and put the zip file there prior to unzipping. You should make sure that all of your data files are in this directory. 

To obtain the human genome blastdb do the following:

```bash
conda install -c conda-forge gdown
```
```bash
gdown https://drive.google.com/uc?id=1LKYFxAovMTcSPk54ZfSoL7kHhTlnCKm_
```
I suggest you run `gdown` from that newly created directory to obtain the hg38.nsq file.

All analyses will be conducted using the bash command line in this directory. 

Once you have all the files you can consult the README.txt file to learn what you have downloaded. 



### 2.3 Prepare the data

You will need to decompress the Swissprot fasta file:

```bash
gunzip uniprot_sprot.fasta.gz
```
Once you have done this, you will need to index/format it so that BLAST can search it.

```bash
makeblastdb -in uniprot_sprot.fasta -dbtype 'prot' -out sprot
```

This will generate three files that BLASTX uses.

You will then need to index/format the human repeat fasta consensus sequences so that BLASTN can search them.
makeblastdb -in chr15.fa -dbtype 'nucl' -out chr15
```bash
makeblastdb -in humrep.ref -dbtype 'nucl' -out humrep
```
Now you are ready to use BLAST.

### 2.4 Have a BLAST

For quick BLASTN help you can type:

```bash
blastn -help
```

For BLASTX:

```bash
blastx -help
```

Try these to see what the allowed syntax, flags and parameters are.  

For detailed documentation for all things BLAST see: https://www.ncbi.nlm.nih.gov/books/NBK1762/


#### 2.4.1 Alignments to Swissprot proteins

General syntax for BLAST searches is as follows:\

```bash
blastn -query [file.fasta] -task [blastn] -db [database file]  -outfmt [0 through 17] -out [outputfile]
```

- I suggest using outfmt 7 and 17, 7 gives you a tab delimited file, 17 gives you a .sam file. 

I used the following. *You will need to reduce the num_threads to something that will run on your VM*. 

```bash
blastx -query Human15gene.fasta -task blastx -db sprot -num_threads 6 -out H15_blastx_sprot.txt -outfmt 7
```
Call your output file whatever you like, as long as it makes sense to you. 

Once BLASTX has completed you can look at your output using "head", "less", "more" or "cat" or open it with a text editor. 

There will be quite a few hits. You can reduce them to a manageable level by parsing the output to find only the alignments with human proteins.

```bash
grep HUMAN H15_blastx_sprot.txt | less
```

#### 2.4.2 Alignments to human repeat consensus sequences

```bash
blastn -query Human15gene.fasta -task blastn -db humrep -out H15_blastn_humrep.txt -outfmt 7
```
You can then run this again to get the .sam output

```bash
blastn -query Human15gene.fasta -task blastn -db humrep -out H15_blastn_humrep.sam -outfmt 17
```

In order to obtain a human repeat sub-sequence from Human15gen.fasta you will need to use samtools-faidx https://www.htslib.org/doc/samtools-faidx.html. You will need to identify the coordinates of the repeat interval that you will use to retrieve the sequence. Do this by inspecting the text output file from above and selecting an interval from a robust looking alignment for the most dominant type of repeat in your output. 

I have used arbitrary coordinate values in the example below, you will need to use your own coordinates.

```bash 
samtools faidx Human15gene.fasta hg38:12345-12345 > hg38_12345_12345.fasta
```

#### 2.4.3 Alignment of your human repeat sub-sequence to the human genome

This may take a while to run, so be patient. Use as many threads as you can get away with for this in order to make it run as fast as possible. 

```bash
blastn -query hg38_12345_13345.fasta -task blastn -db hg38 -num_threads 6 -out hg38_repeats.txt -outfmt 7
```

The output will be very large, so do not open with the text editor. You can see how many hits you have by using:

```bash
head -n5 hg38_repeats.txt
```

