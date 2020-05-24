# Practicals: Ancient DNA and population genomics

Bastien Llamas \(bastien.llamas@adelaide.edu.au\)
2020-05-26 and 2020-05-29

The two tutorials will be separated into:
1. Data handling (Tuesday 2020-05-26)
2. Population genomics applications (Friday 2020-05-29)


## Population genomics
In a nutshell, population genomics is the study of the genomic composition of populations and how evolution shaped it. Questions in population genomics classically target demographic history (population size through time), gene flow between populations, populations ancestry, or identification of conservation biology units.

While population genetics is usually restricted to a small set of genetic loci, population genomics leverages the large genomic datasets that have become available in recent years and uses up to millions of genetic loci at once.

We are not going to focus on the mathematical aspects of population genomics, but rather on how to manipulate genomic datasets and learn about a few popular tools. I encourage you to read Graham Coop's [course notes](https://github.com/cooplab/popgen-notes/blob/master/popgen_notes.pdf) if you are curious about the underlying mathematical theories.

## VCF format: a reminder
You have previously learnt about several raw or processed high throughput sequencing data formats (e.g., FASTQ, SAM/BAM, VCF). At the end of today's tutorial, you will know how to convert the information contained in a VCF file into other formats compatible with widely-used population genomics programs.

We will use a VCF file of human chromosome 22 from the 1000 Genomes Project (1kGP) that we will save into a working directory in your home directory:
```
# Create working directory
mkdir ~/BIOINF_Tuesday
cd ~/BIOINF_Tuesday
# Download data from the 1kGP public FTP site
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz	> 1kGP_chr22.vcf.gz	
```

Have a look at the VCF file using `zless`. Header lines in a VCF file start with `#` and contain metadata. The last line of the header starts with `#CHROM` and is tab separated. It contains 9 columns of information about the variant sites and then sample names:

||Name|Brief description (from https://en.wikipedia.org/wiki/Variant_Call_Format#The_columns_of_a_VCF)|
|:-|:-|:-|
1|	CHROM|	The name of the sequence (typically a chromosome) on which the variation is being called. This sequence is usually known as 'the reference sequence', i.e. the sequence against which the given sample varies.
2|	POS|	The 1-based position of the variation on the given sequence.
3|	ID|	The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.
4|	REF|	The reference base (or bases in the case of an indel) at the given position on the given reference sequence.
5|	ALT|	The list of alternative alleles at this position.
6|	QUAL|	A quality score associated with the inference of the given alleles.
7|	FILTER|	A flag indicating which of a given set of filters the variation has passed.
8|	INFO|    	An extensible list of key-value pairs (fields) describing the variation. See below for some common fields. Multiple fields are separated by semicolons with optional values in the format: <key>=<data>[,data].
9|	FORMAT|	An (optional) extensible list of fields for describing the samples. See below for some common fields.
+|	SAMPLE|	For each (optional) sample described in the file, values are given for the fields listed in FORMAT

Each line in the body of the VCF file is tab separated and represents a variant site.

Before we move forward, let's see if you can retrieve basic information from a VCF file that will be useful for population genomic analyses. You can use bash commands to answer the below questions, or you can explore the options of `bcftools view` and `bcftools query`. To install `bcftools`, use `conda`:
```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda bcftools
```
#### *Questions*
1. How many variant sites are recorded in the VCF file?
2. How many samples are recorded in the VCF file?

## Converting VCF files to population genomics formats


