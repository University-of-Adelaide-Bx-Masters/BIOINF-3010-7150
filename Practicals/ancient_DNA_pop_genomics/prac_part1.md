# Practicals: Ancient DNA and population genomics

Bastien Llamas \(bastien.llamas@adelaide.edu.au\)
2020-05-26 and 2020-05-29

---
The two tutorials will be separated into:
1. Data handling (Tuesday 2020-05-26)
2. Population genomics applications (Friday 2020-05-29)

---
# Day 1: Data handling (Tuesday 2020-05-26)

## Population genomics
In a nutshell, population genomics is the study of the genomic composition of populations and how evolution shaped it. Questions in population genomics classically target demographic history (population size through time), gene flow between populations, populations ancestry, or identification of conservation biology units.

While population genetics is usually restricted to a small set of genetic loci, population genomics leverages the large genomic datasets that have become available in recent years and uses up to millions of genetic loci at once.

We are not going to focus on the mathematical aspects of population genomics, but rather on how to manipulate genomic datasets and learn about a few popular approaches and tools. I encourage you to read Graham Coop's [course notes](https://github.com/cooplab/popgen-notes/blob/master/popgen_notes.pdf) if you are curious about the underlying mathematical theories.


## VCF format: a reminder
You have previously learnt about several raw or processed high throughput sequencing data formats (e.g., FASTQ, SAM/BAM, VCF). In particular, you should now know that [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) files contain information about variants found at specific positions in a reference genome.

At the end of today's tutorial, you will know how to convert the information contained in a VCF file into other formats compatible with widely-used population genomics programs.

We will use a VCF file of human chromosome 22 from the 1000 Genomes Project (1kGP) that we will save into a working directory in your home directory:
```
# Create working directory
mkdir ~/BIOINF_Tuesday
cd ~/BIOINF_Tuesday
# Download data from the 1kGP public FTP site (file size: 214453750 bytes)
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz	> 1kGP_chr22.vcf.gz	
```
Although you could use your own scripts to parse VCF files and analyse variant calls, several tools have already been developed for your convenience. In particular, [BCFtools](http://samtools.github.io/bcftools/bcftools.html) is a set of useful utilities to manipulate variant calls in VCF files. You can install it easily with the conda package management system:
```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda bcftools
```

#### VCF meta-information and header lines
Have a look at the VCF file using `zless`. Meta-information lines start with `##` and contain various metadata. The header line starts with `#` and is tab separated. It contains 9 columns of information about the variants, and then one column per sample name:

||Name|Brief description (from [Wikipedia](https://en.wikipedia.org/wiki/Variant_Call_Format#The_columns_of_a_VCF))|
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

Have a closer look at how the information in the [INFO](https://en.wikipedia.org/wiki/Variant_Call_Format#Common_INFO_fields) and [FORMAT](https://en.wikipedia.org/wiki/Variant_Call_Format#Common_FORMAT_fields) fields is commonly coded. The 1kGP VCF datasets also contain some project-specific keys explained in a file that can be downloaded.
```
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/README_vcf_info_annotation.20141104
```

#### VCF body
The body of the VCF file is tab separated. Each line represents a unique variant site.

#### Other useful 1kGP metadata
You can  download sample details from the 1kGP FTP site to learn about population of origin and sex of each individual.
```
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

---
Before we move forward, let's see if you can retrieve basic information from a 1kGP VCF file that will be useful for population genomic analyses. You can use `bcftools view`, `bcftools query`, or bash commands to answer the below questions.

#### *Questions*
1. Using `bcftools view` or bash commands, determine how many variant sites are recorded in the VCF file.
2. Using `bcftools query` or bash commands, determine how many samples are recorded in the VCF file.
3. The INFO fields contain a lot of information. In particular for the first variant: determine how many samples have data, how many ALT alleles are reported, and what the frequency of the ALT allele is.
4. Same as above for position 16051249 (see the BCFtools manual for region or target formatting).
5. Looking at the information contained in the `FORMAT` field in the body of the VCF file,

---
## Converting VCF files into population genomics formats
The VCF file contains a lot of information that can be very useful for clinical genomics. However, population genomics applications only need a subset of the information in VCF file, i.e., variants genomic coordinates, variants ID, reference (REF) and alternative (ALT) alleles, and sample genotypes.

Let's have a look at the first variant in our VCF file:
```
bcftools view -H 1kGP_chr22.vcf.gz | head -1
```
#### *Questions*
3. What are the REF and ALT alleles?
4. Given REF and ALT alleles, and knowing that the genotypes are phased (use of `|` instead of `/` to separate alleles), what are the possible genotypes?
4. How many samples are recorded in the VCF file?


