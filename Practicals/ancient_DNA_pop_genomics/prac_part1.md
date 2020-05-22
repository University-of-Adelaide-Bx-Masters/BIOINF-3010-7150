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

## Converting VCF files to population genomics formats
You have previously learnt about several raw or processed high throughput sequencing data formats (e.g., FASTQ, SAM/BAM, VCF). At the end of today's tutorial, you will know how to convert the information contained in a VCF file into other formats compatible with widely-used population genomics programs.

We will use the same VCF file as during the Clinical Genomics practicals. But first we need to install bcftools again.

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda bcftools
```

Create working directory in your home directory
mkdir ~/BIOINF_Tuesday

What is Conda? (from https://docs.conda.io/en/latest/)
Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language.

Add the bioconda channel
cd
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

http://10.150.8.155:8787/
