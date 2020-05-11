
# Clinical Genomics & Pedigree Analysis

High-throughput sequencing is currently used ubiquitously in identifying the cause of a large range of genetic diseases.
While single gene and well-known Mendelian genetic disorders, such as sickle-cell anemia, Tay–Sachs disease and cystic fibrosis, can be identified with simple diagnostic techniques, whole genome (WGS) and exome (WES/WXS) sequencing can be used to identify and study a wide variety of inherited traits.
Cost used to be a barrier for using high-throughput sequencing approaches, but now it is possible to sequence a patient in [under 27 hours for less than ~US$1,000 per sample](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0221-8). 

The resolution to which variants can be assessed with WGS means that it is also powerful approach to identify new genetic variants that cause disease.
Using WGS or WES on thousands of samples, we can now establish fine-grained association maps for more and more complex diseases that are likely to be caused by the action of hundreds of genes.

The current clinical workflow works a lot like this:

![Priest. (2017). _Curr Opin Pediatr_. 29(5): 513–519.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5590671/bin/nihms899674f1.jpg)

In Tuesday's tutorial, we discussed annotation of identified variants in three samples, and today we will be looking at family inheritance patterns.
As you may have noticed, the three samples sequenced in our data are related, and form a "trio" (mother-father-daughter). 
Trios and other members of an affected individual's family are often sequenced in clinical genetics, allowing clinicians to establish the inheritance pattern of the trait or identify new _de novo_ mutations that may have arisen independently of the parents.

## This week's tutorial

This week's tutorial is liberally taken from two tutorials written by Aaron Quinlan & his group at University of Utah.

- [Identifying dominant gene candidates with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-Dominant-Tutorial.pdf)
- [Identifying recessive gene candidates with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-Recessive-Tutorial.pdf)

Both of these tutorials do a really good job at introducing the program `gemini`, which is used quite a bit across clinical genetics studies.
[Gemini](https://gemini.readthedocs.io/en/latest/) is a database system that can read in VCF information and family/pedigree information, to enable database querying and clinical genetics analyses.
Information in gemini is stored in database system called SQL.
SQL (pronounced "ess-que-el") stands for Structured Query Language, and is a popular database system in many industries and enables the store of organised information that can be accessed by queries.
It comes in many flavours that you might have heard before, including `MySQL`, `SQLite` and `PostgreSQL`.

## Cohort databases

Lets make some databases!
Gemini can take the VCF file and sample information in the form of a ped file (short for pedigree).
The ped file is actually a standard metadata information file that was developed in the [genetics application `PLINK`](http://zzz.bwh.harvard.edu/plink/data.shtml#ped).
This program is used extensively for genome-wide association studies, and was developed in the era of genotyping arrays rather than WGS approaches. 

```
# Create the gemini db for recessive study
gemini load --cores	4 -v trio.trim.vep.vcf.gz -t VEP \
        --skip-gene-tables -p recessive.ped trio.trim.vep.recessive.db

# Create the gemini db for dominant study
gemini load --cores 4 -v trio.trim.vep.vcf.gz -t VEP \
        --skip-gene-tables -p dominant.ped trio.trim.vep.dominant.db
```


