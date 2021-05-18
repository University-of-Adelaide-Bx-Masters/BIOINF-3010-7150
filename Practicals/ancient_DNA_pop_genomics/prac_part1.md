# Practicals: Ancient DNA and population genomics

Bastien Llamas \(bastien.llamas@adelaide.edu.au\)
2021-05-28 and 2021-06-01

---
The two tutorials will be separated into:
1. Data handling (Friday 2021-05-28)
2. Population genomics applications (Tuesday 2021-06-01)

Icons are used to highlight sections of the tutorials:

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Information  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Hands-on tasks  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> Questions  

---
# Day 1: Data handling (Friday 2021-05-28)

## Tutorial outcomes

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> At the end of today's tutorial, you will know how to convert the information contained in a VCF file into other formats compatible with widely-used population genomics programs.

## Population genomics
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> In a nutshell, population genomics is the study of the genomic composition of populations and how evolution shaped it. Questions in population genomics classically target demographic history (population size through time), gene flow between populations, populations ancestry, or identification of conservation biology units.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> While population genetics is usually restricted to a small set of genetic loci, population genomics leverages the large genomic datasets that have become available in recent years and uses up to millions of genetic loci at once.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> We are not going to focus on the mathematical aspects of population genomics, but rather on how to manipulate genomic datasets and learn about a few popular approaches and tools. I encourage you to read Graham Coop's [course notes](https://github.com/cooplab/popgen-notes/blob/master/popgen_notes.pdf) if you are curious about the underlying mathematical theories.


## VCF format: a reminder
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> You have previously learnt about several raw or processed high throughput sequencing data formats (e.g., FASTQ, SAM/BAM, VCF). In particular, you should now know that [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) files contain information about variant calls found at specific positions in a reference genome.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We will use a VCF file of human chromosome 22 from the 1000 Genomes Project (1kGP) that we will save into a working directory in your home directory:
```bash
# Create working directory
mkdir ~/BIOINF_Tuesday
cd ~/BIOINF_Tuesday

# Download compressed VCF file and its index from the 1kGP public FTP site (VCF file size: 214453750 bytes)
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz	> 1kGP_chr22.vcf.gz
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi	> 1kGP_chr22.vcf.gz.tbi
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Although you could use your own scripts to parse VCF files and analyse variant calls, several tools have already been developed for your convenience. In particular, [`BCFtools`](http://samtools.github.io/bcftools/bcftools.html) is a set of useful utilities to manipulate variant calls in VCF files. Install it easily with the conda package management system.  

```bash
# Activate a conda environment that contains other software we will use today
conda activate variation
# Install bcftools (if not already installed)
conda config --add channels bioconda
conda install -c bioconda bcftools
```

#### VCF meta-information and header lines
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the compressed VCF file using `zless`. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Meta-information lines start with `##` and contain various metadata. The header line starts with `#` and is tab separated. It contains 9 columns of information about the variant calls, and then one column per sample name:

||Name|Brief description (from [Wikipedia](https://en.wikipedia.org/wiki/Variant_Call_Format#The_columns_of_a_VCF))|
|:-|:-|:-|
1|	`CHROM`|	The name of the sequence (typically a chromosome) on which the variation is being called. This sequence is usually known as 'the reference sequence', i.e. the sequence against which the given sample varies.
2|	`POS`|	The 1-based position of the variation on the given sequence.
3|	`ID`|	The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.
4|	`REF`|	The reference base (or bases in the case of an indel) at the given position on the given reference sequence.
5|	`ALT`|	The list of alternative alleles at this position.
6|	`QUAL`|	A quality score associated with the inference of the given alleles.
7|	`FILTER`|	A flag indicating which of a given set of filters the variation has passed.
8|	`INFO`|	An extensible list of key-value pairs (fields) describing the variation. See below for some common fields. Multiple fields are separated by semicolons with optional values in the format: <key>=<data>[,data].
9|	`FORMAT`|	An (optional) extensible list of fields for describing the samples. See below for some common fields.
+|	`SAMPLE`|	For each (optional) sample described in the file, values are given for the fields listed in `FORMAT`

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a closer look at how the information in the [INFO](https://en.wikipedia.org/wiki/Variant_Call_Format#Common_INFO_fields) and [FORMAT](https://en.wikipedia.org/wiki/Variant_Call_Format#Common_FORMAT_fields) fields is commonly coded. The 1kGP VCF datasets also contain some project-specific keys explained in a file that can be downloaded.
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/README_vcf_info_annotation.20141104
```

#### VCF body
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The body of the VCF file is tab separated. Each line represents a unique variant site.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Before we move forward, let's see if you can retrieve basic information from a 1kGP VCF file that will be useful for population genomic analyses.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*
- 1) Using `bcftools view` or bash commands, determine how many variant sites are recorded in the VCF file.
- 2) Using `bcftools query` or bash commands, determine how many samples are recorded in the VCF file.
- 3) The `INFO` fields contain a lot of information. In particular for the first variant position in the file: determine how many samples have data, how many ALT alleles are reported,  what the frequency of the ALT allele is globally, and what the frequency of the ALT allele is in East Asians.
- 4) Same as question 3 for variant position 16051249 (see the [BCFtools manual](http://samtools.github.io/bcftools/bcftools.html) for region or target formatting).
- 5) How many alternative alleles are observed at position 16050654?    
- 6) Looking at the information contained in the `FORMAT` field in the body of the VCF file, what kind of data is stored in the VCF file for each sample?  

---

#### Other useful 1kGP metadata
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Download sample details from the 1kGP FTP site to learn about population of origin and sex of each individual.
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*  

- 7) Using bash commands on the panel file you just downloaded, determine how many different populations and super-populations are represented in the 1kGP dataset.  
- 8) How many individuals are in each super-population?
---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> You can learn more about the populations in the 1kGP [here](https://www.internationalgenome.org/faq/which-populations-are-part-your-study/).
## Converting VCF files into population genomics formats
### Rationale
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> A VCF file may contain a lot of information (e.g. variant annotation) that can be very useful for clinical genomics. This was the case when you looked at a trio data with Jimmy Breen. However, population genomics applications only need a subset of the information in VCF file, i.e., variant genomic coordinates, variant ID, reference (REF) and alternative (ALT) alleles, and sample genotypes. This is what you typically find in the 1kGP data.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Genotypes can be coded differently depending on ploidy (in humans, diploid for autosomes or chrX in females, haploid for mitochondrial genomes and chrY), the number of alternate alleles, whether the genomes are phased (i.e., alleles on maternal and paternal chromosomes are identified) or not, and homo/heterozygosity. For convenience, 0 is used to code REF, and 1, 2, 3, etc are used to code ALT.
||Example|
|:-|:-|
|**Haploid, 2 alleles**|2 genotypes: 0, 1|
|**Diploid, not phased, 2 alleles**|3 genotypes: 0/0, 0/1, 1/1|
|**Diploid, phased, 2 alleles**|4 genotypes: 0\|0, 0\|1, 1\|0, 1\|1|
|**Haploid, n alleles**|n genotypes: e.g., 0, 1, 2, 7|
|**Diploid, not phased, n alleles**|n! genotypes: e.g., 0/0, 1/4, 0/2, 3/3|
|**Diploid, phased, n alleles**|n<sup>2</sup> genotypes: e.g., 0/0, 1/4, 0/2, 3/3|

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the first variant in the VCF file.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*  

- 9) What are the REF and ALT alleles?
- 10) Given REF and ALT alleles found when answering question 9, and knowing that the genotypes are phased, what are the possible genotypes with nucleotides and 1kGP coding?
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> You saw that the VCF genotype information can be very detailed. However, all we need usually for population genomics is a table of samples and variant calls, where the genotype information is coded so it can be parsed easily and file size remains as small as possible (imagine storing and parsing whole genome variation data for >100k individuals). 

### The `PLINK` format
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> [`PLINK`](https://www.cog-genomics.org/plink/1.9/) is a set of utilities that allows converting VCF files into more practical formats and manipulating variant calls. It also performs many operations on variant calls, such as calculating basic statistics or linkage desiquilibrium, for example.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The `PLINK` online manual is extremely detailed but also not easy to follow. Alternatively, you may want to have a look (not now though) at a [`PLINK` tutorial](http://zzz.bwh.harvard.edu/plink/tutorial.shtml) from Harvard University or a recent [book chapter](https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3) by Christopher Chang.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Although `PLINK` can generate many different [file outputs](https://www.cog-genomics.org/plink/1.9/formats), the default outputs are as follows:
* `.bed`: binary file that contains genotype information.
* `.bim`: tab-delimited text file that always accompanies a `.bed` genotype file. It contains variant information, has no header line, and one line per variant with the following six fields:
  * Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
  * Variant identifier
  * Position in morgans or centimorgans (safe to use dummy value of '0')
  * Base-pair coordinate (1-based)
  * Allele 1 (usually minor)
  * Allele 2 (usually major)
* `.fam`: tab-delimited text file that always accompanies a `.bed` genotype file. It contains sample information, has no header line, and one line per sample with the following six fields:
  * Family ID ('FID')
  * Within-family ID ('IID'; cannot be '0')
  * Within-family ID of father ('0' if father isn't in dataset)
  * Within-family ID of mother ('0' if mother isn't in dataset)
  * Sex code ('1' = male, '2' = female, '0' = unknown)
  * Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> OPTIONAL: Install `PLINK` if not already installed.  

```bash
conda install -c bioconda plink
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Convert the VCF file to PLINK files.  

```bash
plink \
  --vcf 1kGP_chr22.vcf.gz \
  --out plink_temp
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*
- 11) How many files have been generated, and what are their extensions?
- 12) How many variants are stored in the variant file? How does it compare with the number of variants in the VCF file?  
- 13) If you look at the content of the `PLINK` variant file, you will notice that some variants are not bi-allelic SNPs. Provide an example of at most 2 other types of variations (tell what variations you observe and report the whole line for each example).  
- 14) Is the information stored in the panel file (`integrated_call_samples_v3.20130502.ALL.panel`) downloaded from the 1kGP FTP site reported in the `PLINK` sample file?
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The VCF file does not contain information about each sample's population of origin or sex. That information is stored in the panel file. Thus we need to build a file that will be used to update the `.fam` output when we convert the VCF file into `PLINK` files. For this, we have to follow instructions from the [`PLINK` online manual](http://www.cog-genomics.org/plink/1.9/data#update_indiv) to build the input file. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> By default, `PLINK` assigns the sample ID from the VCF file as both FID and IID. What we want instead is keep track of the population (pop) and super-population (super_pop) information stored in the panel file for future analyses. We will simply store `"pop"_"super_pop"` in the FID field, and store sample ID in the IID field. Also by default, `PLINK` assigns missing sex information (`0`) to all samples, unless you provide it. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Create a tab-delimited updated ID file with one line per sample and the following 5 fields:
* Old family ID (VCF sample ID)
* Old within-family ID (VCF sample ID)
* New family ID (`"pop"_"super_pop"`)
* New within-family ID (VCF sample ID)
* sex information (1 or M = male, 2 or F = female, 0 = missing) 

```bash
#Check that the panel file only contains "male" or "female" in the sex field, and does not have missing sex information (total should be 2504)  
tail -n+2 integrated_call_samples_v3.20130502.ALL.panel | cut -f4 | sort | uniq -c
```

# Generate updateFields file containing the 5 fields described above
```bash
awk -v \
 'OFS=\t' \
 'NR>1 {print $1, $1, $3"_"$2, $1, toupper(substr($4, 1, 1))}' \
 integrated_call_samples_v3.20130502.ALL.panel \
 > updateFields
```
# Check that the updateFields file contains 2504 lines
`wc -l updateFields`


<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Now we have everything to create the PLINK files. We also want to weed out variants that will not be useful for population genomics analyses, so we will keep bi-allelic SNPs only (`--snps-only just-acgt --biallelic-only strict`), keep high quality variant calls only (`--vcf-filter`), and discard rare alleles (`--maf 0.10`). Note that `PLINK` automatically re-orders alleles in minor/major. If you want to preserve the order of REF and ALT alleles as in the VCF file, then use `--keep-allele-order`.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Convert the VCF file into `PLINK` files.
```bash
# Update sex first (sex and sample IDs cannot be updated at the same time)
plink \
  --vcf 1kGP_chr22.vcf.gz \
  --snps-only just-acgt \
  --biallelic-only strict \
  --vcf-filter \
  --maf 0.10 \
  --update-sex updateFields 3 \
  --make-bed \
  --out 1kGP_chr22

# Remove the .nosex file
rm 1kGP_chr22.nosex

# Then update sample IDs in the .fam file
plink \
  --bfile 1kGP_chr22 \
  --update-ids updateFields \
  --make-just-fam \
  --out 1kGP_chr22
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*
- 15) Does the `.fam` file contain updated information? What fields have been updated when compared to `plink_temp.fam`?
- 16) How many variants are stored in the `.bim` file? How does it compare with the number of variants in `plink_temp.bim`?
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Some population genomics analyses that focus on population demographic history and structure perform better if variants are in relative genome-wide linkage equilibrium, meaning that alleles at different SNP loci must be randomly associated. Indeed, non-random association between alleles (a.k.a. linkage disequilibrium, or LD) would mean redundancy in the data, which would increase computing unnecessarily. For other applications that focus on genomic regions (e.g., natural selection), loci in LD are highly informative. `plink --indep-pairwise` calculates the square of the correlation (*r*<sup>2</sup>) between allele counts in adjacent SNP loci and stores loci that are below (`.prune.in` output file) and above (`.prune.out` output file) a user-defined threshold. *r*<sup>2</sup>=0 when two loci are in perfect equilibrium, *r*<sup>2</sup>=1 when two loci provide redundant information.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Calculate pairwise *r*<sup>2</sup> and create lists of SNP loci in LD or not.
```bash
plink \
  --bfile 1kGP_chr22 \
  --indep-pairwise 200kb 1 0.5 \
  --out ld_snps
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*
- 17) How many variants in the `.prune.in` and `.prune.out` output files?
- 18) How does it compare to the number of variants in `1kGP_chr22.bim`?
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> You can now build `PLINK` files with just the LD-pruned data.
```bash
plink \
 --bfile 1kGP_chr22 \
 --extract ld_snps.prune.in \
 --make-bed \
 --out 1kGP_chr22.ldpruned
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*
- 19) In terms of file size, what do you notice when you look at the `.bed`, `.bim` and `.fam` files before and after LD pruning?
- 20) How do you explain the changes, or lack thereof?
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Let's see how the non-LD-pruned and LD-pruned data behave in a PCA plot.
```bash
# PCA on non-LD-pruned data
plink \
 --bfile 1kGP_chr22 \
 --pca 5 \
 --out 1kGP_chr22.pca_results

# PCA on LD-pruned data
plink \
 --bfile 1kGP_chr22.ldpruned \
 --pca 5 \
 --out 1kGP_chr22.ldpruned.pca_results
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> `PLINK` PCA has generated two outputs with suffixes `.eigenvec` (the PC coordinates for each sample) and `.eigenval` (all the eigenvalues). Go to the `R` console and create screeplots and PCA plots.
```R
library(tidyr)
library(ggplot2)
library(cowplot)

# Non-LD-pruned data for scree plot
adat.scree <- read.table("1kGP_chr22.pca_results.eigenval", header = FALSE)
# Add a column with row number (only needed to be able to do a bar plot)
adat.scree$Name = 1:nrow(adat.scree)
# Rename columns
colnames(adat.scree) <- c("Scree","Name")
# Do the same steps for the LD-pruned data
bdat.scree <- read.table("1kGP_chr22.ldpruned.pca_results.eigenval", header = FALSE)
bdat.scree$Name = 1:nrow(bdat.scree)
colnames(bdat.scree) <- c("Scree","Name")
# Plot non-LD-pruned data
adat.screep <- ggplot(adat.scree,aes(Name,Scree)) +
               geom_bar(stat="identity") + # heights of the bars represent values in the data
               theme(text = element_text(size = 20)) +
               theme(axis.text.x=element_blank(), # no text for the x-axis
                     axis.ticks.x=element_blank()) + # no ticks for the x-axis
               xlab("component") +
               ylab("eigenvalue") +
               ggtitle("non-LD-pruned scree plot")
# Do the same steps for the LD-pruned data
bdat.screep <-ggplot(bdat.scree,aes(Name,Scree)) +
              geom_bar(stat="identity") +
              theme(text = element_text(size = 20)) +
              theme(axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()) +
              xlab("component") +
              ylab("eigenvalue") +
              ggtitle("LD-pruned scree plot")
# Combine scree plots using plot_grid from cowplot
plot_grid(adat.screep,
          bdat.screep,
          align = 'vh',
          hjust = -1,
          nrow = 1)

# Create dataframe for PCA for the non-LD-pruned data
adat <- read.table("1kGP_chr22.pca_results.eigenvec", header = FALSE)
# Rename columns
colnames(adat) <- c("POP", "SAMPLE", "PC1", "PC2", "PC3", "PC4", "PC5")
# Split POP column into super-population SUPERPOP and population POP
adat <- separate(data = adat, col = POP, into = c("SUPERPOP", "POP"), sep = "_")
# Create plot with populations in different colours and super-populations with different point shapes
adat.pc12 <- ggplot(adat, aes(x = PC1, y = PC2, colour = POP, shape = SUPERPOP)) + 
             geom_point() +
             ggtitle("Non-LD-pruned")
# Do the same steps for the LD-pruned data
bdat <- read.table("1kGP_chr22.ldpruned.pca_results.eigenvec", header = FALSE)
colnames(bdat) <- c("POP", "SAMPLE", "PC1", "PC2", "PC3", "PC4", "PC5")
bdat <- separate(data = bdat, col = POP, into = c("SUPERPOP", "POP"), sep = "_")
bdat.pc12 <- ggplot(bdat, aes(x = PC1, y = PC2, colour = POP, shape = SUPERPOP)) + 
             geom_point() +
             ggtitle("LD-pruned")
# Combine plots using plot_grid from cowplot
prow <- plot_grid(adat.pc12 + theme(legend.position="none"), # remove the legend
                  bdat.pc12 + theme(legend.position="none"), # remove the legend
                  align = 'vh', # plots are aligned vertically and horizontally
                  nrow = 1) # 2 plots on one row
# Prepare common legend for the two plots
legend <- get_legend(adat.pc12 + 
                     guides(color = guide_legend(nrow = 4)) + # legend spans 4 lines
                     theme(legend.position = "bottom")) # legend is displayed at the bottom of the plots (i.e., horizontal not vertical)
# Combine plots and legend using plot_grid from cowplot
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .3)) # ratio between plots and legend is 1:0.3
```

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*
- 21) Do you observe any obvious differences between the two plots?
- 22) What patterns do you observe?
---


### The Eigensoft format
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> PLINK was initially developed for GWAS studies and similar largescale medical genomics studies. Another suite of utilities ([`Eigensoft`](https://github.com/DReichLab/EIG)) was developed for population genomics, and as is often the case, the file formats remained different between the two suites of utilities. However, Eigensoft can convert `PLINK` files into many file formats (including `EIGENSTRAT` files that we will use in this tutorial) using [`CONVERTF`](https://github.com/DReichLab/EIG/tree/master/CONVERTF).

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The `EIGENSTRAT` files contain more or less the same information as the `PLINK` files, just in a different format:
* `.eigenstratgeno`: tab-delimited genotype file with one line per SNP and and genotypes in non-separated columns, with the following genotype coding:
  * `0`: no copies of reference allele
  * `1`: one copy of reference allele
  * `2`: two copies of reference allele
  * `9`: missing data
* `.snp`: tab-delimited SNP file with one line per SNP and the following 6 columns (last 2 optional):
  * SNP name
  * Chromosome (X is encoded as 23, Y as 24, mtDNA as 90, and XY as 91)
  * Genetic position (in Morgans). 0 if unknown
  * Physical position (in bases)
  * Optional 5th and 6th columns are reference and variant alleles. For monomorphic SNPs, the variant allele can be encoded as X (unknown)
* `.ind`: tab-delimited sample file with one line per individual and the following 3 columns:
  * sample ID
  * gender (M or F). U for Unknown
  * Case or Control status, or population group label. If this entry is set to "Ignore", then that individual and all genotype data from that individual will be removed from the data set in all `CONVERTF` output.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Build parameter files that will be the inputs for `CONVERTF`. The content of the parameter files is as follows:
* `par.PACKEDPED.EIGENSTRAT.1kGP_chr22`:
```bash
genotypename:    1kGP_chr22.bed
snpname:         1kGP_chr22.bim
indivname:       1kGP_chr22.fam
outputformat:    EIGENSTRAT
genotypeoutname: 1kGP_chr22.eigenstratgeno
snpoutname:      1kGP_chr22.snp
indivoutname:    1kGP_chr22.ind
```
* `par.PACKEDPED.EIGENSTRAT.1kGP_chr22.ldpruned`:
```bash
genotypename:    1kGP_chr22.ldpruned.bed
snpname:         1kGP_chr22.ldpruned.bim
indivname:       1kGP_chr22.ldpruned.fam
outputformat:    EIGENSTRAT
genotypeoutname: 1kGP_chr22.ldpruned.eigenstratgeno
snpoutname:      1kGP_chr22.ldpruned.snp
indivoutname:    1kGP_chr22.ldpruned.ind
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Run `CONVERTF`.
```bash
convertf -p par.PACKEDPED.EIGENSTRAT.1kGP_chr22
convertf -p par.PACKEDPED.EIGENSTRAT.1kGP_chr22.ldpruned
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Build parameter files that will be the inputs for [`SMARTPCA`](https://github.com/DReichLab/EIG/tree/master/POPGEN). You need to build new parameter files as follows:
* `par.1kGP_chr22`:
```bash
genotypename:    1kGP_chr22.eigenstratgeno
snpname:         1kGP_chr22.snp
indivname:       1kGP_chr22.ind
evecoutname:     1kGP_chr22.smartpca_results.evec
evaloutname:     1kGP_chr22.smartpca_results.eval
numoutevec:      5
```
* `par.1kGP_chr22.ldpruned`:
```bash
genotypename:    1kGP_chr22.ldpruned.eigenstratgeno
snpname:         1kGP_chr22.ldpruned.snp
indivname:       1kGP_chr22.ldpruned.ind
evecoutname:     1kGP_chr22.ldpruned.smartpca_results.evec
evaloutname:     1kGP_chr22.ldpruned.smartpca_results.eval
numoutevec:      5
```
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Run `SMARTPCA`.
```bash
smartpca -p par.1kGP_chr22
smartpca -p par.1kGP_chr22.ldpruned
```
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> `SMARTPCA` has generated two output files with the suffixes `.evec` (first row is the eigenvalues for the first 5 PCs, and all further rows contain the PC coordinates for each sample) and `.evac` (all the eigenvalues). Go to the `R` console and create PCA plots.  
```R
library(tidyr)
library(ggplot2)
library(cowplot)

#Non-LD-pruned data for scree plot  
adat.scree <- read.table("1kGP_chr22.smartpca_results.eval", header = FALSE)  
#Add a column with row number (only needed to be able to do a bar plot)  
adat.scree$Name = 1:nrow(adat.scree)  
#Rename columns  
colnames(adat.scree) <- c("Scree","Name")
#Restrict to the first 10 lines (first 10 eigenvalues)  
adat.scree <- head(adat.scree, 10)  
#Do the same steps for the LD-pruned data  
bdat.scree <- read.table("1kGP_chr22.ldpruned.smartpca_results.eval", header = FALSE)  
bdat.scree$Name = 1:nrow(bdat.scree)  
colnames(bdat.scree) <- c("Scree", "Name")  
bdat.scree <- head(bdat.scree, 10)  
#Plot non-LD-pruned data  
adat.screep <- ggplot(adat.scree,aes(x = Name, y = Scree)) +  
               geom_bar(stat = "identity") + # heights of the bars represent values in the data  
               theme(text = element_text(size = 20)) +  
               theme(axis.text.x = element_blank(), # no text for the x-axis  
                     axis.ticks.x = element_blank()) + # no ticks for the x-axis  
               xlab("component") + # x-axis label  
               ylab("eigenvalue") + # y-axis label  
               ggtitle("non-LD-pruned scree plot")  
#Do the same steps for the LD-pruned data  
bdat.screep <-ggplot(bdat.scree,aes(x = Name, y = Scree)) +  
              geom_bar(stat = "identity") +  
              theme(text = element_text(size = 20)) +  
              theme(axis.text.x = element_blank(),  
                    axis.ticks.x = element_blank()) +  
              xlab("component") +  
              ylab("eigenvalue") +  
              ggtitle("LD-pruned scree plot")  
#Combine scree plots using plot_grid from cowplot  
plot_grid(adat.screep,  
          bdat.screep,  
          align = 'vh', # plots are aligned vertically and horizontally  
          nrow = 1) # 2 plots on one row  
          
#Create dataframe for the non-LD-pruned data  
adat <- read.table("1kGP_chr22.smartpca_results.evec", header = FALSE)  
#Reduce table to just the sample name and the first 3 PCs  
adat <- adat[,c(1:4)]
#Rename columns  
colnames(adat) <- c("POPSAMPLE", "PC1", "PC2", "PC3")  
#Split POPSAMPLE column into super-population SUPERPOP, population POP, and SAMPLE  
adat <- separate(data = adat, col = POPSAMPLE, into = c("SUPERPOP", "POP", "SAMPLE"), sep = "_|:")  
#Create plot with populations in different colours and super-populations with different point shapes  
adat.pc12 <- ggplot(adat, aes(x = PC1, y = PC2, colour = POP, shape = SUPERPOP)) +   
             geom_point() +  
             ggtitle("Non-LD-pruned")  
#Do the same steps for the LD-pruned data  
bdat <- read.table("1kGP_chr22.ldpruned.smartpca_results.evec", header = FALSE)  
bdat <- bdat[,c(1:4)]  
colnames(bdat) <- c("POPSAMPLE", "PC1", "PC2", "PC3")  
bdat <- separate(data = bdat, col = POPSAMPLE, into = c("SUPERPOP", "POP", "SAMPLE"), sep = "_|:")  
bdat.pc12 <- ggplot(bdat, aes(x = PC1, y = PC2, colour = POP, shape = SUPERPOP)) +   
             geom_point() +  
             ggtitle("LD-pruned")  
#Combine plots using plot_grid from cowplot  
prow <- plot_grid(adat.pc12 + theme(legend.position="none"), # remove the legend  
                  bdat.pc12 + theme(legend.position="none"), # remove the legend  
                  align = 'vh', # plots are aligned vertically and horizontally  
                  nrow = 1) # 2 plots on one row  
#Prepare common legend for the two plots  
legend <- get_legend(adat.pc12 +   
                     guides(color = guide_legend(nrow = 4)) + # legend spans 4 lines  
                     theme(legend.position = "bottom")) # legend is displayed at the bottom of the plots (i.e., horizontal not vertical)
#Combine plots and legend using plot_grid from cowplot  
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, 0.3)) # ratio between plots and legend is 1:0.3  

```

---  

#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*   

- 23) Are the `SMARTPCA` results fundamentally different from `PLINK` PCA results?
---
