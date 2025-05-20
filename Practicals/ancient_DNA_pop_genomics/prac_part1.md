# Ancient DNA and population genomics practical: Part 1 - Bastien Llamas

---
The two practicals will be separated into:
1. Data handling
2. Population genomics applications

Icons are used to highlight sections of the practicals:

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Information  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Hands-on tasks  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> Questions  

---
# Day 1: Data handling

## Practical outcomes

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> At the end of today's practical, you will know how to convert the information contained in a VCF file into other formats compatible with widely-used population genomics programs.

## Population genomics
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> In a nutshell, population genomics is the study of the genomic composition of populations and how evolution shaped it. Questions in population genomics classically target demographic history (population size through time), gene flow between populations, populations ancestry, or identification of conserved biology units.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> While population genetics is usually restricted to a small set of genetic loci, population genomics leverages the large genomic datasets that have become available in recent years and uses up to millions of genetic loci at once.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> We are not going to focus on the mathematical aspects of population genomics, but rather on how to manipulate genomic datasets and learn about a few popular approaches and tools. I encourage you to read Graham Coop's [course notes](https://github.com/cooplab/popgen-notes/blob/master/popgen_notes.pdf) if you are curious about the underlying mathematical theories.


## VCF format: a reminder
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> You have previously learnt about several raw or processed high throughput sequencing data formats (e.g., FASTQ, SAM/BAM, VCF). In particular, you should now know that [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) files contain information about variant calls found at specific positions in a reference genome.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We will use a VCF file of human chromosome 22 from the 1000 Genomes Project (1kGP) that we will save into a working directory in your home directory:

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Create working directories:
```bash
mkdir -p ~/Project_12_1/{data,scripts,results}
cd ~/Project_12_1/
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Download compressed VCF file and its index from the 1kGP public FTP site (VCF file size: 214453750 bytes):
```bash
wget -O data/1kGP_chr22.vcf.gz 'ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
wget -O data/1kGP_chr22.vcf.gz.tbi 'ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi'
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Although you could use your own scripts to parse VCF files and analyse variant calls, several tools have already been developed for your convenience. We will be using the following tools for the next two practicals
- [`bcftools`](http://samtools.github.io/bcftools/bcftools.html) is a set of utilities to manipulate variant calls in VCF files. 
- [`plink`](https://www.cog-genomics.org/plink/) is a tool kit for population genomics and genome-wide association studies (GWAS).
- [EIGENSOFT](https://github.com/DReichLab/EIG) is a set of tools for population genomic analysis such as principal component analysis (PCA).
- [AdmixTools](https://github.com/DReichLab/AdmixTools) is used to investigate admixture and population history through ancestry estimation and admixture graph modelling.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We will also be using some custom R scripts to visualise the data. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/>  Activate the `popgen` environment and copy R scripts:
```bash
source activate popgen
cp ~/data/ancient/prac_1/* ~/Project_12_1/scripts/
```

#### VCF meta-information and header lines
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the compressed VCF file using `zless`. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Meta-information lines start with `##` and contain various metadata. The header line starts with `#` and is tab separated. It contains 9 columns of information about the variant calls, and then one column per sample:

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
wget --directory-prefix data 'ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/README_vcf_info_annotation.20141104'
```

#### VCF body
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The body of the VCF file is tab separated. Each line represents a unique variant site.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Before we move forward, let's see if you can retrieve basic information from a 1kGP VCF file that will be useful for population genomic analysis.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*
1. Determine how many variant sites are recorded in the VCF file. You can use `bcftools stats`, or `bcftools view` and bash commands.
2. Determine how many samples are recorded in the VCF file. You can use `bcftools stats`, or `bcftools query` and bash commands.
3. The `INFO` fields contain a lot of information. In particular for the first variant position in the file: determine how many samples have data, how many ALT alleles are reported,  what the frequency of the ALT allele is globally, and what the frequency of the ALT allele is in South Asians.
4. Same as question 3 for variant position 16051249 (see the [BCFtools manual](http://samtools.github.io/bcftools/bcftools.html) for region or target formatting).
5. How many alternative alleles are observed at position 16050654?
6. Looking at the information contained in the `FORMAT` field in the body of the VCF file, what kind of data is stored in the VCF file for each sample?

<details>
  <summary>Answers</summary>
  
  Q1: 1,103,547 variants<br>
  `bcftools stats data/1kGP_chr22.vcf.gz | less`<br>
  `bcftools view -H data/1kGP_chr22.vcf.gz | wc -l`<br>

  Q2: 2,504 samples<br>
  `bcftools stats data/1kGP_chr22.vcf.gz | less`<br>
  `bcftools query -l data/1kGP_chr22.vcf.gz | wc -l`<br>

  Q3: AC=1, AF=0.000199681, SAS_AF=0.001<br>
  ```bash
  bcftools view -H data/1kGP_chr22.vcf.gz | head -n1 | awk '{print $1,$2,$8;}'
  ```

  Q4: AC=563, AF=0.11242, SAS_AF=0.2791<br>
  ```bash
  bcftools view -H data/1kGP_chr22.vcf.gz 22:16051249 | awk '{print $1,$2,$8;}'
  ```

  Q5: AC=9,87,599,20 so 4 alleles<br>
  ```bash
  bcftools view -H data/1kGP_chr22.vcf.gz 22:16050654 | awk '{print $1,$2,$8;}'
  ```

  Q6: GT, i.e. genotype<br>
  ```bash
  bcftools view -h data/1kGP_chr22.vcf.gz
  ```
</details>

---

#### Other useful 1kGP metadata
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Download sample details from the 1kGP FTP site to learn about population of origin and sex of each individual.
```bash
wget --directory-prefix data 'ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'
```
It is a tabulated file, with 4 columns corresponding to individual IDs, populations, super-populations, and sex.
 
---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*  

7. Using bash commands on the panel file you just downloaded, determine how many different populations and super-populations are represented in the 1kGP dataset.
8. How many individuals are in each super-population?  

<details>
  <summary>Answers</summary>

  Q7: 26 populations and 5 super-populations
  ```bash
  tail -n+2 data/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $2;}' | sort | uniq | wc -l
  ```
  ```bash
  tail -n+2 data/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3;}' | sort | uniq | wc -l
  ```

  Q8: 661 AFR, 347 AMR, 504 EAS, 503 EUR, 489 SAS
  ```bash
  tail -n+2 data/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3;}' | sort | uniq -c
  ```
</details>

---


<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> You can learn more about the populations in the 1kGP [here](https://www.internationalgenome.org/faq/which-populations-are-part-your-study/).  

## Converting VCF files into population genomics formats
### Rationale
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> A VCF file may contain a lot of information (e.g. variant annotation) that can be very useful for clinical genomics. However, population genomics applications only need a subset of the information in a VCF file, i.e. variant genomic coordinates, variant ID, reference (`REF`) and alternative (`ALT`) alleles, and sample genotypes. This is what you typically find in the 1kGP data.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Genotypes can be coded differently depending on ploidy, the number of alternate alleles, whether the genomes are phased (i.e. alleles on maternal and paternal chromosomes are identified) or not, and homo/heterozygosity. For convenience, 0 is used to code `REF`, and 1, 2, 3, etc are used to code `ALT`.  

||Example|  
|:-|:-|  
|**Haploid, 2 alleles**|2 genotypes: 0, 1|  
|**Diploid, not phased, 2 alleles**|3 genotypes: 0/0, 0/1, 1/1|  
|**Diploid, phased, 2 alleles**|4 genotypes: 0\|0, 0\|1, 1\|0, 1\|1|  
|**Haploid, n alleles**|n genotypes: e.g., 0, 1, 2, 7|  
|**Diploid, not phased, n alleles**|n! genotypes: e.g., 0/0, 1/4, 0/2, 3/3|  
|**Diploid, phased, n alleles**|n<sup>2</sup> genotypes: e.g., 0\|0, 1\|4, 0\|2, 3\|3|  

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the variant in position 16061250 in the VCF file.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*  

9. What are the `REF` and `ALT` alleles?
10. Given `REF` and `ALT` alleles found when answering question 9, and knowing that the genotypes are phased, what are all the possible genotypes?  

<details>
  <summary>Answers</summary>

  Q9: `REF` = T and `ALT` = A,C
  ```bash
  bcftools view -H data/1kGP_chr22.vcf.gz 22:16061250 | awk '{print $1,$2,$4,$5}'
  ```

  Q10: 0|0, 0|1, 1|0, 1|1, 0|2, 2|0, 1|2, 2|1, 2|2
  
</details>
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

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Convert the VCF file to PLINK files.  

```bash
plink \
  --vcf data/1kGP_chr22.vcf.gz \
  --out results/plink_temp
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*

11. How many files have been generated, and what are their extensions?
12. How many variants are stored in the variant file? How does it compare with the number of variants in the VCF file?
13. If you look at the content of the `PLINK` variant file, you will notice that some variants are not bi-allelic SNPs. Provide an example of at most 2 other types of variations (tell what variations you observe and report the whole line for each example).
14. Is the information stored in the panel file (`integrated_call_samples_v3.20130502.ALL.panel`) downloaded from the 1kGP FTP site reported in the `PLINK` sample file? *Hint: look at the `.fam` file*

<details>
  <summary>Answers</summary>

  Q11: 4 files with extensions `.bed`, `.bim`, `.fam`, `.log`
 
  Q12: 1,103,547 variants (same number of variants as in VCF file)
  ```bash
  wc -l results/plink_temp.bim
  ```

  Q13: Multi allelic variation (CNV): `22      esv3647175;esv3647176;esv3647177;esv3647178     0       16050654        <CN3>   A`\
  Indel: `22      rs587747231     0       16050739        T       TA`

  Q14: No
  
</details>
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The VCF file does not contain information about each sample's population of origin or sex. That information is stored in the panel file. Thus we need to build a file that will be used to update the `.fam` output when we convert the VCF file into `PLINK` files. For this, we have to follow instructions from the [`PLINK` online manual](http://www.cog-genomics.org/plink/1.9/data#update_indiv) to build the input file. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> By default, `PLINK` assigns the sample ID from the VCF file as both FID and IID. What we want instead is keep track of the population (`pop`) and super-population (`super-pop`) information stored in the panel file for future analyses. We will simply store `pop_super-pop` in the FID field, and store sample ID in the IID field. Also by default, `PLINK` assigns missing sex information (`0`) to all samples, unless you provide it. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Create a tab-delimited updated ID file with one line per sample and the following 5 fields:
* Old family ID (VCF sample ID)
* Old within-family ID (VCF sample ID)
* New family ID (`pop_super-pop`)
* New within-family ID (VCF sample ID)
* sex information (`1` or `M` = male, `2` or `F` = female, `0` = missing) 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Check that the panel file only contains "male" or "female" in the sex field, and does not have missing sex information (total should be 2,504):
```bash
tail -n+2 data/integrated_call_samples_v3.20130502.ALL.panel | cut -f4 | sort | uniq -c
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Generate updateFields file containing the 5 fields described above:
```bash
awk -v \
 'OFS=\t' \
 'NR>1 {print $1, $1, $3"_"$2, $1, toupper(substr($4, 1, 1))}' \
 data/integrated_call_samples_v3.20130502.ALL.panel \
 > data/updateFields
```
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Check that the updateFields file contains 2,504 lines:
```bash
wc -l data/updateFields
```


<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Now we have everything to create the `PLINK` files. We also want to weed out variants that will not be useful for population genomics analyses, so we will keep bi-allelic SNPs only (`--snps-only just-acgt --biallelic-only strict`), keep high quality variant calls only (`--vcf-filter`), and discard rare alleles (i.e. alleles with a minor allele frequency < 10%; `--maf 0.10`). Note that `PLINK` automatically re-orders alleles in minor/major. You may want to preserve the order of REF and ALT alleles as in the VCF file, so use `--keep-allele-order`.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> To convert the VCF file into `PLINK` files, you need to update sex first (sex and sample IDs cannot be updated at the same time).  

```bash
plink \
  --vcf data/1kGP_chr22.vcf.gz \
  --snps-only just-acgt \
  --biallelic-only strict \
  --vcf-filter \
  --maf 0.10 \
  --update-sex data/updateFields 3 \
  --make-bed \
  --out results/1kGP_chr22
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> You can delete the `.nosex` file.  

```bash
rm results/1kGP_chr22.nosex
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Then update sample IDs in the `.fam` file.  

```bash
plink \
  --bfile results/1kGP_chr22 \
  --update-ids data/updateFields \
  --make-just-fam \
  --out results/1kGP_chr22
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*   

15. Does the `.fam` file contain updated information? What fields have been updated when compared to `plink_temp.fam`?
16. How many variants are stored in the `.bim` file? How does it compare with the number of variants in `plink_temp.bim`?  

<details>
  <summary>Answers</summary>

  Q15: Yes, fields 1 and 5
 
  Q16: 73,246 variants
  ```bash
  wc -l results/1kGP_chr22.bim
  ```
  ```bash
  wc -l results/plink_temp.bim
  ```

</details>
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Some population genomics analyses that focus on population demographic history and structure perform better if variants are in relative genome-wide linkage equilibrium, meaning that alleles at different SNP loci must not be linked, or in other words must be randomly associated. Indeed, non-random association between alleles (a.k.a. linkage disequilibrium, or LD) would mean redundancy in the data, which would increase computing unnecessarily. For other applications that focus on genomic regions (e.g. natural selection), loci in LD are highly informative. `plink --indep-pairwise` calculates the square of the correlation (*r*<sup>2</sup>) between allele counts in adjacent SNP loci and stores loci that are below (`.prune.in` output file) and above (`.prune.out` output file) a user-defined threshold. *r*<sup>2</sup>=0 when two loci are in perfect equilibrium, *r*<sup>2</sup>=1 when two loci provide redundant information. Here we will use a threshold of 0.5.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Calculate pairwise *r*<sup>2</sup> and create lists of SNP loci in LD or not.  

```bash
plink \
  --bfile results/1kGP_chr22 \
  --indep-pairwise 200kb 1 0.5 \
  --out results/ld_snps
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*   

17. How many variants in the `.prune.in` and `.prune.out` output files?
18. How does it compare to the number of variants in `1kGP_chr22.bim`?  

<details>
  <summary>Answers</summary>

  Q17: 12,142 variants in `.prune.in` and 61,104 variants in `.prune.out`
  ```bash
  wc -l results/ld_snps.prune.in
  ```
  ```bash
  wc -l results/ld_snps.prune.out
  ```
 
  Q18: The sum of the number of variants in `.prune.in` and `.prune.out` is the total number of variants in `1kGP_chr22.bim`

</details>
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> You can now build `PLINK` files with just the LD-pruned data.  

```bash  
plink \
 --bfile results/1kGP_chr22 \
 --extract results/ld_snps.prune.in \
 --make-bed \
 --out results/1kGP_chr22.ldpruned
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Have a look at the newly created files.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*  

19. In terms of file size, what do you notice when you look at the `.bed`, `.bim` and `.fam` files before and after LD pruning?
20. How do you explain the changes, or lack thereof?  

<details>
  <summary>Answers</summary>

  Q19: smaller `.bed` and `.bim` after LD pruning, no change for `.fam`
 
  Q20: some variants have been pruned but all samples are kept

</details>
---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Let's see how the non-LD-pruned and LD-pruned data behave in a PCA plot. We will consider only 5 eigenvectors in this analysis.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> PCA on non-LD-pruned data:

```bash
plink \
 --bfile results/1kGP_chr22 \
 --pca 5 \
 --out results/1kGP_chr22.pca_results
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> PCA on LD-pruned data:
```bash
plink \
 --bfile results/1kGP_chr22.ldpruned \
 --pca 5 \
 --out results/1kGP_chr22.ldpruned.pca_results
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> `PLINK` PCA has generated two outputs with suffixes `.eigenvec` (the PC coordinates for each sample) and `.eigenval` (all the eigenvalues). Run the custom Rscript to generate screeplots and PCA plots. The plots should be in `.pdf` files in `results/`.

```bash
Rscript scripts/plot_plink_pca.R
```

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*  

21. Do you observe any obvious differences between the two PCA plots?
22. What patterns do you observe?  

<details>
  <summary>Answers</summary>

  Q21: clusters seem more diffused with the non-LD-pruned data
 
  Q22: relatively obvious clustering by super-populations, AMR seem admixed between EUR, AFR and EAS
</details>
---


### The Eigensoft format
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> `PLINK` was initially developed for GWAS and similar largescale medical genomics studies. Another suite of utilities, ([`Eigensoft`](https://github.com/DReichLab/EIG)), was developed for population genomics, and as is often the case, the file formats are different between the two suites of utilities. However, `Eigensoft` can convert `PLINK` files into many file formats (including `EIGENSTRAT` files that we will use in this practical) using [`CONVERTF`](https://github.com/DReichLab/EIG/tree/master/CONVERTF).

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The `EIGENSTRAT` files contain more or less the same information as the `PLINK` files, just in a different format:  

* `.eigenstratgeno`: tab-delimited genotype file with one line per SNP and and genotypes in non-separated columns, with the following genotype coding:
  * `0`: no copies of reference allele
  * `1`: one copy of reference allele
  * `2`: two copies of reference allele
  * `9`: missing data
* `.snp`: tab-delimited SNP file with one line per SNP and the following 6 columns (last 2 optional):
  * SNP name
  * Chromosome (X is encoded as `23`, Y as `24`, mtDNA as `90`, and XY as `91`)
  * Genetic position (in Morgans). `0` if unknown
  * Physical position (in bases)
  * Optional 5th and 6th columns are reference and variant alleles. For monomorphic SNPs, the variant allele can be encoded as `X` (unknown)
* `.ind`: tab-delimited sample file with one line per individual and the following 3 columns:
  * sample ID
  * gender (`M` or `F`). `U` for Unknown
  * Case or Control status, or population group label. If this entry is set to "Ignore", then that individual and all genotype data from that individual will be removed from the data set in all `CONVERTF` output.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Usually we need to build parameter files that will be the inputs for `CONVERTF`. The content of the parameter files would look like this:  

* `results/par.PACKEDPED.EIGENSTRAT.1kGP_chr22`:
```bash
genotypename:    results/1kGP_chr22.bed
snpname:         results/1kGP_chr22.bim
indivname:       results/1kGP_chr22.fam
outputformat:    EIGENSTRAT
genotypeoutname: results/1kGP_chr22.eigenstratgeno
snpoutname:      results/1kGP_chr22.snp
indivoutname:    results/1kGP_chr22.ind
```
* `results/par.PACKEDPED.EIGENSTRAT.1kGP_chr22.ldpruned`:
```bash
genotypename:    results/1kGP_chr22.ldpruned.bed
snpname:         results/1kGP_chr22.ldpruned.bim
indivname:       results/1kGP_chr22.ldpruned.fam
outputformat:    EIGENSTRAT
genotypeoutname: results/1kGP_chr22.ldpruned.eigenstratgeno
snpoutname:      results/1kGP_chr22.ldpruned.snp
indivoutname:    results/1kGP_chr22.ldpruned.ind
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> However, for the sake of simplicity we will run the commands like this. *Make sure to copy the whole block of code*

```bash
convertf -p <(echo "genotypename:    results/1kGP_chr22.bed
snpname:         results/1kGP_chr22.bim
indivname:       results/1kGP_chr22.fam
outputformat:    EIGENSTRAT
genotypeoutname: results/1kGP_chr22.eigenstratgeno
snpoutname:      results/1kGP_chr22.snp
indivoutname:    results/1kGP_chr22.ind")
```

```bash
convertf -p <(echo "genotypename:    results/1kGP_chr22.ldpruned.bed
snpname:         results/1kGP_chr22.ldpruned.bim
indivname:       results/1kGP_chr22.ldpruned.fam
outputformat:    EIGENSTRAT
genotypeoutname: results/1kGP_chr22.ldpruned.eigenstratgeno
snpoutname:      results/1kGP_chr22.ldpruned.snp
indivoutname:    results/1kGP_chr22.ldpruned.ind")
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Run [`SMARTPCA`](https://github.com/DReichLab/EIG/tree/master/POPGEN), which is the PCA software as part of `Eigensoft`, on the `EIGENSTRAT` files. We will consider only 5 eigenvectors in this analysis.

```bash
smartpca -p <(echo "genotypename:    results/1kGP_chr22.eigenstratgeno
snpname:         results/1kGP_chr22.snp
indivname:       results/1kGP_chr22.ind
evecoutname:     results/1kGP_chr22.smartpca_results.evec
evaloutname:     results/1kGP_chr22.smartpca_results.eval
numoutevec:      5")
```

```bash
smartpca -p <(echo "genotypename:    results/1kGP_chr22.ldpruned.eigenstratgeno
snpname:         results/1kGP_chr22.ldpruned.snp
indivname:       results/1kGP_chr22.ldpruned.ind
evecoutname:     results/1kGP_chr22.ldpruned.smartpca_results.evec
evaloutname:     results/1kGP_chr22.ldpruned.smartpca_results.eval
numoutevec:      5")
```
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> `SMARTPCA` has generated two output files with the suffixes `.evec` (first row is the eigenvalues for the 5 PCs, and all further rows contain the PC coordinates for each sample) and `.evac` (all the eigenvalues). Run the custom Rscript to generate screeplots and PCA plots. The plots should be in `.pdf` files in `results/`.

```bash
Rscript scripts/plot_smartpca.R
```

---  

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/>*Questions*   

23. Are the `SMARTPCA` results fundamentally different from `PLINK` PCA results?  

<details>
  <summary>Answers</summary>

  Q23: No
</details>

---
