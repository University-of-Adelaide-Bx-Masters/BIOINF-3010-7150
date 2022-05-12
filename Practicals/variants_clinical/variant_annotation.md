
# Clinical Genomics Practicals

Practicals for clinical genomics are split into two tutorials:

1. Variant Annotation (Tuesday)
2. Variant Filtering and Pedigree Analyses (Wednesday)

## Some backgound on Clinical Genomics

Clinical genomics and mendelian genetics can be separated into two groups based on their stated analysis goals; Diagnostics and Research.
Many labs in major hospitals, independent research institutes and universities do a mixture of both work.
They have a goal of informing clinical practice (_i.e.,_ can we inform the diagnosis of the clinician in the lab?) through genome sequencing methods, and can further investigate difficult cases that could ultimately lead to developing new tests and techniques.
On a computational level, these environments are where computational methods, workflows and standards are developed.

The complete opposite is true of diagnostic approaches.
The goal is to identify _known_ patterns in sequencing data, with an emphasis of using well established and highly tested "best practices".
Diagnostic tests are also about identifying known issues in the quickest and most reliable way possible, because treatment pathways are reliant on the information, much like other tests such as blood profiles.

A critical thing to considering about high-throughput sequencing is the enormous amount of information that is obtained.
Additionally, a lot of the information that is obtained is difficult to interpret.
This is the reason why we primarily focus on "protein-coding" regions (~1% of the total sequence), where we know there is a good chance that a genetic variation may bring about a phenotypic change.
Non-coding regions of the genome, or the other 99% are much more difficult to interpret, although recent experimental projects such as the [Epigenomics Roadmap](http://www.roadmapepigenomics.org/) or [Encyclopedia of DNA elements (ENCODE)](https://www.encodeproject.org/) are trying to change that, by identifying "functional" regions of the non-coding genome that impact gene expression.

In this tutorial, we are firstly going to look at ways in which we can give context to variant analyses, by running a process called "Variant Annotation".

## 1. Variant Annotation

### This week's data

This week we are going to use three samples from a [gemini database tutorial](https://s3.amazonaws.com/gemini-tutorials/Gemini-Recessive-Tutorial.pdf) written by Aaron Quinlan (University of Utah).
Aaron's group has written some really helpful pieces of software including `bedtools`, `giggle` and `vcfanno`, which are all becoming standard tools in the toolkit of clinical researchers.

To demonstrate a lot of what we're doing this week, we'll be using variants that have already been called using the GenomeAnalysisToolkit (GATK), that was developed by the Broad Institute (Cambridge, USA).
The data was pre-loaded to your VM at ***/home/student/data/Variants/***:
- 1 gzipped compressed Variant Call Format (VCF) = _trio.trim.vep.vcf.gz_
- 2 pedigree files that contain two separate examples looking at dominant and recessive disorders = _recessive.ped / dominant.ped_
- 1 SNP ID annotation file from dbSNP = _hg19.dbSNP.vcf.gz_

We will mostly use **BCFtools** for this tutorial:
https://samtools.github.io/bcftools/bcftools.html

As our VCF file already has annotations attached, let's start by stripping off that information so we can start the process at the start (and hopefully learn a few things along the way!).

```bash
# First let's make a working directory and copy the data over
mkdir -p ~/clinicalGenomics && cd $_
cp /home/student/data/Variants/*  ~/clinicalGenomics/

# Then make a standalone conda environment and install BCFtools through bioconda
conda create --name bcftoolsEnv
conda activate bcftoolsEnv
conda install -c bioconda bcftools

# Finally, use BCFtools to remove fields but keep GT (genotypes)
bcftools annotate -x FILTER,INFO,^FORMAT/GT trio.trim.vep.vcf.gz -Oz -o trio.trim.vcf.gz
```


### Quick primer on VCF files and genotypes

Lets have a look at our VCF file.
Firstly, lets review the Variant Call Format (VCF) file.
This is the standard file for listing variants that are 'called' via a variant calling tools such as `BCFtools`, `GATK` or `freebayes`.
The basic premise for calling variants is to identify both alleles in a diploid genome.
To do this, sequenced DNA fragments ("reads") are aligned to the reference sequence and the base at each position is determined, counted and then run through a number of tests to determine whether the site is different from the reference sequence.
If the site is polymorphic, we count the bases at the position and determine whether it is a homozygous or a heterozygous variant.
In this file, a "0" denotes the reference base and "1" as the alternate base at that position.
So a heterozygous variant is "0/1" and a homozygous alternate variant is "1/1".
Additionally, it is possible to have multi-allelic sites, so additional alternate alleles are coded as greater than 1 (e.g. 2/2 or 0/2).

### VCF headers

If we have a look at the VCF file there is the header with the full VCF information.
The header is denoted by lines that start with two # (i.e. ^##).
The name of the fields for the rest of the file (that contain the actual results) is denoted by lines that start with only one # (i.e. ^#)

```bash
zcat trio.trim.vcf.gz | less
```

Headers are an amazing mass of information that comes from the variant calling process.
It will have metadata regarding the aligned reference genome, steps that were run to make the file and the definitions of the specific fields and tags within the file.

```bash
# View the header using BCFtools
bcftools view -h trio.trim.vcf.gz
```

---
___>>> QUESTIONS <<<___
---

1. What was the name of the reference genome file was used to make the VCF file?
2. What program was used to call variants?
3. What are the names of the 3 individuals that are sampled in the VCF file?
4. The VCF file is a subset to only include variants from specific chromosomes. What are the names of those chromosomes and how many individual variants are found in each?

---

### Adding a FILTER tag

Not all variants are created equal.
Much like the genotype and alignment quality metrics (base and mapping quality) that we learnt in previous lectures, the VCF file contains a QUALITY field that is also phred scaled.

>QUAL phred-scaled quality score for the assertion made in ALT. _i.e.,_ -10log_10 prob(call in ALT is wrong). If ALT is ”.” (no variant) then this is -10log_10 p(variant), and if ALT is not ”.” this is -10log_10p(no variant). High QUAL scores indicate high confidence calls. Although traditionally people use integer phred scores, this field is permitted to be a floating point to enable higher resolution for low confidence calls if desired. If unknown, the missing value should be specified. (Numeric)

So lets say that we want to warn the user that we have some variants that are probably poor quality.
We can add a tag (_i.e.,_ a bit of text) in the FILTER field to indicate that our variant is potentially poor quality.
This is important later on when you start interpreting the value of the variant.

```bash
# Put a the text "LowQual" to the FILTER tag when QUAL<30
bcftools filter -mx -sLowQual -e'%QUAL<30' trio.trim.vcf.gz
```

---
___>>> TASK <<<___
---

Variants that are located close to indels can also indicate poor quality calls, so:
- Use the `bcftools filter` sub-command to tag low quality variants (QUAL<30) that are within 10 base-pairs of an InDel, and count those variants

---

### Adding a variant ID

As you can probably imagine, each variant within the current reference genome has been extensively studied through the continual sampling of patients and individuals from around the world.
Due to this, each variant that is found within an individual sequenced over the last ~10 years has been given an rsId in the [NCBI dbSNP database](https://www.ncbi.nlm.nih.gov/snp/).
These rsIds are helpful because it gives you an extensive list of information about each particular variant.

This is actually not _technically_ correct anymore, as we have sampled so many genomes recently that the rsId numbers couldn't keep up!

There are now databases such as the [genome aggregation database (gnomAD)](https://gnomad.broadinstitute.org/) that samples over 100,000 individuals.

---
___>>> TASK <<<___
---

- Take a little bit of time to explore gnomAD, as it is an important tool in determining a baseline allele-frequency of every variant.
It is also helpful in finding potential loss-of-function (pLoF) variants.
  - Search gnomAD for the gene SATB1
  - Count the number of missense variants and list the pLoF variants in SATB1
- Look up the variant `14-82565377-G-C` (rs75115269) and identify its allele frequency in Ashkenazi Jews
- Can you find its equivalent genomics coordinates and allele frequency for the newer hg38 reference genome?

---

OK, back to annotating Ids.
If we have a database of known Id, we can easily compare to the VCF file and add text to the ID field in the VCF (3rd field after CHROM and POS).
For this we can use any type of tab-delimited file, but for this week I have provided the dbSNP reference VCF which is perfect for this task.

**NOTE:** In order to subset or retrieve data from a tab-delimited file (or any other delimited file for that matter), it is helpful to use an index.
File indexes are like phone-books (if you can remember a time when people used phone books!) which are sorted alpha-numerically to make it easy to find a name/phone number.
A file index, commonly created by the program [`tabix`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/), can be created on most standard bioinformatics files (BED, VCF etc).
Before running this annotation task, we need to also index both the database (_i.e.,_ our hg19.dbSNP.vcf.gz file) and the query (our VCF file that we want to add IDs to).
It is also good practice to create an index every time you make a new VCF file, as most tasks require it to work.
A number of variant toolkit's (`gatk`, `picard`, `sambamb`, etc...) will often create an index automatically for you.
You can either use the `tabix` program using the VCF prefix (`tabix -p vcf`) or use the `bcftools` sub-command `bcftools index`

```bash
# Index our files
bcftools index -t hg19.dbSNP.vcf.gz
bcftools index -t trio.trim.vcf.gz
```

**NOTE:** If you get an _[E::hts_idx_push] Unsorted positions on sequence #1_ error on the index command, you can quickly sort the file again using the `bcftools sort` command and reindexing.

List the files in the directory and see what is produced.

Now we can add rsIDs using the `bcftools annotation` sub-command and output a new files with our IDs attached.

```bash
# Add Ids
bcftools annotate -c CHROM,FROM,ID,REF,ALT \
    -a hg19.dbSNP.vcf.gz \
    -Oz -o trio.trim.dbSNP.vcf.gz trio.trim.vcf.gz
```

---
___>>> QUESTIONS <<<___
---

1. Is the variant we explored in gnomAD (rs75115269) present in this VCF?
2. Can you find the variant rs191680234, get its genomic coordinates, and identify for which individual it is variable?

---

## Full variant annotation

At the moment we are sequentially adding layers of information in order to give context to each of the variants that we identified in the variant calling process.
However, given the databases and information required to fully annotate the VCF, we won't run the full annotation process in this tutorial.
Instead, we downloaded a fully annotated version of our VCF at the start of the tutorial which contains annotations of each variant as executed through Variant Effect Predictor or VEP.
VEP and [snpEff](http://snpeff.sourceforge.net/) are probably the two most popular variant annotation programs and use very similar effect standards.

As the name suggests, VEP uses large variant annotation databases to assign an effect and/or consequence to the particular change, which take the form of Sequence Ontology (SO) terms.

![VEP Variant Types](https://m.ensembl.org/info/genome/variation/prediction/consequences.jpg)

As you can see from the table contained in the above link, each variant type is assigned a [level of variant _impact_ that is associated to the change](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html).
Some variants have a **HIGH** impact, including splice site changes, stop codon introductions and frameshift variants.
These we can interpret as potentially having a loss of function and therefore could influence phenotypes in a negative way.
**MODERATE** impact variants include missense variants and inframe deletions or insertions, that can have a significant change but can also lead to harmless changes.
**LOW** impact variants such as synonymous variants are likely to have a very low impact on overall the function of the protein.
Additional to protein-coding changes, non-coding or regulatory variant sequence ontologies are coded as **MODIFIERS**, although some non-coding regulatory sequence changes can produce a more significant impact.

So how are these included in the actual VCF file? Let's look:

```bash
# View the VCF that contains full annotation (without the header _i.e.,_ -H param)
bcftools view -H trio.trim.vep.vcf.gz | head

2	41647	.	A	G	4495.41	PASS	CSQ=intron_variant&non_coding_transcript_variant|||ENSG00000184731|FAM110C|ENST00000460464|||||processed_transcript|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000184731|FAM110C|ENST00000461026|||||processed_transcript|||||||||,intron_variant|||ENSG00000184731|FAM110C|ENST00000327669||||-/321|protein_coding|YES|CCDS42645.1|||||||	GT:AD:DP:GQ:PL	0/0:56,0:56:99:0,169,2183	0/1:33,35:68:99:1139,0,1044	0/1:119,117:237:99:3356,0,3283
2	45895	.	A	G	463.75	PASS	CSQ=missense_variant|aTc/aCc|I/T|ENSG00000184731|FAM110C|ENST00000327669|1/2|benign(0)|tolerated(0.62)|164/321|protein_coding|YES|CCDS42645.1|||||||,upstream_gene_variant|||ENSG00000184731|FAM110C|ENST00000460464|||||processed_transcript|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000184731|FAM110C|ENST00000461026|||||processed_transcript|||||||||	GT:AD:DP:GQ:PL	1/1:0,6:6:18.05:207,18,0	1/1:0,9:9:24.07:292,24,0	./.:.:.:.:.
2	224970	.	C	T	4241.64	PASS	CSQ=intron_variant|||ENSG00000035115|SH3YL1|ENST00000415006||||-/246|protein_coding||CCDS62842.1|||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000403657||||-/227|protein_coding||CCDS62841.1|||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000403658||||-/227|protein_coding||CCDS62841.1|||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000405430||||-/342|protein_coding|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000473104|||||processed_transcript|||||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000451005||||-/255|protein_coding|||||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000356150||||-/342|protein_coding|YES|CCDS42646.2|||||||,intron_variant&NMD_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000479739||||-/155|nonsense_mediated_decay|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000463865|||||processed_transcript|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000472012|||||processed_transcript|||||||||,downstream_gene_variant|||ENSG00000035115|SH3YL1|ENST00000431160||||-/230|protein_coding|||||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000403712||||-/323|protein_coding||CCDS54332.1|||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000468321|||||processed_transcript|||||||||GT:AD:DP:GQ:PL	0/1:40,26:66:99:789,0,1374	0/1:47,41:88:99:1247,0,1555	0/1:93,80:175:99:2205,0,2918
...
```

My eyes (glaven!).
So.....much......text......
As you can see, there is a mass of information in the INFO field, all of which starts with a consequence tag (CSQ=).
This field has a lot of information separated by pipes (|) and it is also possible to get multiple annotations per variant.
If you look at the header you can get the header information for each of these fields that are separated by |.

```bash
zgrep "^##INFO=<ID=CSQ" trio.trim.vep.vcf.gz

##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL|CCDS|RadialSVM_score|RadialSVM_pred|LR_score|LR_pred|CADD_raw|CADD_phred|Reliability_index">
```

---
___>>> QUESTIONS <<<___
---

1. How many variants have a reported "missense_variant"?
2. How many variants that are greater than quality 30 have a reported "missense_variant"?
3. How many variants are annotated to have gained a stop codon?
4. Can you find the variant we looked at before (rs191680234) and tell what type of variant it is and if it passes QC

---

### Variant scores

Of course, the impact of the variant may not always tell you the definitive pathogenicity of that variant.
LOW impact variants may also have a high impact within some systems, so a number of additional metrics have been established to add additional interpretation power to variant annotation.
These include:

- GERP
- CADD
- SIFT
- PolyPhen
- REVEL

[All the additional fields are here](https://m.ensembl.org/info/genome/variation/prediction/protein_function.html)

### CADD

The Combined Annotation Dependent Depletion (CADD) tool scores the predicted deleteriousness of single nucleotide variants and insertion/deletions variants in the human genome by integrating multiple annotations including conservation and functional information into one metric.
Phred-style CADD raw scores are displayed and variants with higher scores are more likely to be deleterious.
CADD provides a ranking rather than a prediction or default cut-off, with higher scores more likely to be deleterious.
For example, scores above 30 are 'likely deleterious' and scores below as 'likely benign'.
Variants with scores over 30 are predicted to be the 0.1% most deleterious possible substitutions in the human genome.

