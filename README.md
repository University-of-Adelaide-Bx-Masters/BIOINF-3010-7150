* TOC
{:toc}

## BIOINF3010/7150: Genomics Applications
{:.no_toc}

Semester 1 2020

### Practical Timetable

| **Week** | **Monday** | **Practical**                                     |
|----------|------------|---------------------------------------------------|
| **1**    | 2/3        | [Introduction to Bash] (Dan)                      |
| **2**    | 9/3        | [Read Quality Control] (Nathan) [SAMTools and alignments] (Jimmy) |
| **3**    | 16/3       | [SARS-CoV-2 Resequencing] (Nathan) [SARS-CoV-2 Short Read Assembly] (Nathan)         |
| **4**    | 23/3       | [Short and long read alignment] (Nathan) [E. coli K-12 Hybrid Genome Assembly] (Nathan)         |
| **5**    | 30/3       | [Bacterial genome assembly] (Lloyd)               |
| **6**    | 6/4        | [HiC analysis] (Lloyd/Ning)                       |
| **-**    |            |                                                   |
| **7**    | 27/4       | [Genome graphs1] (Yassine) and [Genome graphs2] (Yassine)                         |
| **8**    | 4/5        | [BLAST analysis and databases] (Dave)             |
| **9**    | 11/5       | [Clinical genomics] (Jimmy)                       |
| **10**   | 18/5       | [Agricultural genomics] (?)                       |
| **11**   | 25/5       | [Population genetics] (Bastien)                   |
| **12**   | 1/6        | [Metagenomics 16S profiling] (Raphael)            |

[Introduction to Bash]: Practicals/Bash_Practicals/1_IntroBash.md
[Read Quality Control]: Practicals/Read_QC/read-qc.md
[SAMTools and alignments]: Practicals/Alignments_Practicals/alignment-cram.md
[SARS-CoV-2 Resequencing]: Practicals/resequencing/resequencing.md
[SARS-CoV-2 Short Read Assembly]: Practicals/short_read_assembly/short-read-assembly.md
[Short and long read alignment]: Practicals/short_long_alignment/short_long_alignment.md
[E. coli K-12 Hybrid Genome Assembly]: Practicals/hybrid_genome_assembly/index.md
[Bacterial genome assembly]: Practicals/
[HiC analysis]: Practicals/HiC/hicprac.md
[Genome graphs1]: Practicals/Graph_Genomes/prac_part1.md
[Genome graphs2]: Practicals/Graph_Genomes/prac_part2.md
[BLAST analysis and databases]: Practicals/
[Clinical genomics]: Practicals/
[Agricultural genomics]: Practicals/
[Population genetics]: Practicals/
[Metagenomics 16S profiling]: Practicals/

## Assessment

### Assessment Tasks

| **Assessment**                                            | **Subject**         |
|-----------------------------------------------------------|---------------------|
| [Assessment 0](Assignments/Assignment_0/bash_questions.md)| Bash                |
| [Assessment 1](Assignments/Assignment_1/genome_sequencing.md)                                          | Genome sequencing     |
| [Assessment 2]()                                          | Experimental design |
| [Assessment 3]()                                          | Publishing a genome |
| [Assessment 4]()                                          | Ancient DNA         |
| [Assessment 5]()                                          | Metagenomics        |
| [Project](#project) (PG only)                                     | Complete Dataset    |

### Project

In this course, the following next-generation sequencing (NGS) datasets/protocols will be examined in detail:

- Whole genome sequencing/Resequencing
- SNV variant/structural variation analysis
- Enrichment/Capture sequencing (Methyl-capture, ChIPseq, RIPseq)
- Metagenomics/Microbial profiling

Each of these NGS approaches uses similar programs and analysis approaches, such as quality control (quality and sequencing adapter trimming), genome alignment, and downstream visualisation and statistical methods.
For the project, you will take a published (or otherwise obtained) dataset and complete all the analysis tasks (from raw data to final results) during the course.
You have the freedom to choose any dataset you would like. You will prepare a final report that will be due at the end of the semester.
The report should be prepared using RStudio as an Rmd document including all code needed to perform the analysis, and will include the standard components of a scientific report:

- Introduction (background on the study and identification of the research hypothesis)
- Methods (analysis steps and programs used)
- Results (what you found) and; 
- Discussion (how the results relate to the research hypothesis and the published literature).

The Rmd document and a compiled knitted html will form the submission; marks will be awarded to the code and Rmd that you use.

| Section | Mark |
|---------|-----:|
| Abstract | 5% |
| Introduction and hypothesis |	10% |
| Methods | 20% |
| Results and Discussion | 30% |
| References | 5% |
| Analysis scripts | 30% |

#### Project Data

For the project I was able to download a number of publicly available datasets from the Encylopedia of DNA elements (ENCODE) project, which is a large multi-national study that wrapped up a while ago.
The purpose of the study was to identify any "functional" region of the genome that may not be gene-coding, so the project sequenced a lot of RNA sequencing, ChIP-seq (Transcription Factor-binding), DNA methylation sequencing and arrays etc.

##### Sample Information

GM12878 is a human lymphoblastoid cell-line, a component of the human Lymphoblastic Leukaemias, taken from a large family from Utah (Central European Ancestry) in 1985.
These cell-lines are widely used in genomics as reference sets for large projects and are easy to obtain and use in a research setting.

##### RNA-seq

In the data directory you will find a range of RNA-seq and ChIPseq data from the human cell-line GM12878.
ENCODE datasets were produced back in 2012 by a number of labs in the US.
They include RNA-seq from four different RNA fractions:

- Long PolyA+ enriched RNA from Whole-cells
- Long RNA from Whole-cells without PolyA enrichment
- Short total RNA
- Long total RNA

Short vs Long refers to the size selection of the RNA before making the library.
Short is generally less than 100bp and large is >100bp.

All of the library protocols are available already so you can have a look at the specifics (https://public-docs.crg.eu/rguigo/Data/jlagarde/encode_RNA_dashboard//hg19/).

For differential expression, there is 6 samples from the paper "Cis-Regulatory Circuits Regulating NEK6 Kinase Overexpression in Transformed B Cells Are Super-Enhancer-Independent" by Huang et al. 2017 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5393904/).
These GM12878 cells are the same as above, with one group of 3 clones from normal cells, and the other group of 3 clones with a deleted region.

##### ChIP-seq

If you would like to do something slighly different, I have also included two ChIP-seq datasets that enrich for CTCF transcription factor binding sites (https://en.wikipedia.org/wiki/CTCF).
CTCF is an important TF for structural organisation of the chromosome and is used a lot on 3D chromosome structure analyses (3C/4C/5C/HiC-seq).

Each replicate is also sampled on GM12878.

---

All the data is available from the following link:
https://universityofadelaide.box.com/v/mscProjectData

**Note 1:** The data from this directory is approximately 100GB, meaning that you cannot download the data in one go. I would suggest choosing specific libraries you would like to work on and download those separately onto your VM so you don't fill up the VM's allocated space.

**Note 2:** Some of the data is from 2012-2014, so some of the sequencing technology is quite old!

Good luck!


### Assessment Checklist

Have you:

- [ ] Answered all the questions?
- [ ] Followed naming conventions for Assessments?
- [ ] Checked that you have not breached the [Academic Honesty Policy](http://www.adelaide.edu.au/policies/230/).
- [ ] Identified the work as yours?
	- Emails should have the course and assessment task names.
	- Documents should be named with your name, the course name and the assessment task.
	- Printed documents should have you name and the course and assessment task in the text/footer/header.
- [ ] Used appropriate electronic communication with assessors?
	- Emails should have a meaningful subject.
- [ ] Handed in the assignment before the due time (see MyUni)?

## Useful Links

[How To Ask Questions The Smart Way](http://www.catb.org/esr/faqs/smart-questions.html)

[How to write a good bug report](https://musescore.org/en/developers-handbook/how-write-good-bug-report-step-step-instructions)
