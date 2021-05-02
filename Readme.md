* TOC
{:toc}

## BIOINF3010/7150: Genomics Applications
{:.no_toc}

Semester 1 2021 - Provisional Timetable

### Lecture and Practical Timetable

| **Week** | **Monday** | **Lecture**                                       | **Practical**                                     |
|----------|------------|---------------------------------------------------|---------------------------------------------------|
| **1**    | 1/3        | Introduction to sequencing Molecular basis of polymer extension sequencing (Sanger)High throughput sequencing: Illumina Sequence quality (Dave) | [Introduction to Bash 1] (Dave)   [Introduction to Bash 2] (Dave)                  |
| **2**    | 8/3        | Resequencing (Exome/WGS)     (Dave)                    |[Read Quality Control] (Dave) [SAMTools and alignments] (Joe/Dave) |
| **3**    | 15/3       | Short read assembly; approaches and issues (Dave |[SARS-CoV-2 Resequencing] (Dave) [SARS-CoV-2 Short Read Assembly] (Dave)         |
| **4**    | 22/3       | Single molecule sequencing PacBio/Nanopore Uses/Characteristics/Error profiles (Dave)|[Short and long read alignment] (Dave) [E. coli K-12 Hybrid Genome Assembly] (Dave)         |
| **5**    | 29/3       | Friday public holiday, no lecture |  [Tuesday open practical session] (Dave) [Friday public holiday no practical]         |
| **6**    | 5/4        | _De novo_ assembly genome size estimate (k-mers), coverage. (Lloyd) | [Assembly practical pt 1] (Lloyd) [Assembly practical pt 2] (Lloyd)                        |
| **-**    |            |  |    MID-SEMESTER BREAK                                             |
| **7**    | 26/4       | Genome Graphs (Yassine) |  [HiC analysis] (Callum)   [Genome graphs1] (Yassine)                    |
| **8**    | 3/5        | Annotation - Gene finding, Repeat identification/classification/masking, comparative genomics (Dave) | [Genome graphs2] (Yassine) [Intro to BLAST] (Dave)            |
| **9**    | 10/5       | Variant calling and high-throughput genotyping (Julien) | [BLAST practical] (Dave) [Clinical genomics1] (Julien)                     |
| **10**   | 17/5       | High-throughput genotyping technologies and applications (Rick)|[Clinical genomics2] (Julien) [Agricultural genomics] (Rick)                       |
| **11**   | 24/5       | Population genomics (Bastien) | Open Practical Session Tuesday (Dave) [Population genetics1] (Bastien)              |
| **12**   | 31/5       | Wrap up lecture (Dave) |   [Population genetics2] (Bastien) [Open  Prac session - Friday] (Dave)           |
| **13**   | 7/6        |   TBD         |   TBD        |

[Introduction to Bash 1]: Practicals/Bash_Practicals/1_IntroBash.md
[Introduction to Bash 2]: Practicals/Bash_Practicals/2_BashScripting.md
[Read Quality Control]: Practicals/Read_QC/read-qc.md
[SAMTools and alignments]: Practicals/Alignments_Practicals/alignment-cram.md
[SARS-CoV-2 Resequencing]: Practicals/resequencing/resequencing.md
[SARS-CoV-2 Short Read Assembly]: Practicals/short_read_assembly/short-read-assembly.md
[Short and long read alignment]: Practicals/short_long_alignment/short_long_alignment.md
[E. coli K-12 Hybrid Genome Assembly]: Practicals/hybrid_genome_assembly/index.md
[Bacterial genome assembly]: Practicals/
[HiC analysis]: Practicals/HiC/hic-pro_prac.md
[Genome graphs1]: Practicals/Graph_Genomes/prac_part1.md
[Genome graphs2]: Practicals/Graph_Genomes/prac_part2.md
[Intro to BLAST]: Practicals/BLAST_practical/BLAST_intro.md
[BLAST practical]: Practicals/BLAST_practical/BLAST_practical_v2.md
[Clinical genomics1]: Practicals/variants_clinical/variant_annotation.md
[Clinical genomics2]: Practicals/variants_clinical/variant_filtering.md
[Agricultural genomics]: Practicals/
[Population genetics1]: Practicals/ancient_DNA_pop_genomics/prac_part1.md
[Population genetics2]: Practicals/ancient_DNA_pop_genomics/prac_part2.md


## Assessment

### Assessment Tasks

| **Assessment**                                            | **Subject**         |
|-----------------------------------------------------------|---------------------|
| [Assessment 0]()                                          | Bash                |
| [Assessment 1]()                                          | Genome sequencing     |
| [Assessment 2]()                                          | CANU genome assembly|
| [Assessment 3]()                                          | Genome Annotation|
| [Assessment 4]()                                          | Ancient DNA  (No Link)       |
| [Project](#project) (PG only)                             | Complete Dataset    |

### Major Project (Post Grad only)

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

#### Major Project Data

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

All of the library protocols are available already so you can have a look at the [specifics](https://public-docs.crg.eu/rguigo/Data/jlagarde/encode_RNA_dashboard//hg19/).

For differential expression, there is 6 samples from the paper "Cis-Regulatory Circuits Regulating NEK6 Kinase Overexpression in Transformed B Cells Are Super-Enhancer-Independent" by [Huang et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5393904/).
These GM12878 cells are the same as above, with one group of 3 clones from normal cells, and the other group of 3 clones with a deleted region.

##### ChIP-seq

If you would like to do something slighly different, I have also included two ChIP-seq datasets that enrich for [CTCF transcription factor binding sites](https://en.wikipedia.org/wiki/CTCF).
CTCF is an important TF for structural organisation of the chromosome and is used a lot on 3D chromosome structure analyses (3C/4C/5C/HiC-seq).

Each replicate is also sampled on GM12878.

---

All the data is available [here]https://universityofadelaide.box.com/v/mscProjectData

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

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
