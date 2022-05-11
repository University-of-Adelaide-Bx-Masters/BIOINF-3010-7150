* TOC 
{:toc}

## BIOINF3010/7150: Genomics Applications
{:.no_toc}

Semester 1 2022 - *Provisional Timetable*

### Lecture and Practical Timetable

| **Week** | **Monday** | **Lecture**                                       | **Practical**                                     |
|----------|------------|---------------------------------------------------|---------------------------------------------------|
| **1**    | 28/2       | Introduction to sequencing Molecular basis of polymer extension sequencing (Sanger)High throughput sequencing: Illumina Sequence quality (Dave) | [Introduction to Bash 1] (Dave)   [Introduction to Bash 2] (Dave)                  |
| **2**    | 7/3        | Resequencing (Exome/WGS)     (Dave)                    |[Read Quality Control] (Dave)[Read Quality Control contd] |
| **3**    | 14/3       | **No lecture - Public Holiday**  |[SARS-CoV-2 Resequencing] (Dave)  [TBA]       |
| **4**    | 21/3       | Short read assembly; approaches and issues (Dave) | [SAMTools and alignments] (Dave)  [SARS-CoV-2 Short Read Assembly] (Dave)          |
| **5**    | 28/3       | Single molecule sequencing PacBio/Nanopore Uses/Characteristics/Error profiles (Anna) |  [Short and long read alignment] (Anna) [E. coli K-12 Hybrid Genome Assembly] (Anna)         |
| **6**    | 4/4        | _De novo_ assembly genome size estimate (k-mers), coverage. (Chelsea) | [Genome Assembly I] (Chelsea) [Genome Assembly II] (Chelsea)                        |
| **-**    |            |  |    MID-SEMESTER BREAK                                             |
| **7**    | 25/4       | **No lecture - Public Holiday** _Lecture in Tuesday prac session ->_ |  **Genome Graphs Lecture - Tuesday  (Yassine)**   [Genome graphs1] (Yassine)                    |
| **8**    | 2/5        | Annotation - Gene finding, Repeat identification/classification/masking, comparative genomics (Dave) | [Genome graphs2] (Yassine) [Intro to BLAST] (Dave)            |
| **9**    | 9/5        | Structural variation and cancer genomics (Paul)  |   [BLAST practical] (Dave) [Structural variation] (Paul)                       |
| **10**   | 16/5       | Variant calling and high-throughput genotyping (Julien) |  [Clinical genomics1] (Julien)  [Clinical genomics2] (Julien)                   |
| **11**   | 23/5       | Population genomics (Bastien) | [Population genomics1] (Bastien)  [Population genomics2] (Bastien)              |
| **12**   | 30/5       | Wrap up lecture (Dave) |  Open Practical Session Tuesday (Dave)  [Open  Prac session - Friday] (Dave)           |
| **13**   | 6/6        |   TBD         |   TBD        |

[Introduction to Bash 1]: Practicals/Bash_Practicals/1_IntroBash.md
[Introduction to Bash 2]: Practicals/Bash_Practicals/2_BashScripting.md
[Read Quality Control]: Practicals/Read_QC/read-qc.md
[SAMTools and alignments]: Practicals/Alignments_Practicals/alignment-cram.md
[SARS-CoV-2 Resequencing]: Practicals/resequencing/resequencing.md
[SARS-CoV-2 Short Read Assembly]: Practicals/short_read_assembly/short-read-assembly.md
[Short and long read alignment]: Practicals/short_long_alignment/short_long_alignment.md
[E. coli K-12 Hybrid Genome Assembly]: Practicals/hybrid_genome_assembly/index.md
[Bacterial genome assembly]: Practicals/
[Genome Assembly I]: Practicals/Genome_assembly/genome_assembly_prac_1.md
[Genome Assembly II]: Practicals/Genome_assembly/genome_assembly_prac_2.md
[Genome graphs1]: Practicals/Graph_Genomes/prac_part1.md
[Genome graphs2]: Practicals/Graph_Genomes/prac_part2.md
[Intro to BLAST]: Practicals/BLAST_practical/BLAST_intro.md
[BLAST practical]: Practicals/BLAST_practical/BLAST_practical_v2.md
[Structural variation]:Practicals/Structural variation/SV_practical.md
[Clinical genomics1]: Practicals/variants_clinical/variant_annotation.md
[Clinical genomics2]: Practicals/variants_clinical/variant_filtering.md
[Agricultural genomics]: Practicals/
[Population genomics1]: Practicals/ancient_DNA_pop_genomics/prac_part1.md
[Population genomics2]: Practicals/ancient_DNA_pop_genomics/prac_part2.md


## Assessment

### Assessment Tasks

| **Assignment**                                            | **Subject**         |
|-----------------------------------------------------------|---------------------|
| [Assignment 0]()                                          | Bash                |
| [Assignment 1]()                                          | Genome sequencing     |
| [Assignment 2]())                                          | Genome assembly|
| [Assignment 3]()                                          | Genome Annotation|
| [Assignment 4]()                                          | Population Genomics  (No Link)       |
| [Project]() (PG only)                          | Complete Dataset   |

### Major Project (Post Grad only)

In this course, the following next-generation sequencing (NGS) datasets/protocols will be examined in detail:

- Whole genome sequencing/Resequencing
- Whole genome assembly
- Sequence alignment
- Ancient DNA/Population Genomics
- Genome graphs
- SNV variant/structural variation analysis

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

Look [here](./Assignments/Major_Project/major_project.md) for major project data. 

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
