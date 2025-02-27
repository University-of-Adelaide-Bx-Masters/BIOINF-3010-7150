* TOC 
{:toc}

## BIOINF3010/7150: Genomics Applications
{:.no_toc}

Semester 1 2025 - *Provisional Timetable*

### Coordinator:

- Anna Sheppard

#### Instructors:

- Anna Sheppard
- Jessica Webb
- Kat Pillman
- Melanie Smith
- Lloyd Low
- Zhipeng Qu
- Chelsea Matthews
- Jacqueline Rehn
- Julien Soubrier
- Bastien Llamas

#### Course contact Anna Sheppard
- email: anna.sheppard@adelaide.edu.au
- Phone: 08 8313 0232
- Office: Molecular Life Sciences Rm 4.44

#### Contact hours
Lectures: Monday 9am, Benham G10
- 3 Mar --- 2 Jun

Practical 01: Wednesday 9am-11am, Ingkarni Wardlii 218  
- 5 Mar --- 11 Jun

Practical 02: Thursday 9am-11am, Ingkarni Wardlii 218
- 6 Mar --- 12 Jun

### Lecture and Practical Timetable

| **Week** | **Monday** | **Lecture**                                       | **Practical**                                     |
|----------|------------|---------------------------------------------------|---------------------------------------------------|
| **1**    | 3/3        | **No live lecture!**<br>Introduction to sequencing - Molecular basis of polymer extension sequencing (Sanger) High throughput sequencing: Illumina Sequence quality ***(recording)*** | [Introduction to Bash 1] (Anna)   [Introduction to Bash 2] (Anna)                  |
| **2**    | 10/3       | **No lecture - Public Holiday**                   |[Read Quality Control] (Kat) [Read Quality Control] (Kat) |
| **3**    | 17/3       | Resequencing (Exome/WGS)     ***(recording)***  		| [SAMTools and alignments] (Melanie)   [SARS-CoV-2 Resequencing] (Melanie)      |
| **4**    | 24/3       | Annotation - Gene finding, Repeat identification/classification/masking, comparative genomics ***(recording)***  |[Intro to BLAST] (Jess)     [BLAST practical] (Jess)          |
| **5**    | 31/3       | Short read assembly; approaches and issues ***(recording)*** | [Short and long read alignment]  (Anna)   [SARS-CoV-2 Short Read Assembly] (TBD)     |
| **6**    | 7/4        | Single molecule sequencing PacBio/Nanopore Uses/Characteristics/Error profiles ***(recording)***  |  PG Major Project information (Chelsea) [E. coli K-12 Hybrid Genome Assembly] (Lloyd)     |
| **-**    |            |  |    **MID-SEMESTER BREAK**                                             |
| **7**    | 28/4       | _De novo_ assembly genome size estimate (k-mers), coverage (Zhipeng)  |    [Genome Assembly I] (Zhipeng) [Genome Assembly II] (Zhipeng)        |
| **8**    | 5/5        | Variant calling and high-throughput genotyping (Julien)  |  [Clinical genomics1] (Julien)  [Clinical genomics2] (Julien)         |
| **9**    | 12/5       | Genome Graphs (Chelsea)  |  [Genome graphs1] (Chelsea)  [Genome graphs2] (Chelsea)     |
| **10**   | 19/5       | Structural variation and cancer genomics (Jacqueline) |   [Structural variation] (Jacqueline) [Structural variation] (Jacqueline)          |
| **11**   | 26/5       | Population genomics (Bastien) | [Population genomics1] (Bastien)  [Population genomics2] (Bastien)              |
| **12**   | 2/6        | Wrap up lecture (TBD) |  Open Prac session for PG major projects only  Wednesday (Chelsea)  Open Prac session for PG major projects only  - Thursday (Chelsea)           |
| **13**   | 9/6        | **No lecture - Public Holiday**         			|   TBD        |

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
[Structural variation]: Practicals/Structural_variation/SV_practical.md
[Clinical genomics1]: Practicals/variants_clinical/variant_annotation.md
[Clinical genomics2]: Practicals/variants_clinical/variant_filtering.md
[Agricultural genomics]: Practicals/
[Population genomics1]: Practicals/ancient_DNA_pop_genomics/prac_part1.md
[Population genomics2]: Practicals/ancient_DNA_pop_genomics/prac_part2.md


## Assessment

### Assessment Tasks (This section for information only, links are disabled. Links to assessment items are availabe on MyUni)

| **Assignment**                                            | **Subject**         |
|-----------------------------------------------------------|---------------------|
| [Assignment 0]()                                          | Bash                |
| [Assignment 1]()                                          | Genome sequencing   |
| [Assignment 2]()                                          | Genome assembly     |
| [Assignment 3]()                                          | Genome graphs and SV|
| [Assignment 4]()                                          | Clinical and Population Genomics |
| [Project](./Assignments/Major_Project/major_project.md) (PG only)                                     | Complete Dataset    |

- Each assignment is worth 20% of the final mark for Undergraduates.  
- Each assignment is worth 12% of the final mark for Masters students.  
- Major Project (Masters students only) is worth 40% of the final mark.  

### Major Project (Post Grad only)

The following next-generation sequencing (NGS) datasets/protocols may be available for the major project:

- Whole genome sequencing/Resequencing
- Whole genome assembly
- Sequence alignment
- Ancient DNA/Population Genomics
- Genome graphs
- SNV variant/structural variation analysis

Each of these NGS approaches uses similar programs and analysis approaches, such as quality control (quality and sequencing adapter trimming), genome alignment, and downstream visualisation and statistical methods.
For the project, you will take a published (or otherwise obtained) dataset and complete all the analysis tasks (from raw data to final results) during the course.
You have the freedom to choose any dataset you would like from the scientific literature, but you need to select and get approval for your dataset by the end of week 6. You will prepare a final report that will be due at the end of the semester.
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

Look [here](./Assignments/Major_Project/major_project.md) if you want to use provided datasets. 

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
