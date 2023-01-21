* TOC 
{:toc}

## BIOINF3010/7150: Genomics Applications
{:.no_toc}

Semester 1 2023 - *Provisional Timetable*

### Lecture and Practical Timetable

| **Week** | **Monday** | **Lecture**                                       | **Practical**                                     |
|----------|------------|---------------------------------------------------|---------------------------------------------------|
| **1**    | 27/2       | Introduction to sequencing Molecular basis of polymer extension sequencing (Sanger)High throughput sequencing: Illumina Sequence quality (Dave) | [Introduction to Bash 1] (Dave)   [Introduction to Bash 2] (Dave)                  |
| **2**    | 6/3        | Resequencing (Exome/WGS)     (Dave)                    |[Read Quality Control] (Dave)[Read Quality Control contd] |
| **3**    | 13/3       | **No lecture - Public Holiday**  |[SARS-CoV-2 Resequencing] (Dave)  [Intro to BLAST] (Dave)       |
| **4**    | 20/3       | Annotation - Gene finding, Repeat identification/classification/masking, comparative genomics (Dave)  |[BLAST practical] (Dave) [SAMTools and alignments] (Dave)          |
| **5**    | 27/3       | Short read assembly; approaches and issues (Dave)  |  [SARS-CoV-2 Short Read Assembly] (Dave) [Short and long read alignment]  (Anna)         |
| **6**    | 3/4        | Single molecule sequencing PacBio/Nanopore Uses/Characteristics/Error profiles (Anna)  | [E. coli K-12 Hybrid Genome Assembly] (Anna)  **Friday no practical -public holiday**                         |
| **-**    |            |  |    MID-SEMESTER BREAK                                             |
| **7**    | 24/4       | _De novo_ assembly genome size estimate (k-mers), coverage. (Zhipeng)  | **Tuesday no practical - public holiday**   [Genome Assembly I] (Zhipeng)         |
| **8**    | 1/5        | Genome Graphs (Chelsea)  | [Genome Assembly II] (Zhipeng) [Genome graphs1] (Chelsea)           |
| **9**    | 8/5        | Structural variation and cancer genomics (Paul)  |  [Genome graphs2] (Chelsea) [Structural variation] (Paul)      |
| **10**   | 15/5       | Variant calling and high-throughput genotyping (Julien) |  [Clinical genomics1] (Julien)  [Clinical genomics2] (Julien)         |
| **11**   | 22/5       | Population genomics (Yassine) | [Population genomics1] (Yassine)  [Population genomics2] (Yassine)              |
| **12**   | 29/5       | Wrap up lecture (Dave) |  Open Practical Session Tuesday (Dave)  [Open  Prac session - Friday] (Dave)           |
| **13**   | 5/6        |   **Monday public holiday**         |   TBD        |

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

### Assessment Tasks

| **Assignment**                                            | **Subject**         |
|-----------------------------------------------------------|---------------------|
| [Assignment 0]()                                          | Bash                |
| [Assignment 1]()                                          | Genome sequencing   |
| [Assignment 2]()                                          | Genome assembly     |
| [Assignment 3]()                                          | Genome graphs and SV|
| [Assignment 4]()                                          | Clinical and Population Genomics |
| [Project]() (PG only)                                     | Complete Dataset    |

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

Look [here](./Assignments/Major_Project/major_project.md) for an example of major project data. 

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
