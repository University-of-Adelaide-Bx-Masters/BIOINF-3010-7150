# Major Project (40% [10% project selection] + [30% project assessment])

For the major project, you will take a published dataset and complete all the analysis tasks (from raw data to final results) and write up a report. This report must be structured like a small journal article, with abstract (summarising the project), introduction (background on the study and identification of the research hypothesis), methods (analysis steps and programs used), results (what you found) and discusson (how the results relate to the research question) sections. Marks will also be awarded to the bash/R or RMarkdown scripts that you use.

|Section                    |Mark |
|:--------------------------|:----|
|Abstract                   |5%   |
|Introduction + hypothesis  |10%  |
|Methods                    |20%  |
|Results & Discussion       |30%  |
|References                 |5%   |
|Analysis scripts           |30%  |

**You have the freedom to choose any dataset from any research article you would like**, however you must get approval beforehand.

If you cannot find a suitable dataset, we have provided the following datasets:

# De-novo genome assembly and genome graphs

## Introduction

In this major project, you will *de novo* assemble the genome of a strain of Brewer's Yeast using a 30x Third generation sequencing dataset (Nanopore). The sequencing dataset for genome assembly is from the *Saccharomyces cerevisiae* reference assembly panel [(ScRAP) project](https://www.nature.com/articles/s41588-023-01459-y). After you assemble the genome of your designated strain, you will need to assess the assembly quality by using the approaches that you have learned in this course (e.g. align to reference genome, QUAST, BUSCO). Then you will need to call variants using two given Next Generation Sequencing datasets (Illumina) from [two other strains](https://www.ncbi.nlm.nih.gov/bioproject/379572). After you identify the variants from those two different strains, you will need to build genome graphs for the assembled genome of your designated strain using a small subset of detected variants (see following instructions for more details on getting a suitably small subset), and then compare the read mapping of a test Illumina dataset when using genome graphs with integrated variants and without integrated variants. The principal research question is **How good is your *de novo* assembled reference genome and whether using genome graphs with integrated variants will improve the read mapping?**  You will note that the instructions for this assignment are quite general, this is because you are expected to use the knowledge you have gained in the practicals and from previous assignments to select specific methods/tools for your analysis. This is what makes this a Major Project. 


## Datasets

| Dataset                                            | Platform            | Sequence type | Usage                   | Location                                               |
|----------------------------------------------------|---------------------|---------------|-------------------------|--------------------------------------------------------|
| ERR85624XX_30x.fastq.gz                            | Nanopore PromethION | Single end    | De novo genome assembly | ~/data/major_project_2024/01_raw_data/01_TGS/axxxxxxx/ |
| ERR1938683_1.fastq.gz, ERR1938683_2.fastq.gz       | Illumina MiSeq      | PE150         | Variant calling         | ~/data/major_project_2024/01_raw_data/02_NGS/          |
| ERR1938686_1.fastq.gz, ERR1938686_2.fastq.gz       | Illumina MiSeq      | PE150         | Variant calling         | ~/data/major_project_2024/01_raw_data/02_NGS/          |
| Illumina_R1.10x.fastq.gz, Illumina_R2.10x.fastq.gz | Illumina MiSeq      | PE150         | Genome graphs mapping   | ~/data/major_project_2024/01_raw_data/02_NGS/          |

*[Important] Each student will be provided with a uniqe Nanopore sequencing dataset from a unique strain of S. cerevisiae for de novo genome assembly. You can find your unique input dataset in the folder named with your "student a number" in `~/data/major_project_2024/01_raw_data/01_TGS/`. Please double-check that you do copy/use your **correct** unique dataset before you proceed with any further analysis.*

You will also be provided with the genome sequence of the reference S288c strain in the folder of `~/data/major_project_2024/02_databases/` for comparison with your *de novo* assembled genome.

## Recommended analysis pipeline

### Step 1, de-novo genome assembly

*[Reference Prac, `Genome Assembly I` and `Genome Assembly II`]* 

Use the knowledge that you have learned in Prac `Genome Assembly I` and `Genome Assembly II` to *de novo* assemble the reference genome of your designated strain. Assess your assembly quality  (If you want to run QUAST with reference, you can use the provided S288c genome as reference). 

### Step 2, variant calling

*[Reference Prac, `Read Quality Control`, `SAMTools and alignments`, `Short and Long Read Alignment`, and `Genome Assembly II`]* 

You will use the assembled genome of your designated strain from step 1 as reference, and call variants using the Illumina sequencing reads from two other strains (ERR1938683 for strain S288C, and ERR1938686 for strain SK1). Please refer to Prac `Read Quality Control` and `SAMTools and alignments` for how to check the quality of Illumina sequencing reads and how to align the Illumina short reads to a genome. Then you can use `bcftools` to call variants for these other two strains based on your designated strain reference genome (see Prac `Genome Assembly II` for how to call variants). Report the summary of called variants.

### Step 3, mapping using genome graphs

*[Reference Prac, `Genome graph1`, `Genome graph2`]* 

You will need to construct two types of genome graphs for the reference genome of your designated strain using `vg`. One genome graph is without variants integrated, and the other one is with variants integrated (using the variants you detected in step 2). **Important: due to the computing limitation of our teaching VM, we can only use a small subset of variants called in Step 2 to construct the genome graphs (You can use the command `bcftools filter -i 'DP>500 && QUAL>100' your_called_variants.vcf.gz -O z -o variants_for_genome_graphs.vcf.gz` to get a small subset of all detected variants in step 2).** Then map the provided test Illumina sequencing dataset (Illumina_R1.10x.fastq.gz and Illumina_R2.10x.fastq.gz, from strain S288C) using `vg` to the constructed genome graphs to compare the read mapping, and answer the question of "whether genome graphs (with variants included) can improve the read mapping" (Hint: you can use `jq` and `.identity` that you have learned to check the read mapping identity for genome graphs with variants and without variants). 

## Some comments

1. Set up your project folder structure in an organised way.

2. You are welcome to fine-tune different parameters in your analysis, for example, when *de novo* assembling your genome using `flye`, or doing short read alignment using `BWA`. However, it is fine if you just stick with the default settings.

3. You are welcome to include additional analysis in your project, for example, gene annotation using `BLAST` from the assembled genome to identify genes/proteins of interest (see Prac `BLAST`). However this is optional.

4. You may want to include figures/tables in your report, however, please don't try to include every figure/table from your analysis in your report. The general rule for figures is no more than 6-8 figures (you can have multiple sub-panels in one figure) in a scientific article.

5. There is no hard word limit for your report. However, please remember sometimes **Less is more**.
