# Assignment 1: Genome Sequencing

1. Describe the 4 main components of a FASTQ read/record.
2. Illumina short reads suffer from a deterioration in quality towards the 3' end. Describe the process which causes this.
3. Illumina short reads may contain portions of adapter sequences at their 3' end. Describe how and why some reads may contain parts of an adapter while others may not.
4. Compare and contrast SAM and BAM files
5. Index files are regularly encountered in bioinformatics. For example, `.bai` (and `.csi`) is an index file for BAM files and `.fai` are index files of FASTA files. Describe, in general terms, what index files are and what they facilitate.
6. Oxford Nanopore (ONT) and Pacific Biosciences (PacBio) are both single molecule, long read sequencing technologies. Describe their error profile and how this differs from that of Illumina short reads.
7. Following the generation of WGS sequencing reads you may choose to a) align reads to a reference genome or b) perform a de novo genome assembly. Compare and contrast these two different approaches and describe why you might choose one over the other.
8. Why are paired-end reads considered superior to single-end reads?
9. The contig N50 is often reported for a genome assembly. Describe how the N50 is calculated and why it is not a measure of assembly quality/accuracy.
10. Describe how PacBio reads can achieve higher quality/accuracy compared to ONT reads. What is the trade-off for getting higher quality reads?
11. ONT sequencing relies on detecting step-changes in current as DNA passes through a nanopore embeded within a membrane. Describe a situation when it may not be possible to detect a step-change in current as a DNA strand passes through the nanopore.
12. ONT sequencing suffers from high error rates as well as systemmatic errors. Describe approaches that are used to try and reduce these errors.
13. In the Hybrid Genome Assembly practical, you began to explore the effects of using different amounts of input Illumina and PacBio data. Continue this investigation to generate 3 genome assemblies with different amounts of input data and report the following for each:
 * Amount of Illumina data used
 * Amount of PacBio data used
 * How many contigs were generated?
 * How many scaffolds were generated?
 * Compare each assembly to the Reference genome using MUMmer and plot the `.delta` file using R or Assemblytics. Include these in your submission.