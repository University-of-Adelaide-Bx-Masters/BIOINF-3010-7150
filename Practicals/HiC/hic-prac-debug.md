# HiC-Pro Run through.

Please email error messages to: a1767591@adelaide.edu.au

***Please note, the `$` means that line is input. If copying and pasting commands, please omit the `$`.***

## 1. Create Anaconda Environment

``` bash
$ conda create -n hic python=3 biopython numpy
$ conda activate hic
$ pip install fithic
```

## 2. Download necessary files
``` bash
$ cd

$ mkdir utils && cd utils

$ wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.19.02.jar

$ wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/hicpro2fithic.py

$ wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/hicpro2juicebox.sh

$ wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/digest_genome.py

$ cd

$ wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata.tar.gz

$ wget http://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/hicpro_2.11.4_ubuntu.img

$ cd
```

## 3. Ensure we have the required files.
- Generate genome fragment file. Double check hg19 is uncompressed. It should be in the data directory: `/home/student/data`. If not, uncompress it and proceed.
``` bash
$ cd utils

$ python digest_genome.py -r A^AGCTT -o ../data/hg19.hindIII.bed ../data/hg19.fa
```

- Ensure the singularity image is in the home directory.

```bash
$ cd
$ pwd
/home/student

$ ls -la

drwxr-xr-x 15 student student      4096 May  3 05:50 .
drwxr-xr-x  4 root    root         4096 Feb 11  2020 ..
-rw-------  1 student student     16592 Apr 27 01:38 .bash_history
-rw-r--r--  1 student student       220 Feb 11  2020 .bash_logout
-rw-r--r--  1 student student      4258 Feb 11  2020 .bashrc
drwxrwxr-x  4 student student      4096 Apr 23 07:01 .cache
drwxrwxr-x  2 student student      4096 Feb 11  2020 .conda
drwxrwxr-x  4 student student      4096 Feb 17  2020 .config
drwxr-xr-x  3 student student      4096 Feb 17  2020 .local
-rw-r--r--  1 student student       807 Feb 11  2020 .profile
drwxr-xr-x 14 student student      4096 Feb 17  2020 .rstudio
-rw-rw-r--  1 student student        75 Feb 17  2020 .selected_editor
-rw-------  1 student student     12738 Apr 27 01:25 .viminfo
-rw-rw-r--  1 student student       180 May  3 05:43 .wget-hsts
-rw-rw-r--  1 student student  35351076 Mar  9  2016 HiCPro_testdata.tar.gz
-rw-r--r--  1 root    root     30460361 Feb 17  2020 Juicebox_1.11.08.jar
drwxrwxr-x  3 student student      4096 Feb 11  2020 R
drwxr-xr-x  2 student student      4096 Mar 12 00:24 alignment_prac
drwxrwxr-x  2 student student      4096 Feb 17  2020 conda.d
-rw-r--r--  1 student student      2793 Apr 27 01:18 config_test_latest.txt
drwxrwxr-x  2 student student      4096 May  3 05:58 data
drwxrwxr-x  2 student student      4096 Apr 27 01:16 hicpro
-rwxr-xr-x  1 root    root    838254623 Feb 12  2020 hicpro.img
-rw-rw-r--  1 student student 692744223 Apr 27  2020 hicpro_2.11.4_ubuntu.img # <-- What we're looking for.
-rwxr-xr-x  1 student student      3189 May  2  2012 make_hg19.sh
drwxrwxr-x 15 student student      4096 Feb 11  2020 miniconda3
drwxrwxr-x  2 student student      4096 Apr 26 22:57 tmp
drwxrwxr-x  2 student student      4096 May  3 05:43 utils
```

- Unzip the Bowtie 2 index folder. 
``` bash
$ cd data

$ unzip hg19.zip

$ ls -la

drwxrwxr-x  2 student student       4096 May  3 06:08 .
drwxr-xr-x 15 student student       4096 May  3 05:50 ..
-rw-r--r--  1 student student         26 Feb 17  2020 README
-rw-r--r--  1 student student  960018873 May  2  2012 hg19.1.bt2 # <-- What we want
-rw-r--r--  1 student student  716863572 May  2  2012 hg19.2.bt2 # <-- What we want
-rw-r--r--  1 student student       3833 May  2  2012 hg19.3.bt2 # <-- What we want
-rw-r--r--  1 student student  716863565 May  2  2012 hg19.4.bt2 # <-- What we want
-rw-r--r--  1 student student       1971 Apr 24 00:40 hg19.chrom.sizes
-rw-r--r--  1 student student 3199905909 Apr 24 00:47 hg19.fa
-rw-r--r--  1 student student  948731419 Apr 24 00:47 hg19.fa.gz
-rw-rw-r--  1 student student   36459776 May  3 05:58 hg19.hindIII.bed
-rw-r--r--  1 student student  960018873 May  3  2012 hg19.rev.1.bt2 # <-- What we want
-rw-r--r--  1 student student  716863572 May  3  2012 hg19.rev.2.bt2 # <-- What we want
-rw-r--r--  1 student student 3694403333 Apr 24 02:19 hg19.zip
-rwxr-xr-x  1 student student       3189 May  2  2012 make_hg19.sh
```
- Generate the chrom.sizes files. Open a text editor and copy the following script into it and save it as `get_chrom_sizes.py` 
``` python
import sys
from Bio import SeqIO

for record in SeqIO.parse(str(sys.argv[1]), 'fasta'):
    print(str(record.id)+'\t'+str(len(record.seq)))
```
> - Then from the directory that has `hg19.fa`, run `python get_chrom_sizes.py > hg19.chrom.sizes`

- Uncompress the test data
``` bash
$ cd 
$ tar -xvzf HiCPro_testdata.tar.gz
```

We should now have everything we need to get HiC-Pro to run. To make it easier on ourselves we will make a new directory for the HiC-Pro run and move all necessary files there.

``` bash
$ cd
$ pwd
/home/student

$ mkdir hicpro-run && cd hicpro-run

$ mv /home/student/data/hg19.hindIII.bed ./
$ mv /home/student/data/hg19.chrom.sizes ./
$ ln -s /home/student/data/hg19.1.bt2 ./
$ ln -s /home/student/data/hg19.2.bt2 ./
$ ln -s /home/student/data/hg19.3.bt2 ./
$ ln -s /home/student/data/hg19.4.bt2 ./
$ ln -s /home/student/data/hg19.rev.1.bt2 ./
$ ln -s /home/student/data/hg19.rev.2.bt2 ./

$ ls -la

drwxrwxr-x  2 student student     4096 May  3 06:33 .
drwxr-xr-x 15 student student     4096 May  3 06:31 ..
lrwxrwxrwx  1 student student       29 May  3 06:32 hg19.1.bt2 -> /home/student/data/hg19.1.bt2
lrwxrwxrwx  1 student student       29 May  3 06:32 hg19.2.bt2 -> /home/student/data/hg19.2.bt2
lrwxrwxrwx  1 student student       29 May  3 06:32 hg19.3.bt2 -> /home/student/data/hg19.3.bt2
lrwxrwxrwx  1 student student       29 May  3 06:33 hg19.4.bt2 -> /home/student/data/hg19.4.bt2
-rw-r--r--  1 student student     1971 Apr 24 00:40 hg19.chrom.sizes
-rw-rw-r--  1 student student 36459776 May  3 05:58 hg19.hindIII.bed
lrwxrwxrwx  1 student student       33 May  3 06:33 hg19.rev.1.bt2 -> /home/student/data/hg19.rev.1.bt2
lrwxrwxrwx  1 student student       33 May  3 06:33 hg19.rev.2.bt2 -> /home/student/data/hg19.rev.2.bt2
```
## 4. Set Up the config file.
``` bash
$ cd 

$ vim config_test_latest.txt # Or your preferred text editor.
```
Now we can set up our config file correctly. Lay it out as follows:
``` bash
# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM - PBS - Start Editing Here !!
#######################################################################
N_CPU = 2
LOGFILE = hicpro.log

JOB_NAME = IMR90_split
JOB_MEM = 10gb
JOB_WALLTIME = 6:00:00
JOB_QUEUE = batch
JOB_MAIL = nservant@curie.fr

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _R1
PAIR2_EXT = _R2

#######################################################################
## Alignment options
#######################################################################

FORMAT = phred33
MIN_MAPQ = 0

BOWTIE2_IDX_PATH = /home/student/hicpro-run # <-- Important
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = hg19 # <-- Important
GENOME_SIZE = /home/student/hicpro-run/hg19.chrom.sizes # <-- Important

#######################################################################
## Allele specific
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /home/student/hicpro-run/hg19.hindIII.bed # <-- Important
LIGATION_SITE = AAGCTAGCTT
MIN_FRAG_SIZE = 100
MAX_FRAG_SIZE = 100000
MIN_INSERT_SIZE = 100
MAX_INSERT_SIZE = 600

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 1
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## ICE Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
```

## 5. Run HiC-Pro
Now we can run HiC-Pro
``` bash
$ singularity shell --bind /home/student:/mnt hicpro_2.11.4_ubuntu.img

# Assuming this has configured correctly, ls -la /mnt and ls -la ~ should return the same output.

$ HiC-Pro -c /mnt/config_test_latest.txt -i /mnt/test_data/ -o /mnt/test_out
```

### I know many of you got errors at this point. Please email me the output as I have been **unable** to recreate it on my VM.  

Assuming this has worked for you proceed to the next step.

## 6. Convert HiC-Pro output to FitHiC input.
``` bash
$ cd utils
$ ls -la
drwxrwxr-x  2 student student     4096 May  3 05:43 .
drwxr-xr-x 17 student student     4096 May  3 06:51 ..
-rw-rw-r--  1 student student     6843 May  3 05:43 digest_genome.py
-rw-rw-r--  1 student student     5519 May  3 05:43 hicpro2fithic.py # <-- What we'll use.
-rw-rw-r--  1 student student     5740 May  3 05:43 hicpro2juicebox.sh
-rw-rw-r--  1 student student 36139200 Mar  5  2020 juicer_tools_1.19.02.jar

$ cd
$ pwd
/home/student

$ mkdir fithic_input && cd fithic_input

$ sample=("dixon_2M" "dixon_2M_2")

$ res=("500000" "1000000")

$ for i in ${sample[@]}; do for j in ${res[@]}; do mkdir ${i}_${j}; python ../utils/hicpro2fithic.py -i ../test_out/hic_results/matrix/${i}/raw/${j}/*.matrix -b ../test_out/hic_results/matrix/${i}/raw/${j}/*bed -s ../test_out/hic_results/matrix/${i}/iced/${j}/*biases -o ${i}_${j} -r ${j}; done; done
```

## 7. Run FitHiC
```bash
$ for j in ${res[@]}; do for i in dixon*${j}/; do name=$(basename $i); fithic -i ${i}/fithic.interactionCounts.gz -f ${i}/fithic.fragmentMappability.gz -r ${j} -t ${i}/fithic.biases.gz -o ./ -l ${name}; done; done
```

At this point you may get an error saying `ValueError: max() arg is an empty sequence`. If so, this is is most likely due to the presence of unplaced scaffolds in the fithic files. To remove these we need to clean up our input files. Using each of the sed commands below in a for loop we can clean them up nicely.
- Note: we are not running these individually so only copy and paste the for loop.
>- `sed -r -n -e '/^chr[X0-9]{1,2}\t/p' fithic.biases > clean_fithic.biases`
>- `sed -r -n -e '/^chr[X0-9]{1,2}\t/p' fithic.fragmentMappability > clean_fragmentMappability`
>- `sed '/gl/d' fithic.interactionCounts > clean_fithic.interactionCounts`

 Rather than do this individually, we'll loop this cleaning process with:
``` bash
for i in ./*; do gunzip $i/*; sed -r -n -e '/^chr[X0-9]{1,2}\t/p' $i/fithic.biases > $i/clean_fithic.biases; sed -r -n -e '/^chr[X0-9]{1,2}\t/p' $i/fithic.fragmentMappability > $i/clean_fithic.fragmentMappability; sed '/gl/d' $i/fithic.interactionCounts > $i/clean_fithic.interactionCounts; gzip $i/* ; done
```

Now when we run FitHiC it should work with no issues.
``` bash
for j in ${res[@]}; do for i in dixon*${j}/; do name=$(basename $i); fithic -i ${i}/clean_fithic.interactionCounts.gz -f ${i}/clean_fithic.fragmentMappability.gz -r ${j} -t ${i}/clean_fithic.biases.gz -o ./ -l ${name}; done; done
```

You may get some errors saying chromosome chrY not found. Don't worry about that. It's just saying the mitochondrial genome is not present in the bias file. We're not interested in the mitochondria here.

The files of interest will end in `*.significances.txt.gz`

## 8. Convert HiC-Pro to `.hic` for Juicebox Visualisation
``` bash
$ cd
$ mkdir juicebox_input
$ bash utils/hicpro2juicebox.sh -i test_out/hic_results/data/dixon_2M/dixon_2M.allValidPairs -g hicpro-run/hg19.chrom.sizes -j utils/juicer_tools_1.19.02.jar -r hicpro-run/hg19.hindIII.bed -o juicebox_input/
```

This example will be too low resolution to see any of the really cool stuff. Have a look at the Juicebox web-browser version at <https://aidenlab.org/juicebox/>