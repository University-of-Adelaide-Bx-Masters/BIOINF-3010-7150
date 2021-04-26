# BIOINF 7150 Hi-C data analysis practical

## Updated 26 April 2021 by Callum McPhillamy

Before we start, we need to download some files.

login to your vm.

``` bash
conda create --name hic
conda activate hic
conda install -c bioconda biopython numpy
pip install fithic
cd
mkdir utils
cd utils

wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.19.02.jar

wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/hicpro2fithic.py

wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/hicpro2juicebox.sh

wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/digest_genome.py

cd

wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata.tar.gz

wget http://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/hicpro_2.11.4_ubuntu.img

cd
```

In this practical, we will learn the standard procedure of handling Hi-C data by processing them from raw sequencing data to statistically significant interactions.

First of all, watch this video to familiarize with the biological background of 3D genome: how does DNA loops are formed: <https://www.youtube.com/watch?v=dES-ozV65u4>

So now you know what the nature of DNA folding.

This DNA folding mechanism can further lead to different biological events, for example, DNA loops are formed and regions that located far away from each other can have close proximity three dimensionally. More importantly, such event can affect gene regulation if the DNA loops bring distal enhancer and promoter into close proximity.

Nowadays, the most popular way to investigate the 3D genome is by analysing Hi-C data.

Hi-C stands for high resolution chromosome conformation capture assay, it was invented in 2009 (<https://www.ncbi.nlm.nih.gov/pubmed/19815776>). The procedure to prepare a Hi-C library will be covered in this week's lecture. ![HiC-prep](https://science.sciencemag.org/content/326/5950/289/F1.large.jpg)

In this practical, we will learn a standard workflow to analyse Hi-C data from raw sequencing reads to significant interactions, see Figure 1\. We will use one of the most popular pipeline called HiC-Pro (<https://github.com/nservant/HiC-Pro>) to analyse a small dataset from Rao et al. 2014\. We'll also learn how to use a tool called Juicebox to visualise the interaction matrix and use Fit-Hi-C2 to identify significant interactions.

![HiC-pipeline](https://raw.githubusercontent.com/UofABioinformaticsHub/genomics_applications/master/images/hic-pipeline.png) Figure 1: Standard workflow for analysing Hi-C data.

1. in `~/`, make sure `HiCPro_testdata.tar.gz` exists.
2. Generate files required by HiC-Pro  

**Genome Fragment file:**

```bash
cd 

cd utils

python digest_genome.py -r A^AGCTT -o <outfile_path.bed> <input_reference_genome.fa>
```
E.g. `python digest_genome.py -r A^AGCTT -o ./data/hg19.hindIII.bed ./data/hg19.fa`

**Singularity image**
- The hicpro singularity image, `hicpro_2.11.4_ubuntu.img`, should already be in the home directory.

**Reference genome bowtie index**
- The bowtie2 files will be in `~/data/hg19.zip`.
- `unzip hg19.zip` to extract the files.

**Chromosome sizes file**
- To generate the `chrom.sizes` file, open a text editor and type the following script to get the length of each chromosome in the hg19 reference file.
``` python
import sys # Allows Python to take command line arguments.
from Bio import SeqIO # Imports the SeqIO module from BioPython. This makes for super easy fasta file parsing.

for record in SeqIO.parse(str(sys.argv[1]), 'parse'): # Iterate over the fasta file, printing the record ID and the length of each fasta record
  print(str(record.id)+'\t'+str(record.seq))
```
 Save the file as `get_chrom_sizes.py`.  
 Then run `python get_chrom_sizes.py hg19.fa > hg19.chrom.sizes`.


3. Recall, the `utils` directory is where we downloaded some scripts to at the start of the prac. Inside this directory, double check these files are inside it:

```bash
cd

cd utils

ls
digest_genome.py
hicpro2fithic.py
hicpro2juicebox.sh
juicer_tools_1.19.02.jar
```
## Check out HiC-Pro

Here we used the singularity container to use HiC-Pro.

To shell into the HiC-Pro image, type `singularity shell hicpro_2.11.4_ubuntu.img` into your terminal

you should see `Singularity hicpro_2.11.4_ubuntu.img:~>`

Try `HiC-Pro -h`, see if you can access the help page for HiC-Pro.

Based on the help page, we can see that there are three required parameters for running HiC-Pro, `-i|--input INPUT`, `-o|--output OUTPUT` and `-c|--conf CONFIG`.

## Explain and edit the config file

The configure file is specific design for HiC-Pro to recognise all the parameters it need to run the pipeline.

First we get the configure file by doing:

```bash
tar -zxvf HiCPro_testdata.tar.gz
less config_test_latest.txt
```

Now we need to edit the configure file so HiC-Pro know where it can find all the files it needs, including:

1. reference genome
2. genome size file
3. restriction fragment file

## Run HiC-Pro and explain every steps of it

In order to run HiC-Pro, first get into the singularity shell, but this time, we need files from local so we need to bind the files to it by doing:

```
singularity shell --bind /home/student:/mnt hicpro_2.11.3_ubuntu.img
```

And when we use the `--bind` parameter, we can find all the local files in `/mnt`, try `ls -lh /mnt`.

Now we can run HiC-Pro with the edited configure file.

```
HiC-Pro -c /mnt/config_test_latest.txt -i /mnt/test_data/ -o /mnt/test_out
```

![hic-pro process](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-015-0831-x/MediaObjects/13059_2015_831_Fig4_HTML.gif?as=webp)

## Inspect the results, answer these questions

1. What's the uniquely mapping rates of of dixon_2M and dixon_2M_2?

2. Are all the uniquely mapping read pairs become valid interactions? If not, how many percentage of it becomes valid interactions?

3. What kinds of uninformative reads have been removed, please list at least three kinds? What's the biggest contamination?

4. How many binned interactions have over 5 read pairs mapped to the bins at each resolution?

## Identify significant interactions using FitHiC

- Why do we need to identify significant interactions?

Before we use FithiC, we need to check one thing by:

```
grep 'ningbioinfostruggling' miniconda3/envs/hic/lib/python3.6/site-packages/fithic/fithic.py
```

Source to hic environment and check out the help page of Fithic

```
source activate miniconda3/envs/hic/
fithic -h
```

To run Fithic, first we need to get the input for fithic, and it required python 2.7 while now we are in the environment of python 3

```
python2 utils/hicpro2fithic.py -h
```

This script can link the output of HiC-Pro to FitHiC

so we have 4 target to identify significant interactions, to avoid writing the command over and over again, we use loops:

```
mkdir fithic_input

cd fithic_input

sample=("dixon_2M" "dixon_2M_2")

res=("500000" 1000000)

for i in ${sample[@]}; do for j in ${res[@]}; do mkdir ${i}_${j}; python2 ../utils/hicpro2fithic.py -i ../test_out/hic_results/matrix/${i}/raw/${j}/*.matrix -b ../test_out/hic_results/matrix/${i}/raw/${j}/*bed -s ../test_out/hic_results/matrix/${i}/iced/${j}/*biases -o ${i}_${j} -r ${j}; done; done
```

Now we can again source to hic environment.

```
for j in ${res[@]}; do for i in dixon*${j}/; do name=$(basename $i); fithic -i ${i}/fithic.interactionCounts.gz -f ${i}/fithic.fragmentMappability.gz -r ${j} -t ${i}/fithic.biases.gz -o ./ -l ${name}; done; done
```

At this point you may get an error along the lines of `ValueError: max() arg is an empty sequence`. If so, this is is most likely due to the presence of unplaced scaffolds in the fithic files. To remove these we need to clean up our input files. Using each of the sed commands below in a for loop we can clean them up nicely.
- `sed -r -n -e '/^chr[X0-9]{1,2}\t/p' fithic.biases > clean_fithic.biases`
- `sed -r -n -e '/^chr[X0-9]{1,2}\t/p' fithic.fragmentMappability > clean_fragmentMappability`
- `sed '/gl/d' fithic.interactionCounts > clean_fithic.interactionCounts`
- Rather than do this individually, we'll loop this cleaning process with:
`for i in ./*; do gunzip $i/*; sed -r -n -e '/^chr[X0-9]{1,2}\t/p' $i/fithic.biases > $i/clean_fithic.biases; sed -r -n -e '/^chr[X0-9]{1,2}\t/p' $i/fithic.fragmentMappability > $i/clean_fithic.fragmentMappability; sed '/gl/d' $i/fithic.interactionCounts > $i/clean_fithic.interactionCounts; gzip $i/* ; done`

Now when we run FitHiC it should work with no issues.
for j in ${res[@]}; do for i in dixon*${j}/; do name=$(basename $i); fithic -i ${i}/clean_fithic.interactionCounts.gz -f ${i}/clean_fithic.fragmentMappability.gz -r ${j} -t ${i}/clean_fithic.biases.gz -o ./ -l ${name}; done; done


1. At resolution of 500kb, how many significant interactions are identified with p-value <= 0.05?
2. Think about what we can use with this interactions, i.e. how to give them biological meaning?

## Visualise HiC interactions using heatmaps, introduce TADs and loops

Lastly we'll try to visualise the interaction matrix.

Take dixon_2M at 500000 as example.

Like FitHiC, first we need to make the output of HiC-Pro compatible.

```
cd

mkdir juicebox_input

bash utils/hicpro2juicebox.sh -h

bash utils/hicpro2juicebox.sh -i test_out/hic_results/data/dixon_2M/dixon_2M.allValidPairs -g hg19-ref/chrom_hg19.sizes -j utils/juicer_tools_1.19.02.jar -r hg19-ref/HindIII_resfrag_hg19.bed -o juicebox_input/
```

Show the Juicebox web-browser version at <https://aidenlab.org/juicebox/>
