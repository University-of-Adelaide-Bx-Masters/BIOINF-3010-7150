# Ancient DNA and population genomics practical: Part 2 - Bastien Llamas


Icons are used to highlight sections of the practicals:

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Information  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Hands-on tasks  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> Questions 

---
# Day 2: Population genomics applications

## Practical outcomes

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> At the end of today's practical, you will know how to explore contemporary and ancient genomic diversity to infer population history. The practical is loosely based on Monday's lecture about the population history of Indigenous peoples of the Americas, in particular the [Posth *et al.*](https://www.sciencedirect.com/science/article/pii/S0092867418313801) *Cell* paper that was published in 2018.

## Reconstructing the Deep Population History of Central and South America [(Posth *et al.* 2018 Cell)](https://www.sciencedirect.com/science/article/pii/S0092867418313801)
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> This study generated 49 new genome-wide datasets that consist of enriched Illumina high-throughput sequencing of 1.2M SNPs from ancient DNA samples. All sampled individuals are from archaeological sites in Central (Belize) and South (Brazil, Peru, Chile, Argentina) America. The skeletal remains are dated between 10,900–700 BP (years before present), with the large majority of remains older than 1,000 BP.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> We know from previous genetic studies that Indigenous peoples of the Americas were isolated from the rest of the world since the peopling of the Americas until European colonisation during the 16<sup>th</sup> century. Thus we can safely assume that our ancient individual genomic datasets should not harbour signs of recent genetic admixture with non-Indigenous Americans.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Finally, previous work on two particular ancient individuals from North America informed us about the expected ancestry in Central and South America:
* Anzick-1: Skeletal remains of a young boy who lived ~12,800 years ago in what is currently Montana, USA. Anzick-1 belonged to a population ancestral to all contemporary Indigenous Americans [(Rasmussen *et al.* 2014 Nature)](https://www.nature.com/articles/nature13025).
* USR1: Skeletal remains of an infant recovered at Upward Sun River in what is currently Alaska, USA. The remains have been dated to ~11,400 BP. USR1 represents an ancient Beringian population that diverged from the ancestors of Indigenous Americans (to which Anzick-1 belongs) ~20,000 BP [(Moreno-Mayar *et al.* 2018 Nature)](https://www.nature.com/articles/nature25173). Therefore, USR1 is representative of an outgroup population for all Central and South Indigenous Americans.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The research question we asked was relatively simple: what was the processus for the earliest population movements in South America? In particular, we were interested in determining if all contemporary populations descend from one migration event, or successive waves separated in time.


## Let's explore the datasets

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/>  Activate the `popgen` environment:
```bash
source activate popgen
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> I provided 2 datasets in `EIGENSTRAT` format. As a reminder, the `EIGENSTRAT` format consists of 3 files:
* `.geno`: tab-delimited genotype file with one line per SNP and genotypes in non-separated columns, with the following genotype coding:
  * `0`: no copies of reference allele
  * `1`: one copy of reference allele
  * `2`: two copies of reference allele
  * `9`: missing data
* `.snp`: tab-delimited SNP file with one line per SNP and the following 6 columns (last 2 optional):
  * SNP name
  * Chromosome (X is encoded as 23, Y as 24, mtDNA as 90, and XY as 91)
  * Genetic position (in Morgans). 0 if unknown
  * Physical position (in bases)
  * Optional 5th and 6th columns are reference and variant alleles. For monomorphic SNPs, the variant allele can be encoded as X (unknown)
* `.ind`: tab-delimited sample file with one line per individual and the following 3 columns:
  * sample ID
  * gender (M or F). U for Unknown
  * Case or Control status, or population group label. If this entry is set to "Ignore", then that individual and all genotype data from that individual will be removed from the data set in all `CONVERTF` output.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Copy `R` scripts and unarchive the practical data (stored in `~/data/genomics/ancient/`) in your working directory.
```bash
mkdir -p ~/Project_12_2/{data,results,scripts}
cd ~/Project_12_2/
cp ~/data/ancient/prac_2/*.R ~/Project_12_2/scripts/
tar xvzf ~/data/ancient/tutorial_popgen.tar.gz -C data/
ll data/
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The files with the `AllAmerica_Ancient.eigenstrat` prefix contain the 49 ancient Central and South Americans and  data from modern South American individuals, as well as Anzick-1 and USR1. The files with the `AllAmerica_Ancient_YRI.eigenstrat` prefix contain the same data with the addition of one African individual (YRI: Yoruba). 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Using bash commands, answer the following questions.

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*<br>
Q1. How many individuals are in the `AllAmerica_Ancient.eigenstrat.ind` dataset?<br>
Q2. Is there missing data in the ancient dataset `AllAmerica_Ancient.eigenstrat.geno`?<br>
Q3. How many SNPs in each dataset? *Hint: look at the `.snp` files*<br>

<details>
  <summary>Answers</summary>
  
  Q1: 213<br>
  `wc -l data/AllAmerica_Ancient.eigenstrat.ind`<br>
  
  Q2: yes, there is a lot of `9`<br>
  `grep -c "9" data/AllAmerica_Ancient.eigenstrat.geno`<br>

  Q3: 1,196,673 SNPs in each dataset<br>
  `for i in data/*.snp; do wc -l $i; done`<br>

</details>

---

## PCA

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We are going to use `SMARTPCA` (as part of the `EIGENSOFT` utilities, see the end of Tuesday's practical) and an implementation of `ADMIXTOOLS` in an `R` package called admixR (it needs `ADMIXTOOLS` to be already installed). These programs are installed in the conda environment `popgen`.
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Make sure to copy the `$PATH` variable to `.Renviron`.
```bash
echo "PATH=$PATH" >> ~/.Renviron
```
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Go to Session > Restart R. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We will use [SMARTPCA](https://github.com/DReichLab/EIG/tree/master/POPGEN) to build a PCA of the data. Because ancient data contain a lot of missing data, we are going to force `SMARTPCA` to construct the eigenvectors based on the contemporary populations (listed in [`poplistname`](https://github.com/DReichLab/EIG/tree/master/POPGEN)) and then project the ancient samples onto the PCA ([`lsqproject`](https://github.com/DReichLab/EIG/blob/master/POPGEN/lsqproject.pdf)).  

```bash
smartpca -p <(echo "genotypename:    data/AllAmerica_Ancient.eigenstrat.geno
snpname:         data/AllAmerica_Ancient.eigenstrat.snp
indivname:       data/AllAmerica_Ancient.eigenstrat.ind
evecoutname:     results/AllAmerica_Ancient.smartpca_results.evec
evaloutname:     results/AllAmerica_Ancient.smartpca_results.eval
numoutevec:      10
lsqproject:      YES
poplistname:     data/poplistPCA")
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> `SMARTPCA` has generated two output files with the suffixes `.evec` (first row is the eigenvalues for the first 10 PCs, and all further rows contain the PC coordinates for each sample) and `.evac` (all the eigenvalues). Run the custom Rscript to generate the scree and PCA plots. The combined plot should be in a `.pdf` file in `results/`.

```bash
Rscript scripts/plot_smartpca.R
```

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*<br>

Q4. The scree plot represents the value for each eigenvector, i.e. the variance in the data explained by the eigenvector. In your opinion, does the first eigenvector explain much variance compared to other vectors?<br>
Q5. PC1 seems to capture the variation observed between eskimos and modern Peruvian (PEL), while PC2 seems to capture the variation just within PEL. Knowing that PEL is individuals from Lima, the capital city of Peru, why would the PEL population be so diverse?<br>
Q6. Where do the ancient samples cluster in regards to the PCA coordinates? And where in regards to contemporary populations?<br>

<details>
  <summary>Answers</summary>
  
  Q4: No<br>
  
  Q5: Mixed ancestry between Indigenous South Americans and Europeans (colonial history), Africans (slave trade), and East Asians (20th century migrations)<br>

  Q6: around 0-0, on top of contemporary populations<br>

</details>

---


## *F*3 statistics

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> We want to infer the relative divergence times between pairs of populations, and in which order they split from each other. We can use an outgroup *F*3 statistic by fixing the outgroup as YRI (Yoruba, in Africa), and calculating pairwise *F*3 statistics between populations. The higher the *F*3 value, the more shared drift between the two test populations, i.e. the more related they are.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Using `ADMIXTOOLS` to compute *F* and *D* statistics can be very time consuming because the programs are not user friendly, and building the parameter files is usually an opportunity for human error. Instead, we can use the `R` implementation [`admixr`](https://github.com/bodkan/admixr) of `ADMIXTOOLS` by Martin Petr. You may want to read the very short [*Bioinformatics* Applications Note](https://academic.oup.com/bioinformatics/article/35/17/3194/5298728), or better, explore the comprehensive [tutorial](https://bodkan.net/admixr/articles/tutorial.html) on your own time. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> In the `R` console (not the terminal!), set the working directory and run *F*3 statistics on a subset of populations using the script `scripts/run_F3.R`. Stop where it says to stop...

```R
setwd("~/Project_12_2/")
file.edit("scripts/run_F3.R")
#run the script line by line in the console
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The output table contains a lot of information that we can unpack:
* `A`, `B`, `C`: populations used in the test. Note that we keep YRI as the outgroup `C`
* `f3`: *F*3 statistic value
* `stderr`: standard error of the *F*3 statistic calculated using the block jackknife
* `Zscore`: *Z*-score value, which is the number of standard errors the *F*3 is from 0 (i.e. how strongly do we reject the null hypothesis of no admixture)
* `nsnps`: number of SNPs used
 
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Resume the `R` script in the `R` console and plot the results in a heatmap to better visualise the pairwise comparisons. Stop where it says to stop...

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*<br>

Q7. What two populations/individuals seem to diverge earlier than the others?<br>

<details>
  <summary>Answers</summary>
  
  Q7: Peru_Lauricocha_5800BP and Chile_LosRieles_10900BP<br>

</details>

---


## *D* statistics

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Now we want to know how populations compare in terms of ancestry from Anzick-1. For this, we can compare either Anzick-1 or USR1 with the different populations, using YRI as an outgroup. Since we know already from the paper that USR1 did not contribute any ancestry to South Americans, we are basically testing the proportion of Anzick-1 ancestry only.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Resume the `R` script in the `R` console and run *D* statistics on a subset of populations. Stop where it says to stop...

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Again, the output table contains a lot of information that we can unpack:
* `W`, `X`, `Y`, `Z`: populations used in the test. Note that we keep YRI as the outgroup `Z`, USR1 and Anzick-1 as the sources of admixture `X` and `Y`, and we rotate all other populations as `W`
* `D`: *D* statistic value
* `stderr`: standard error of the *D* statistic calculated using the block jackknife
* `Zscore`: *Z*-score value, which is the number of standard errors the *D* is from 0 (i.e. how strongly we reject the null hypothesis of no admixture)
* `BABA`, `ABBA`: counts of observed BABA and ABBA site patterns
* `nsnps`: number of SNPs used

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> However, a graphic representation is always easier to interpret. Resume the `R` script in the `R` console and plot the results.  

---
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*<br>

Q8. Is there any test population/individual for which *D* is not different from 0? What does it mean in terms of admixture?<br>
Q9. Is there any test population/individual for which *D* is different from 0? Any particular pattern to report?<br>

<details>
  <summary>Answers</summary>
  
  Q8: Yes, the Eskimo population. Anzick-1 did not contribute ancestry to the Eskimos<br>

  Q9: Yes, all ancient and contemporary South Americans. There seems to be variable amount of Anzick-1 ancestry<br>

</details>

---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> To estimate the minimum number of streams of ancestry contributing to Central and South American populations, we have used in our study the software `qpWave` (also implemented in `admixr`). `qpWave` assesses whether *F*4-statistics of the form *F*4(A = South American 1, B = South American 2; X = outgroup 1, Y = outgroup 2) form a matrix that is consistent with different ranks: rank 0 is consistent with a single stream of ancestry relative to the outgroups, rank 1 means 2 streams of ancestry, etc. This is how we could identify at least 3 streams of ancestry: one related to Anzick-1, 2 others related to other North American populations and never reported before.


