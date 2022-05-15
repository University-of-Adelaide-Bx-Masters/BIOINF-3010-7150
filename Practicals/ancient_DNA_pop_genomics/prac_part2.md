# Practicals: Ancient DNA and population genomics

Bastien Llamas \(bastien.llamas@adelaide.edu.au\)
2022-05-24 and 2022-05-27

---
The two practicals will be separated into:
1. Data handling (Tuesday 2022-05-24)
2. Population genomics applications (Friday 2022-05-27)

Icons are used to highlight sections of the practicals:

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Information  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Hands-on tasks  
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> Questions 

---
# Day 2: Population genomics applications (Friday 2022-05-27)

## Practical outcomes

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> At the end of today's practical, you will know how to explore contemporary and ancient genomic diversity to infer population history. The practical is loosely based on Friday's lecture about the population history of Indigenous peoples of the Americas, in particular the [Posth *et al.*](https://www.sciencedirect.com/science/article/pii/S0092867418313801) *Cell* paper we published in 2018.

## Reconstructing the Deep Population History of Central and South America [(Posth *et al.* 2018 Cell)](https://www.sciencedirect.com/science/article/pii/S0092867418313801)
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> For this study, we generated 49 new genome-wide datasets that consist of an enrichment and Illumina high-throughput sequencing of 1.2M SNPs from ancient DNA samples. All sampled individuals are from Central (Belize) and South (Brazil, Peru, Chile, Argentina) American individuals. The skeletal remains are dated between 10,900–700 BP (years before present), with the large majority of remains older than 1,000 BP.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> We know from previous genetic studies that Indigenous peoples of the Americas were isolated from the rest of the world since the peopling of the Americas until European colonisation during the 16<sup>th</sup> century. Thus we can safely assume that our ancient individual genomic datasets should not harbour signs of recent genetic admixture with non-Indigenous Americans.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Finally, previous work on two particular ancient individuals from North America informed us about the expected ancestry in Central and South America:
* Anzick-1: Skeletal remains of a young boy who lived ~12,800 years ago in what is currently Montana, USA. Anzick-1 belonged to a population ancestral to all contemporary Indigenous Americans [(Rasmussen *et al.* 2014 Nature)](https://www.nature.com/articles/nature13025).
* USR1: Skeletal remains of an infant recovered at Upward Sun River in what is currently Alaska, USA. The remains have been dated to ~11,400 BP. USR1 represents an ancient Beringian population that diverged from the ancestors of Indigenous Americans (to which Anzick-1 belongs) ~20,000 BP [(Moreno-Mayar *et al.* 2018 Nature)](https://www.nature.com/articles/nature25173). Therefore, USR1 is representative of an outgroup population for all Central and South Indigenous Americans.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The research question we asked was relatively simple: what was the processus for the earliest population movements in South America? In particular, we were interested in determining if all contemporary populations descend from one migration event, or successive waves separated in time.


## Let's explore the datasets

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Activate the environment `variation`.
```bash
conda activate variation
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

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Unarchive the practical data (stored in `~/data/genomics/ancient/`) in your working directory.
```bash
mkdir -p ~/BIOINF_Tuesday
cd ~/BIOINF_Tuesday
tar xvzf ~/data/genomics/ancient/practical_tuesday.tar.gz -C .
ll
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The files with the `AllAmerica_Ancient.eigenstrat` prefix contain the 49 ancient Central and South Americans and  data from modern South American individuals, as well as Anzick-1 and USR1. The files with the `AllAmerica_Ancient_YRI.eigenstrat` prefix contain the same data with the addition of one African individual (YRI: Yoruba). 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Using bash commands, answer the following questions.

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*
- 1) How many individuals are in the `AllAmerica_Ancient.eigenstrat.ind` dataset?
- 2) Is there missing data in the ancient dataset `AllAmerica_Ancient.eigenstrat.geno`?
- 3) How many SNPs in each dataset? Hint: look at the `.snp` files.  

---

## PCA

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We are going to use `SMARTPCA` (as part of the `EIGENSOFT` utilities, see the end of Tuesday's practical) and an implementation of `ADMIXTOOLS` in an `R` package called admixR (it needs `ADMIXTOOLS` to be already installed). OPTIONAL: If any of the programs are not installed, do this:
```bash
conda install -c bioconda eigensoft
conda install -c bioconda admixtools
echo "PATH=$PATH" > ~/.Renviron
```
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Go to Session > Restart R. OPTIONAL: Install `R` packages if they are not readily available (use the `R` console).

```R
library(devtools)
devtools::install_github("bodkan/admixr")
install.packages("tidyverse")
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Build a parameter file named `par.AllAmerica_Ancient.smartpca` that will be one of the inputs for [SMARTPCA](https://github.com/DReichLab/EIG/tree/master/POPGEN). Because ancient data contain a lot of missing data, we are going to force `SMARTPCA` to construct the eigenvectors based on the contemporary populations (listed in [`poplistname`](https://github.com/DReichLab/EIG/tree/master/POPGEN)) and then project the ancient samples onto the PCA ([`lsqproject`](https://github.com/DReichLab/EIG/blob/master/POPGEN/lsqproject.pdf)).  

```bash
genotypename:    AllAmerica_Ancient.eigenstrat.geno
snpname:         AllAmerica_Ancient.eigenstrat.snp
indivname:       AllAmerica_Ancient.eigenstrat.ind
evecoutname:     AllAmerica_Ancient.smartpca_results.evec
evaloutname:     AllAmerica_Ancient.smartpca_results.eval
numoutevec:      5
lsqproject:      YES
poplistname:     poplistPCA
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Run `SMARTPCA`

```bash
smartpca -p par.AllAmerica_Ancient.smartpca
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> `SMARTPCA` has generated two output files with the suffixes `.evec` (first row is the eigenvalues for the first 5 PCs, and all further rows contain the PC coordinates for each sample) and `.evac` (all the eigenvalues). Go to the `R` console and create plots.

```R
library(stringr)
library(ggplot2)
library(cowplot)

# Set your working directory
setwd("~/BIOINF_Tuesday")

# data for scree plot
adat.scree <- read.table("AllAmerica_Ancient.smartpca_results.eval", header = FALSE)
# Add a column with row number (only needed to be able to do a bar plot)
adat.scree$Name = 1:nrow(adat.scree)
# Rename columns
colnames(adat.scree) <- c("Scree","Name")
# Restrict to the first 10 lines (first 10 eigenvalues)
adat.scree <- head(adat.scree, 10)
# Plot scree plot data
adat.screep <- ggplot(adat.scree,aes(x = Name, y = Scree)) +
               geom_bar(stat = "identity") + # heights of the bars represent values in the data
               theme(axis.text.x = element_blank(), # no text for the x-axis
                     axis.ticks.x = element_blank()) + # no ticks for the x-axis
               xlab("component") + # x-axis label
               ylab("eigenvalue") # y-axis label
          
# Create dataframe for the data
adat <- read.table("AllAmerica_Ancient.smartpca_results.evec", header = FALSE)
# Reduce table to just the sample name, the first 3 PCs, and the population name
adat <- adat[,c(1:4,7)]
# Rename columns
colnames(adat) <- c("SAMPLE", "PC1", "PC2", "PC3", "POP")
# Create column for ancient vs contemporary, conditioning on the string "BP" (included in the dates) or the name of the ancient sample (Anzick and USR1)
adat$DATE <- ifelse(grepl('BP', adat$POP), "Ancient",
                    ifelse(grepl('Anzick', adat$POP), "Anzick",
                           ifelse(grepl('USR', adat$SAMPLE), "USR", "Contemporary")))
# Create column with simplified names
#adat$GROUP <- word(adat$POP, 1, sep = "_|\\.")
# Create plot with populations in different colours and super-populations with different point shapes
adat.pc12 <- ggplot(adat, aes(x = PC1, y = PC2)) + 
             geom_point(aes(fill = POP, colour = POP, shape = DATE), size = 4) +
             scale_shape_manual(values=c(21, 22, 24, 25))

# Combine plots using plot_grid from cowplot
prow <- plot_grid(adat.pc12 + theme(legend.position="none"), # remove the legend
                  adat.screep + theme(legend.position="none"), # remove the legend
                  align = 'vh', # plots are aligned vertically and horizontally
                  nrow = 1, # 2 plots on one row
                  rel_widths = c(1, .5)) # width ratio between plots is 1:0.5

# Prepare legend for the PCA plot
legend <- get_legend(adat.pc12 + 
                     guides(color = guide_legend(nrow = 5)) + # legend spans 5 lines
                     theme(legend.position = "bottom")) # legend is displayed at the bottom of the plots (i.e., horizontal not vertical)
                     
# Combine plots and legend using plot_grid from cowplot
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .3)) # height ratio between plots and legend is 1:0.3
```

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*
- 4) The scree plot represents the value for each eigenvector, i.e., the variance in the data explained by the eigenvector. In your opinion, does the first eigenvector explain much variance compared to other vectors?  
- 5) PC1 seems to capture the variation observed between eskimos and modern Peruvian (PEL), while PC2 seems to capture the variation just within PEL. Knowing that PEL is individuals from Lima, the capital city of Peru, why would the PEL population be so diverse?  
- 6) Where do the ancient samples cluster in regards to the PCA coordinates? And where in regards to contemporary populations?  

---


## *F*3 statistics

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> We want to infer the relative divergence times between pairs of populations, and in which order they split from each other. We can use an outgroup *F*3 statistic by fixing the outgroup as YRI, and calculating pairwise *F*3 statistics between populations. The higher the *F*3 value, the more shared drift between the two test populations, i.e. the more related they are.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Using `ADMIXTOOLS` to compute *F* and *D* statistics can be very time consuming because the programs are not user friendly, and building the parameter files can be time consuming. Instead, we can use the `R` implementation [`admixr`](https://github.com/bodkan/admixr) of `ADMIXTOOLS` by Martin Petr. You may want to read the very short [*Bioinformatics* Applications Note](https://academic.oup.com/bioinformatics/article/35/17/3194/5298728), or better, explore the comprehensive [tutorial](https://bodkan.net/admixr/articles/tutorial.html) on your own time. 

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Load the libraries needed to run `admixr` (use the `R` console) and run *F*3 statistics on a subset of populations.

```R
library(admixr)
library(tidyverse)

# Load the dataset that includes the African individual
snpsAmerica <- eigenstrat(prefix = "AllAmerica_Ancient_YRI.eigenstrat")

# Create a list of population we want to test (just a subset of the 
pops <- c("Eskimo", "Aymara", "Anzick", "USR1", "Peru_Lauricocha_8600BP", "Peru_Lauricocha_5800BP", "Brazil_LapaDoSanto_9600BP", "Chile_LosRieles_10900BP")
result <- f3(A = pops, B = pops, C = "YRI", data = snpsAmerica)

head(result)
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> The output table contains a lot of information that we can unpack:
* `F3`: *F*3 statistic value
* `stderr`: standard error of the *F*3 statistic calculated using the block jackknife
* `Zscore`: *Z*-score value, which is the number of standard errors the *F*3 is from 0 (i.e. how strongly do we reject the null hypothesis of no admixture)
* `nsnps`: number of SNPs used
 
<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> We can also plot a heatmap to better visualise the pairwise comparisons:
```R
# Sort the population labels according to an increasing F3 value relative to Aymara
ordered <- filter(result, A == "Aymara", B != "Aymara") %>% arrange(f3) %>% .[["B"]] %>% c("Aymara")

# Plot heatmap of pairwise F3 values
result %>%
  filter(A != B) %>%
  mutate(A = factor(A, levels = ordered),
         B = factor(B, levels = ordered)) %>%
  ggplot(aes(A, B)) + geom_tile(aes(fill = f3))
```

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*  

- 7) What two populations/individuals seem to diverge earlier than the others?  

---


## *D* statistics

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Now we want to know how populations compare in terms of ancestry from Anzick-1. For this, we can consider of comparing either Anzick-1 or USR1 with the different populations, using YRI as an outgroup. Since we know already that USR1 did not contribute any ancestry to South Americans, we are basically testing the proportion of Anzick-1 ancestry only.

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> Run *D* statistics on a subset of populations.

```R
# Create a list of populations we want to test (just a subset of the populations in the dataset)
pops2 <- c("Eskimo", "Aymara", "Peru_Lauricocha_5800BP", "Brazil_LapaDoSanto_9600BP", "Chile_LosRieles_10900BP")
result2 <- d(W = pops, X = "USR1", Y = "Anzick", Z = "YRI", data = snpsAmerica)

head(result2)
```

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> Again, the output table contains a lot of information that we can unpack:
* `D`: *D* statistic value
* `stderr`: standard error of the *D* statistic calculated using the block jackknife
* `Zscore`: *Z*-score value, which is the number of standard errors the *D* is from 0 (i.e. how strongly we reject the null hypothesis of no admixture)
* `BABA`, `ABBA`: counts of observed site patterns
* `nsnps`: number of SNPs used

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/computer_black_24dp.png" alt="Computer"/> However, a graphic representation is always easier to interpret:  

```R
# Sort the population labels according to an increasing D value and plot average and Z-score
ggplot(result2, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr))
```

---
#### <img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/quiz_black_24dp.png" alt="Questions"/> *Questions*  

- 8) Is there any test population/individual for which *D* is not different from 0? What does it mean in terms of admixture?  
- 9) Is there any test population/individual for which *D* is different from 0? Any particular pattern to report?  

---

<img src="https://raw.githubusercontent.com/University-of-Adelaide-Bx-Masters/BIOINF-3010-7150/master/images/book_black_24dp.png" alt="Book"/> To estimate the minimum number of streams of ancestry contributing to Central and South American populations, we have used in our study the software `qpWave` (also implemented in `admixr`). `qpWave` assesses whether *F*4-statistics of the form *F*4(A = South American 1, B = South American 2; X = outgroup 1, Y = outgroup 2) form a matrix that is consistent with different ranks: rank 0 is consistent with a single stream of ancestry relative to the outgroups, rank 1 means 2 streams of ancestry, etc. This is how we could identify at least 3 streams of ancestry: one related to Anzick-1, 2 others related to other North American populations and never reported before.


