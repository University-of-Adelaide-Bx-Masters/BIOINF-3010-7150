# Practicals: Ancient DNA and population genomics

Bastien Llamas \(bastien.llamas@adelaide.edu.au\)
2020-05-26 and 2020-05-29

---
The two tutorials will be separated into:
1. Data handling (Tuesday 2020-05-26)
2. Population genomics applications (Friday 2020-05-29)

Icons are used to highlight sections of the tutorials:

:blue_book: Information

:computer: Hands-on tasks

:question: Questions

---
# Day 2: Population genomics applications (Friday 2020-05-29)

## Tutorial outcomes

:blue_book: At the end of today's tutorial, you will know how to explore contemporary and ancient genomic diversity to infer population history. The tutorial is loosely based on yesterday's lecture about the population history of Indigenous peoples of the Americas, in particular the [Posth *et al.*](https://www.sciencedirect.com/science/article/pii/S0092867418313801) *Cell* paper we published in 2018.

## Reconstructing the Deep Population History of Central and South America [(Posth *et al.* 2018 Cell)](https://www.sciencedirect.com/science/article/pii/S0092867418313801)
:blue_book: For this study, we generated 49 new genome-wide datasets that consist of an enrichment and Illumina high-throughput sequencing of 1.2M SNPs from ancient DNA samples. All sampled individuals are from Central (Belize) and South (Brazil, Peru, Chile, Argentina) American individuals. The skeletal remains are dated between 10,900–700 BP (years before present), with the large majority of remains older than 1,000 BP.

:blue_book: We know from previous genetic studies that Indigenous peoples of the Americas were isolated from the rest of the world since the peopling of the Americas until European colonisation during the 16<sup>th</sup> century. Thus we can safely assume that our ancient individual genomic datasets should not harbour signs of recent genetic admixture with non-Indigenous Americans.

:blue_book: Finally, previous work on two particular ancient individuals from North America informed us about the expected ancestry in Central and South America:
* Anzick-1: Skeletal remains of a young boy who lived ~12,800 years ago in what is currently Montana, USA. Anzick-1 belonged to a population ancestral to all contemporary Indigenous Americans [(Rasmussen *et al.* 2014 Nature)](https://www.nature.com/articles/nature13025).
* USR1: Skeletal remains of an infant recovered at Upward Sun River in what is currently Alaska, USA. The remains have been dated to ~11,400 BP. USR1 represents an ancient Beringian population that diverged from the ancestors of Indigenous Americans (to which Anzick-1 belongs) ~20,000 BP [(Moreno-Mayar *et al.* 2018 Nature)](https://www.nature.com/articles/nature25173). Therefore, USR1 is representative of an outgroup population for all Central and South Indigenous Americans.

:blue_book: The research question we asked was relatively simple: what was the processus for the earliest population movements in South America? In particular, we were interested in determining if all contemporary populations descend from one migration event, or successive waves separated in time.


## Let's explore the datasets

:computer: Activate the environment `variation`.
```bash
conda activate variation
```

:blue_book: I provided 2 datasets in `EIGENSTRAT` format. As a reminder, the `EIGENSTRAT` format consists of 3 files:
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

:computer: Unarchive the tutorial data (stored in `~/data/genomics/ancient/`) in your working directory.
```bash
mkdir -p ~/BIOINF_Friday
cd ~/BIOINF_Friday
tar xvzf ~/data/genomics/ancient/tutorial_friday.tar.gz -C .
ll
```

:blue_book: The files with the `AllAmerica_Ancient.eigenstrat` prefix contain the 49 ancient Central and South Americans and  data from modern South American individuals, as well as Anzick-1 and USR1. The files with the `AllAmerica_Ancient_YRI.eigenstrat` prefix contain the same data with the addition of one African individual (Yoruba). 

:computer: Using bash commands, answer the following questions.

---
#### :question: *Questions*
1. How many individuals are in the `AllAmerica_Ancient.eigenstrat.eigenstrat.ind` dataset?
2. Is there missing data in the ancient dataset `AllAmerica_Ancient.eigenstrat.geno`?
3. How many SNPs in each dataset? Hint: look at the `.snp` files.
---

## PCA

:computer: We are going to use `SMARTPCA` (as part of the `EIGENSOFT` utilities, see the end of Tuesday's tutorial) and an implementation of `ADMIXTOOLS` in an `R` package called admixR (it needs `ADMIXTOOLS` to be already installed). OPTIONAL: If any of the programs are not installed, do this:
```bash
conda install -c bioconda eigensoft
conda install -c bioconda admixtools
```
:computer: OPTIONAL: Install `R` packages if they are not readily available (use the `R` console).
```R
install.packages("admixr")
install.packages("tidyverse")
```

:computer: Build a parameter file named `par.AllAmerica_Ancient.smartpca` that will be one of the inputs for [SMARTPCA](https://github.com/DReichLab/EIG/tree/master/POPGEN). Because ancient data contain a lot of missing data, we are going to force `SMARTPCA` to construct the eigenvectors based on the contemporary populations (listed in [`poplistname`](https://github.com/DReichLab/EIG/tree/master/POPGEN)) and then project the ancient samples onto the PCA ([`lsqproject`](https://github.com/DReichLab/EIG/blob/master/POPGEN/lsqproject.pdf)).
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

:computer: `SMARTPCA` has generated two output files with the suffixes `.evec` (first row is the eigenvalues for the first 5 PCs, and all further rows contain the PC coordinates for each sample) and `.evac` (all the eigenvalues). Go to the `R` console and create PCA plots.
```R
library(stringr)
library(ggplot2)
library(cowplot)

# Set your working directory
setwd("~/BIOINF_Friday")

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
                  rel_widths = c(1, .5)) # ratio between plots is 1:0.5

# Prepare legend for the PCA plot
legend <- get_legend(adat.pc12 + 
                     guides(color = guide_legend(nrow = 5)) + # legend spans 4 lines
                     theme(legend.position = "bottom")) # legend is displayed at the bottom of the plots (i.e., horizontal not vertical)
                     
# Combine plots and legend using plot_grid from cowplot
plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .3)) # ratio between plots and legend is 1:0.3
```

:blue_book: A few observations:
* The contemporary Peruvians from the 1kGP (PEL) show the biggest diversity, possibly due to admixture with non-Indigenous American groups.
* All ancient samples cluster with contemporary South Americans. 


## *F*3 statistics

:blue_book: Using Eigensoft to compute *F* and *D* statistics can be very time consuming because the programs are not user friendly. Instead, we can use the `R` implementation [`admixr`](https://github.com/bodkan/admixr) by Martin Petr (article [here](https://academic.oup.com/bioinformatics/article/35/17/3194/5298728)). There is a comprehensive [tutorial](https://bodkan.net/admixr/articles/tutorial.html) that you can explore on your own time. 

:blue_book: We want to infer the relative divergence times between pairs of populations, and in which order they split of from each other. We can use an outgroup *F*3 statistic by fixing the outgroup as YRI, and calculating pairwise *F*3 statistics between populations. The higher the *F*3 value, the more shared drift between the two test populations, i.e. the more related they are.

:computer: Load the libraries needed to run `admixr` (use the `R` console) and run *F*3 statistics on a subset of populations.
```R
library(admixr)
library(tidyverse)

:computer: 
# Load the dataset that includes the African individual
snpsAmerica <- eigenstrat(prefix = "AllAmerica_Ancient_YRI.eigenstrat")

# Create a list of population we want to test (just a subset of the 
pops <- c("Surui", "Aymara", "Anzick", "USR1", "Peru_Lauricocha_5800BP", "Brazil_LapaDoSanto_9600BP", "Chile_LosRieles_10900BP")
result <- f3(A = pops, B = pops, C = "YRI", data = snpsAmerica)

head(result)
```

:blue_book: The output table contains a lot of information that we can unpack:
* `F3`: *D* statistic value
* `stderr`: standard error of the *D* statistic calculated using the block jackknife
* `Zscore`: Z-score value, which is the number of standard errors the *F*3 is from 0 (i.e. how strongly do we reject the null hypothesis of no admixture)
* `nsnps`: number of SNPs used
 
:computer: We can also represent the table into a heatmap to better visualise the pairwise comparisons:
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


## *D* statistics

:blue_book: Now we want to know how populations compare in terms of ancestry from Anzick-1. For this, we can consider of either Anzick-1 or USR1 with the different populations, using YRI as an outgroup. Since we know already that USR1 did not contribute any ancestry to South Americans, We are basically testing the proportion of Anzick-1 ancestry only.

:computer: Run *D* statistics on a subset of populations.
```R
# Load the dataset that includes the African individual
snpsAmerica <- eigenstrat(prefix = "AllAmerica_Ancient_YRI.eigenstrat")

# Create a list of population we want to test (just a subset of the 
pops <- c("Surui", "Aymara", "Peru_Lauricocha_5800BP", "Brazil_LapaDoSanto_9600BP", "Chile_LosRieles_10900BP")
result <- d(W = pops, X = "USR1", Y = "Anzick", Z = "YRI", data = snpsAmerica)

head(result)
```

:blue_book: Again, the output table contains a lot of information that we can unpack:
* `D`: *D* statistic value
* `stderr`: standard error of the *D* statistic calculated using the block jackknife
* `Zscore`: Z-score value, which is the number of standard errors the *D* is from 0 (i.e. how strongly do we reject the null hypothesis of no admixture)
* `BABA`, `ABBA`: counts of observed site patterns
* `nsnps`: number of SNPs used

:computer: However, a graphic representation is always easier to interpret:
```R
# Sort the population labels according to an increasing D value and plot average and Z-score
ggplot(result, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr))
```

:blue_book: We can confidently say that Anzick-1 contributed ancestry to all South Americans. We also observe a trend whereby the older the sample, the more Anzick-1 ancestry is contains.


