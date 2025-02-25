Assignment 7: PGS
================

-   [Assignment Overview](#assignment-overview)
-   [Getting Ready](#getting-ready)
-   [Genotyping Quality Control](#genotyping-quality-control)
    -   [General QC](#general-qc)
    -   [Global Ancestry Investigation](#global-ancestry-investigation)
    -   [a. PCA-specific QC](#a-pca-specific-qc)
    -   [b. PCA computation](#b-pca-computation)
    -   [c. Visualization](#c-visualization)
-   [Imputation](#imputation)
-   [Polygenic Scores (PGS)](#polygenic-scores-pgs)
    -   [PGS accuracy](#pgs-accuracy)
-   [Authors and contributions](#authors-and-contributions)

# Assignment Overview

In this assignment we will learn about population stratification,
imputation of genotypes, and using polygenic scores. Polygenic scores
(PGSs) can be useful for predicting disease susceptibility. In order to
calculate PGSs, we need two things: GWAS summary statistics (including
effect sizes), and genotypes. Most of the time, only a subset of a
person’s genotypes are actually measured (e.g. via SNP array), and so we
must impute the rest using a matched population of fully genotyped
individuals. This is the goal of Assignment 7.

Throughout the assignment we will be using a Mini Cohort that has
genetic data and some phenotypic variables, together with the 1000
Genomes Project samples. Both datasets are in bfile plink format, which
encompasses 3 files: *.bim, .bed and .fam* all the files can be located
under the following path: */projects/bmeg/A7*

# Getting Ready

In this assignment, you will use the plink tool extensively. A plink
tutorial can be found here:
<https://zzz.bwh.harvard.edu/plink/tutorial.shtml>

``` bash
## Install plink1.9 onto your A1 conda environment:
conda install -c bioconda plink
```

# Genotyping Quality Control

## General QC

Before we can start working on the genetic data, we need to ensure that
the quality is adequate. Thus, we are gonna check the following
measuring for our MiniCohort:

1.  **SNP call rate:** The call rate represents the percentage of
    participants with non-missing data for that SNP. Removing variants
    with a call rate lower than 95% avoids potential wrong calls to be
    included in further analysis

2.  **Minor Allele Frequency:** The minor allele frequency (MAF) echoes
    the less common allele frequency across the population. The MAF
    estimates tend to be more accurate for higher MAFs and the
    population sample size the MAF was based on. If there are too few
    samples representing the rare-allele, is hard to distinguish between
    a true rare-allele and sequencing errors.

3.  **Sample call rate:** Similar to SNP call rate, it allows to filter
    out all samples exceeding 98% missing genetic variants out of all
    the calls.

``` bash
## Using only one run of plink 1.9 (with different flags)
## 1. Filter out -SNPs- with more than 5% missingness
## 2. Filter out -variants- with less than 1% MAF
## 3. Filter out -samples- with more than 2% missingness
## 4. Create an output file in bfile format (which contains the bed, fam and bim files) for the MiniCohort QCed data

#?# Type the command you used below: - 3pt
plink --bfile /projects/bmeg/A7/Mini_cohort --geno 0.05 --maf 0.01 --mind 0.02 --make-bed --out filtered
```

## Global Ancestry Investigation

In order to enhance imputation accuracy when dealing with ethnically
diverse cohorts is important to understand the genetic ancestries of the
cohort’s participants. Knowing the ancestral populations will ensure
that the most closely related population is used as a reference for the
imputation. For instance, one would not want to impute haplotypes of an
individual of Yoruban ancestry with a population of East Asians because
many of the haplotypes will differ between the two ancestries, leading
to imputing the wrong variants for the Yoruban person. Hence, we will
analyze the global ancestry of our cohort using Principal Component
Analysis (PCA). PCA is an unsupervised, unbiased way to reduce the
complexity of multidimensional.

## a. PCA-specific QC

We first need to ensure that only the most informative genetic variants
are used in the analysis. To do this, we will:

1.  **Filter out high linkage disequilibrium (LD) regions:** Because
    high LD regions will add redundancy to the PCA (leading to these
    regions dominating top PCs), they need to be removed.

2.  **LD pruning:** Similarly, LD causes redundancy even outside the
    particularly problematic high-LD regions. Thus, we will use
    LD-pruning to identify variants that are in LD, and select one per
    block.

``` bash
## Using only one run of plink 1.9 (with different flags)
## 1. Filter out the high-LD regions contained in the --high_LD_regions_hg19.txt-- file, located in /projects/bmeg/A7/

## 2. Use the --indep-pairwise to do LD prunning with the following parameters:
## - Window size: 200, 
## - Variant Count: 100 
## - r^2: 0.2 
#?# Type the command you use to create the Mini Cohort PCA-QCed bfile below: - 1pt
plink --bfile filtered --r --ld-snp-list /projects/bmeg/A7/high_LD_regions_hg19.txt --indep-pairwise 200 100 0.2 --out filtered2


## Use the output -.prune.in- file to extract only the informative variants and create a new bfile format (bed, fam and bim files) from:
## 1. The General MiniCohort QC bfile created before
## 2. The 1KGP_reference bfile located in /projects/bmeg/A7/
#?# Type the commands you used below: - 3pt
plink --bfile filtered --extract filtered2.prune.in --make-bed --out filtered3
plink --bfile /projects/bmeg/A7/1kgp_reference --extract filtered2.prune.in --make-bed --out filtered4
```

## b. PCA computation

To assess the ethnic diversity in our cohort, we will use One-thousand
Genome Project (1KGP) data as a reference for our PCA analysis. These
dataset has genetic information of major continental populations:
Admixed American (AMR), European (EU), Asian (AS) and African (A).

``` bash
## Merge your pruned bfiles of the Mini_cohort and the 1KGP created on the previous step 
## Remember to create a new bfile (.fam, .bed and .bim files) that contains the merged data.
## IMPORTANT TIME CONSTRAINT: This step can take ~15 minutes, so make sure to check the server status before you run it!
#?# Type the command you used below: - 1pt
plink --bfile filtered3 --bmerge filtered4 --make-bed --out filtered5

#?# Perform a PCA analysis in plink on the merged set - 1 pt
plink --bfile filtered5 --pca 
```

## c. Visualization

``` r
## Copy the PCA .eigenvec file to your computer, together with the samples_info.txt located in /projects/bmeg/A7/

## Load the .eigenvec file onto R, change the column names to: FID, IID, PC1, PC2, PC3, ..., PC20
#?# Type the command you used below: - 1pt
#scp alzhang_bmeg22@orca1.bcgsc.ca:/projects/bmeg/A7/samples_info.txt /Users/alicezhang/400E_Assignments/Assignment_7/samples_info.txt
eigenvec_table <- read.table('plink.eigenvec')
colnames(eigenvec_table) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

## Load the samples_info.txt file onto R, change the column names to: FID, IID, SuperPopulation, Population
#?# Type the commands you used below: - 1pt
my_data <- read.delim("samples_info.txt")
colnames(my_data) <- c("IID", "FID", "SuperPopulation", "Population")

## Merge the .eigenvec and sample_info data.frames together using the IID column
## Tip: Look into the -merge- function!
#?# Type the command you used below: - 1pt
joined_df <- merge(my_data, eigenvec_table, by = "IID")

## Using ggplot create a scatterplot, using: 
## x-axis: PC1
## y-axis: PC2
## color: SuperPopulation - to use the Population information to color the samples and be able to appreciate population structure!
#?# Type the command you used below: 1pt
library(ggplot2)
ggplot(data = joined_df, aes(PC1, PC2, colour = SuperPopulation))+geom_point()
```

![](Assignment_7_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#?# Where do the cohort samples fall? Are they all clustered together? - 1 pt
#The cohort samples are clustered together on the lefthand side of the graph, near the SAS and AMR superpopulation samples.

#?# Which Population would you use as a reference for imputation?, Why? - 1 pt
#We can use the population EAS, that is east asian as a reference for imputation. 


#?# Do you think looking at the top two PCs is sufficient to tell what population is best? Why/why not? - 2 pt
#Looking at the top two PCs is sufficient to tell what population is best because PCAs are able to tell which components are able #to demonstrate the given data the best with a minimum number of components.
```

# Imputation

Imputation of genetic data is a very computationally intensive analysis,
that can take a long time. So we have performed it for you. Using the
chromosome 17 imputation information located in */projects/bmeg/A7/*
under the *Mini\_cohort\_chr17\_imputation\_results.info.gz* we will
calculate some post-imputation metrics.

``` r
## Load the Mini_cohort_chr17_imputation_results.info.gz file to your Rstudio environment 
chr_info <- read.table(file = "Mini_cohort_chr17_imputation_results.info.gz", header = TRUE)
## Use the information in the file to answer the following questions. Accompany each of the answers with the code you used to get to them and a brief explanation of your thought process behind.
#?# What is the percentage of imputed SNPs? 0.5 pt
print(length( which( (chr_info$Genotyped == "Imputed") ))/length(chr_info$Genotyped))
```

    ## [1] 0.9929754

``` r
# 99.3%

## The metric of imputation quality is Rsq, this is the estimated value of the squared correlation between imputed and true genotypes. Since true genotypes are not available, this calculation is based on the idea that poorly imputed genotype counts will shrink towards their expectations based on allele frequencies observed in the population (https://genome.sph.umich.edu/wiki/Minimac3_Info_File#Rsq).  An Rsq < 0.3 is often used to flag poorly imputed SNPs. 
#?# What is the percentage of poorly imputed SNPs?
print(length( which( (chr_info$Rsq < 0.3) ))/length(chr_info$Rsq))
```

    ## [1] 0.6341327

``` r
#63.4%

#?# Create a histogram to visualize the distribution of the MAF - 1 pt
hist(chr_info$MAF)
```

![](Assignment_7_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#?# Which MAF is most frequent? What does that mean? - 1 pt
# 0-0.02. This means that there is a low frequency of the second most common allele occuring in the sample population. This could mean the SNP's major allele is conserved. 

#?# What is the maximum MAF? Why is that? - 1 pt
max(chr_info$MAF)
```

    ## [1] 0.5

``` r
#The max MAF is 0.5. Because the MAF is the minor allele frequency, it cannot possibly have a higher frequency than the major allele frequency (the most common allele). Or else it would not be the minor allele frequency anymore!
```

# Polygenic Scores (PGS)

A GWAS for affinity for tapas (the Spanish appetizer) was performed and
199 SNPs were found significantly associated. The significant SNPs and
their assigned effect sizes are described in the
*Tapas\_enjoyability\_GWAS\_sumStats.txt* file. Thanks to the imputation
performed in our MiniCohort, we were able to obtain the dosages (double
risk alleles=2, one risk allele=1, no risk alleles=0) for each one of
the SNPs associated to the Tapas ‘enjoyability’, described in the
*MiniCohort\_Tapas\_SNPdosages.txt*.

PGS are calculated by multiplying the effect sizes of each SNP by the
dosage of an individual for those SNP and then adding together all the
effectSize x dosage. The formula is outlined below, where:

-   i: individual of which you are calculating the PGS

-   j: SNP that has been found to be associated to the trait (Tapas
    enjoyability in this case)

-   Beta: Effect size

-   dosage: number of risk alleles the *individual i* has of the *risk
    allele j*? (2,1 or 0)

![](PGS_formula.png)

``` r
## Load to your RStudio:
## 1.  -Tapas_enjoyability_GWAS_sumStats.txt-
## 2.  -MiniCohort_Tapas_SNPdosages.txt- 
## Both are located in the A7 directory on github.

## Using the base PRS formula outlined below, calculate the Tapas enjoyability PGS for the individuals in the Mini Cohort 
#?# Include your rationale and the documented code you used - 5pt
mc_snp <- read.table(file = "MiniCohort_Tapas_SNPdosages.txt", header = TRUE)
tapas <- read.table(file="Tapas_enjoyability_GWAS_sumStats.txt", header = TRUE)
#basically multiply mc_snp by tapas with, rows of tapas with the respective mc_snp columns
cut_snp <- mc_snp[,-1]
cut_snp <- cut_snp[,-1]
columns <- colnames(cut_snp)
pgs <- c(1:63) #random values to start with
for (i in 1:nrow(cut_snp)){
  val <- sum(cut_snp[i,]*tapas$Effect_Size)
  pgs[i] <- val
}
print(pgs)
```

    ##  [1]  24.84181  29.43423  12.72804  24.56438  12.78852   8.35869  14.65265
    ##  [8]  18.34110  30.43842  20.87267  31.25186   1.99843  23.20303   3.07489
    ## [15]   7.93791  25.22931   8.70430  31.31115  38.36467  17.80205  24.70952
    ## [22]  31.53787  27.73976  34.45540   5.59570  22.31392   7.58563 -10.40604
    ## [29]  -8.12046  -0.33614  14.24906  24.79131  17.26722  -0.39075  27.09533
    ## [36]  13.18553  24.15010  13.73062  23.79800   9.92227   4.47148  16.86456
    ## [43]  23.26737  32.91110  14.18713  28.50019  28.23938  22.84820  36.94977
    ## [50]   7.32103  18.08052  12.32953  32.65660   5.28892  29.95761  28.97731
    ## [57]  19.32970  16.18814   3.76833   9.64055   5.07950  13.91063   9.74386

``` r
#?# Use ggplot to plot the distribution of the Tapas PGS: - 2 pt
## Include the code and the graph in your analysis! 
## Tip: http://www.cookbook-r.com/Graphs/Plotting_distributions_(ggplot2)/
pgs <- data.frame(pgs)
ggplot(pgs, aes(x=pgs))+geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Assignment_7_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#?# What is the distribution of the tapas PGS? - 1pt
# The distribution looks random. 
```

## PGS accuracy

``` r
## The Tapas enjoyability was measured in a range of 0-1, with 0 being hating tapas and 1 being completely in love with tapas.
## This tapas likability is captured in the "Tapas_enjoyability" column of the -MiniCohort_Tapas_SNPdosages.txt- file. 
#?# Make a scatterplot with a linear regression line, where x is the Tapas-PGS and y is their actual Tapas enjoyability - 2 pt
## Tip: http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
pgs$enjoy <- mc_snp$Tapas_enjoyability
ggplot(pgs, aes(x= pgs, y=enjoy)) + geom_point()+geom_smooth()
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](Assignment_7_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#?# What is the correlation coefficient between the PGS and Tapas enjoyability? Is Spearman or Pearson correlation more appropriate here? Why? - 3 pt
cor(pgs$pgs, pgs$enjoy, method = "pearson")
```

    ## [1] 0.1409448

``` r
cor(pgs$pgs, pgs$enjoy, method = "spearman")
```

    ## [1] 0.171456

``` r
# The Pearson correlation can only be used to evaluate linear relationships between continuous relationships. Since the Tapas enjoyability is discrete, and since from the graph we see that the relationship between the PGS and Tapas enjoyability does not look linear, a Pearson correlation is not appropriate in this case. The Spearman correlation can be use to evaluate monotonic relationships. From the scatterplot the relationship doesn't look monotonic, since as x increases y does not clearly increase or decrease, and fluctuates. We calculated the correlation coefficients using both and found both to be insignificantly correlated (<0.2). 


#?# How predictive is the PGS for tapas preference? Include in your answer why do you think it is/isn't accurate and what could affect its predicitvity - 2pt 
# The PGS for tapas preference is not very predictive. It is not accurate because the pearson and spearman correlation are found to be insignificant, there also seems to be no correlation appear in the plot. There are a bunch of things that could affect the predictivity like personal preferences, experiences, environmental factors and external mix of genes. All these factors are not necessarily captured in genetics. 
#TODO why
```

# Authors and contributions

Following completion of your assignment, please fill out this section
with the authors and their contributions to the assignment. If you
worked alone, only the author (e.g. your name and student ID) should be
included.

Authors: Alice Zhang (52721552) and Abhishek Dhir (87603866)

Contributions: Alice and Abhishek worked together asynchronously on the
assignment together and then pushed the combined changes in the end,
after a synchronous working session.
