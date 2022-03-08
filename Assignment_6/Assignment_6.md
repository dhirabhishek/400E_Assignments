Assignment 6: ATAC-seq
================

-   [Assignment overview](#assignment-overview)
-   [Part 0: Getting ready](#part-0-getting-ready)
-   [Part 1: understanding the
    experiment](#part-1-understanding-the-experiment)
    -   [`#?#` *Make the above plot. Each point should represent one of
        the samples. - 1
        pt*](#-make-the-above-plot-each-point-should-represent-one-of-the-samples---1-pt)
    -   [`#?#` *Can we compare BRM014 to DMSO across all time points?
        Why/why not? - 1
        pt*](#-can-we-compare-brm014-to-dmso-across-all-time-points-whywhy-not---1-pt)
-   [Part 2: QC](#part-2-qc)
    -   [`#?#` Make a plot with read coverage on the y-axis (total
        number of reads) and the samples on the x-axis. - 3
        pt\*](#-make-a-plot-with-read-coverage-on-the-y-axis-total-number-of-reads-and-the-samples-on-the-x-axis---3-pt)
    -   [`#?#` *Which sample has the most coverage? - 0.5
        pt*](#-which-sample-has-the-most-coverage---05-pt)
    -   [`#?#` *Which sample has the least? - 0.5
        pt*](#-which-sample-has-the-least---05-pt)
    -   [`#?#` *What is the % difference between the max and min
        (relative to the min)? - 0.5
        pt*](#-what-is-the--difference-between-the-max-and-min-relative-to-the-min---05-pt)
    -   [`#?#` *Create a new data.frame containing only the BI_protac
        and control samples - 1
        pt*](#-create-a-new-dataframe-containing-only-the-bi_protac-and-control-samples---1-pt)
    -   [`#?#` *For this subset, calculate the counts per million reads
        (CPM) for each sample - 2
        pt*](#-for-this-subset-calculate-the-counts-per-million-reads-cpm-for-each-sample---2-pt)
    -   [`#?#` *Plot the kernel density estimate for CPM (x axis). 1
        curve per sample, different colours per curve. - 1
        pt*](#-plot-the-kernel-density-estimate-for-cpm-x-axis-1-curve-per-sample-different-colours-per-curve---1-pt)
    -   [`#?#` *Plot the kernel density estimate for log(CPM+1) (x
        axis), coloured as before - 1
        pt*](#-plot-the-kernel-density-estimate-for-logcpm1-x-axis-coloured-as-before---1-pt)
    -   [`#?#` *Why do you think log-transforming is usually performed
        when looking at genomics data? What about adding 1 before log
        transforming? - 2
        pt*](#-why-do-you-think-log-transforming-is-usually-performed-when-looking-at-genomics-data-what-about-adding-1-before-log-transforming---2-pt)
    -   [`#?#` *Some regions have very large CPMs. Inspect the peaks for
        which CPM>400. What do you notice about them? 3
        pt*](#-some-regions-have-very-large-cpms-inspect-the-peaks-for-which-cpm400-what-do-you-notice-about-them-3-pt)
    -   [`#?#` *Calculate the pairwise correlations between log(CPM+1)s
        for the samples and plot them as a heatmap (samples x samples) -
        3
        pt*](#-calculate-the-pairwise-correlations-between-logcpm1s-for-the-samples-and-plot-them-as-a-heatmap-samples-x-samples---3-pt)
    -   [`#?#` *What do you expect the correlations between replicates
        to look like? Is that what you see? - 2
        pt*](#-what-do-you-expect-the-correlations-between-replicates-to-look-like-is-that-what-you-see---2-pt)
    -   [`#?#` *Filter your data, retaining only regions where the
        average counts per sample is greater than 10, and also remove
        mitochondrial regions - 3
        pt*](#-filter-your-data-retaining-only-regions-where-the-average-counts-per-sample-is-greater-than-10-and-also-remove-mitochondrial-regions---3-pt)
    -   [`#?#` *How many peaks did you have before? How many do you have
        now? - 1
        pt*](#-how-many-peaks-did-you-have-before-how-many-do-you-have-now---1-pt)
-   [Part 3: Differential ATAC](#part-3-differential-atac)
    -   [`#?#` *Make a count matrix called `countMatrix` for the
        BI_protac and control samples, including only the peaks we
        retained above - 2
        pt*](#-make-a-count-matrix-called-countmatrix-for-the-bi_protac-and-control-samples-including-only-the-peaks-we-retained-above---2-pt)
    -   [`#?#` *Make an MA plot for allDEStatsPairedTreatControlvsProtac
        -2pt*](#-make-an-ma-plot-for-alldestatspairedtreatcontrolvsprotac--2pt)
    -   [`#?#` *Make an MA plot for allDEStatsPairedTime6vs24 - 1
        pt*](#-make-an-ma-plot-for-alldestatspairedtime6vs24---1-pt)
    -   [`#?#` *Perform the same differential peak analysis using loess
        regularization. - 1
        pt*](#-perform-the-same-differential-peak-analysis-using-loess-regularization---1-pt)
    -   [`#?#` *Make the same two MA plots as before, but this time
        using the loess normalized analysis - 1
        pt*](#-make-the-same-two-ma-plots-as-before-but-this-time-using-the-loess-normalized-analysis---1-pt)
    -   [`#?#` *What was the first normalization method? What changed in
        the MA plots? Which analysis do you think is more reliable and
        why? - 4
        pt*](#-what-was-the-first-normalization-method-what-changed-in-the-ma-plots-which-analysis-do-you-think-is-more-reliable-and-why---4-pt)
-   [Part 4: GC bias](#part-4-gc-bias)
    -   [`#?#` *Convert the region IDs to a GRanges object - 3
        pt*](#-convert-the-region-ids-to-a-granges-object---3-pt)
    -   [`#?#` *Extract the genomic DNA sequences for each peak using
        hg38 - 3
        pt*](#-extract-the-genomic-dna-sequences-for-each-peak-using-hg38---3-pt)
    -   [`#?#` *Create scatter plots (one per sample, e.g. using
        facet_wrap), including lines of best fit (GAM), where each plot
        shows GC content (x axis) vs CPM (y axis) for each peak (points)
        -2pt*](#-create-scatter-plots-one-per-sample-eg-using-facet_wrap-including-lines-of-best-fit-gam-where-each-plot-shows-gc-content-x-axis-vs-cpm-y-axis-for-each-peak-points--2pt)
    -   [`#?#` *Repeat the above, but this time showing only the lines
        of best fit and all on the same plot - 2
        pt*](#-repeat-the-above-but-this-time-showing-only-the-lines-of-best-fit-and-all-on-the-same-plot---2-pt)
    -   [`#?#` *Given this result, predict whether we will see a
        significant relationship between GC content and logFC in our
        differential peak analysis (loess-normalized). Justify your
        prediction. Predicting “wrong” will not be penalized, as long as
        your justification is correct. Don’t retroactively change your
        answer. - 2
        pt*](#-given-this-result-predict-whether-we-will-see-a-significant-relationship-between-gc-content-and-logfc-in-our-differential-peak-analysis-loess-normalized-justify-your-prediction-predicting-wrong-will-not-be-penalized-as-long-as-your-justification-is-correct-dont-retroactively-change-your-answer---2-pt)
    -   [`#?#` *Plot the relationship between GC and logFC for the
        loess-normalized ControlvsProtac analysis. Also include a line
        of best fit (blue) and y=0 (red) - 2
        pt*](#-plot-the-relationship-between-gc-and-logfc-for-the-loess-normalized-controlvsprotac-analysis-also-include-a-line-of-best-fit-blue-and-y0-red---2-pt)
    -   [`#?#` *Now plot the same thing for the NON loess-normalized
        ControlvsProtac analysis. - 1
        pt*](#-now-plot-the-same-thing-for-the-non-loess-normalized-controlvsprotac-analysis---1-pt)
    -   [`#?#` *Was your prediction correct? Do you think we should also
        account for GC normalization in our differential ATAC analysis?
        Why/why not? - 3
        pt*](#-was-your-prediction-correct-do-you-think-we-should-also-account-for-gc-normalization-in-our-differential-atac-analysis-whywhy-not---3-pt)
-   [Part 5: Differential analysis
    results](#part-5-differential-analysis-results)
    -   [`#?#` *Suppose we perform the analyses above, redoing the
        differential analysis once more with GC normalization, and also
        considering that we tested loess and the default normalization
        methods. Did we P-hack? Why or why not? - 2
        pt*](#-suppose-we-perform-the-analyses-above-redoing-the-differential-analysis-once-more-with-gc-normalization-and-also-considering-that-we-tested-loess-and-the-default-normalization-methods-did-we-p-hack-why-or-why-not---2-pt)
    -   [`#?#` *Now considering the two comparisons (6 vs 24 hours, and
        protac vs control). EdgeR performed a correction for MHT, but if
        we want to analyze the results from both comparisons, do we need
        to re-adjust to account for the fact that we tested two
        different hypothesis sets (time and treatment)? Why/not? - 2
        pt*](#-now-considering-the-two-comparisons-6-vs-24-hours-and-protac-vs-control-edger-performed-a-correction-for-mht-but-if-we-want-to-analyze-the-results-from-both-comparisons-do-we-need-to-re-adjust-to-account-for-the-fact-that-we-tested-two-different-hypothesis-sets-time-and-treatment-whynot---2-pt)
    -   [`#?#` *How many differential peaks did you find (FDR\<0.01). -
        1 pt*](#-how-many-differential-peaks-did-you-find-fdr001---1-pt)
    -   [`#?#` *Make a volcano plot of the
        allDEStatsPairedTreatControlvsProtac, with -log10(p-value) on
        the y axis and logFC on the x. Colour points that are
        significant at an FDR\<0.01. - 2
        pt*](#-make-a-volcano-plot-of-the-alldestatspairedtreatcontrolvsprotac-with--log10p-value-on-the-y-axis-and-logfc-on-the-x-colour-points-that-are-significant-at-an-fdr001---2-pt)
    -   [`#?#` *Plot the logCPM (x axis) by -log10(Pvalue) (y axis),
        again colouring by FDR\<0.01. - 2
        pt*](#-plot-the-logcpm-x-axis-by--log10pvalue-y-axis-again-colouring-by-fdr001---2-pt)
    -   [`#?#` *Do you think our initial filtering on peaks with at
        least 10 reads on average per sample was a good choice? Why or
        why
        not?*](#-do-you-think-our-initial-filtering-on-peaks-with-at-least-10-reads-on-average-per-sample-was-a-good-choice-why-or-why-not)
-   [Authors and contributions](#authors-and-contributions)

# Assignment overview

*Today we will be looking at a differential ATAC-seq dataset between
cells treated with an anti BAF protac and control (untreated) cells. The
cell type is HAP1, a cancer cell line with a near-haploid genome. We
will use this dataset to explore differential analysis. *

*The GEO entry is located here, where you can read more about the
experiments:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148175> *

*This is the paper: <https://www.nature.com/articles/s41588-021-00777-3>
*

*“Acute BAF perturbation causes immediate changes in chromatin
accessibility”*

# Part 0: Getting ready

``` r
#install any of these you might not have already

library(ggplot2)
library(edgeR)
```

    ## Loading required package: limma

``` r
library(reshape)
library(GenomicRanges)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:reshape':
    ## 
    ##     expand, rename

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

``` r
library(csaw)
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(Biostrings)
```

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
#download the data
atacSeqData = read.table(textConnection(readLines(gzcon(url("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148175/suppl/GSE148175_count_matrix_raw_atac_BRM014_ACBI1.csv.gz")))), 
                      sep=",", stringsAsFactors = FALSE, header = TRUE)
```

``` r
#create a sample metadata data.frame

samples = data.frame(ID = names(atacSeqData)[2:ncol(atacSeqData)])
samples$replicate = gsub("(R[12])_([0-9]+[minh]+)_(.*)$","\\1",samples$ID)
samples$timeName = gsub("(R[12])_([0-9]+[minh]+)_(.*)$","\\2",samples$ID)
samples$treatment = gsub("(R[12])_([0-9]+[minh]+)_(.*)$","\\3",samples$ID)
samples$treatment[samples$treatment=="N"]="BRM014"
samples$time= as.numeric(gsub("[a-z]*","",samples$timeName))
samples$time[grepl("min",samples$timeName)]=samples$time[grepl("min",samples$timeName)]/60
```

# Part 1: understanding the experiment

*Now using `samples` make a plot showing the experimental design, with
time on the x axis, treatment on the y axis, and one plot on the left
and one on the right for the two replicates (e.g. using `facet_grid`).*

### `#?#` *Make the above plot. Each point should represent one of the samples. - 1 pt*

``` r
#here, if the point is there, it means such a sample exists, if absent it means that there is no such sample
ggplot(data =samples, aes(x=time, y = treatment))+
  geom_point()+
  facet_wrap(~replicate)
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

*In this study, one of the things they were comparing was BRM014 to
DMSO. The drug BRM014 is dissolved in DMSO, so DMSO alone is the
appropriate control to gauge the effect of BRM014.*

### `#?#` *Can we compare BRM014 to DMSO across all time points? Why/why not? - 1 pt*

###No, some sample points only exist for some time points in BRM014, and
not DMSO.

# Part 2: QC

*With most genomics data, it is important both that samples have
sufficient coverage, and that the samples have similar coverage. Either
case can lead to underpowered analysis, or misleading results. Calcualte
the read coverage for each sample. *

### `#?#` Make a plot with read coverage on the y-axis (total number of reads) and the samples on the x-axis. - 3 pt\*

``` r
# there are many ways you could do this; one of which is using the melt/cast functions from reshape

meltedSeqData<-melt(atacSeqData)
```

    ## Using region as id variables

``` r
meltedSeqData$treatment = gsub("(R[12])_([0-9]+[minh]+)_(.*)$","\\3",meltedSeqData$variable)
readCoverage <- aggregate(meltedSeqData$value, list(meltedSeqData$variable), sum)
colnames(readCoverage) <- c("treatment", "reads")
ggplot(data =readCoverage, aes(x=treatment, y = reads))+
  geom_point()
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### `#?#` *Which sample has the most coverage? - 0.5 pt*

R1_24h_DMSO has the most coverage

``` r
readCoverage[which.max(readCoverage$reads),]
```

    ##     treatment   reads
    ## 1 R1_24h_DMSO 4686253

### `#?#` *Which sample has the least? - 0.5 pt*

Control

``` r
readCoverage[which.min(readCoverage$reads),]
```

    ##        treatment   reads
    ## 12 R1_6h_control 2341332

### `#?#` *What is the % difference between the max and min (relative to the min)? - 0.5 pt*

``` r
(max(readCoverage$reads)-min(readCoverage$reads))/min(readCoverage$reads)*100
```

    ## [1] 100.1533

100.1533%

*In cases where samples have vastly different coverage, you can
potentially down-sample the higher-coverage samples. Sometimes, throwing
out the data in this way can also introduce new problems, so we’re going
to stick with the data we have.*

*For this assignment, we will look only at BI_protac vs control data. *

### `#?#` *Create a new data.frame containing only the BI_protac and control samples - 1 pt*

``` r
downsampledSeqData<-meltedSeqData[meltedSeqData$treatment == "BI_protac" | meltedSeqData$treatment == "control", ]
```

### `#?#` *For this subset, calculate the counts per million reads (CPM) for each sample - 2 pt*

``` r
library(edgeR)
protac_reads <- downsampledSeqData[downsampledSeqData$treatment == "BI_protac", ]
control_reads <- downsampledSeqData[downsampledSeqData$treatment == "control", ]

readCoverageProtac <- aggregate(protac_reads$value, list(protac_reads$variable), sum)
readCoverageControl <- aggregate(control_reads$value, list(control_reads$variable), sum)

protac_reads_R1_24 <- downsampledSeqData[downsampledSeqData$variable == "R1_24h_BI_protac", ]
control_reads_R1_24 <- downsampledSeqData[downsampledSeqData$variable == "R1_24h_control", ]
protac_reads_R1_6 <- downsampledSeqData[downsampledSeqData$variable == "R1_6h_BI_protac", ]
control_reads_R1_6 <- downsampledSeqData[downsampledSeqData$variable == "R1_6h_control", ]
protac_reads_R2_24 <- downsampledSeqData[downsampledSeqData$variable == "R2_24h_BI_protac", ]
control_reads_R2_24 <- downsampledSeqData[downsampledSeqData$variable == "R2_24h_control", ]
protac_reads_R2_6 <- downsampledSeqData[downsampledSeqData$variable == "R2_6h_BI_protac", ]
control_reads_R2_6 <- downsampledSeqData[downsampledSeqData$variable == "R2_6h_control", ]

CPM_protac_R1_24 <- cpm(protac_reads_R1_24$value)
CPM_control_R1_24 <- cpm(control_reads_R1_24$value)
CPM_protac_R1_6 <- cpm(protac_reads_R1_6$value)
CPM_control_R1_6 <- cpm(control_reads_R1_6$value)
CPM_protac_R2_24 <- cpm(protac_reads_R2_24$value)
CPM_control_R2_24 <- cpm(control_reads_R2_24$value)
CPM_protac_R2_6 <- cpm(protac_reads_R2_6$value)
CPM_control_R2_6 <- cpm(control_reads_R2_6$value)
print(max(CPM_control_R1_24))
```

    ## [1] 41162.2

``` r
#readCoverageProtac$CPM_protac <- readCoverageProtac$x/1000000
#readCoverageControl$CPM_control <- readCoverageControl$x/1000000
```

### `#?#` *Plot the kernel density estimate for CPM (x axis). 1 curve per sample, different colours per curve. - 1 pt*

``` r
protac_densityR1_24 <- density(CPM_protac_R1_24)
control_densityR1_24 <- density(CPM_control_R1_24)
protac_densityR1_6 <- density(CPM_protac_R1_6)
control_densityR1_6 <- density(CPM_control_R1_6)
protac_densityR2_24 <- density(CPM_protac_R2_24)
control_densityR2_24 <- density(CPM_control_R2_24)
protac_densityR2_6 <- density(CPM_protac_R2_6)
control_densityR2_6 <- density(CPM_control_R2_6)

plot(protac_densityR1_24, type = "l", col = "red")
lines(control_densityR1_24,col="blue")
lines(protac_densityR1_6, col = "green")
lines(control_densityR1_6, col = "orange")
lines(protac_densityR2_24, col = "pink")
lines(control_densityR2_24, col = "yellow")
lines(protac_densityR2_6, col = "purple")
lines(control_densityR2_6, col = "black")
legend( x = "topright", legend=c("control_R1_24", "BI_protac_R1_24", "control_R1_6", "control_R1_6", "control_R2_24", "BI_protac_R2_24", "control_R2_6", "control_R2_6"), 
        col=c("blue", "red", "green", "orange", "pink", "yellow", "purple", "black"), pch=c("-","-", "-","-", "-","-", "-","-"))
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### `#?#` *Plot the kernel density estimate for log(CPM+1) (x axis), coloured as before - 1 pt*

``` r
plot(density(log(CPM_protac_R1_24+1)), type = "l", col = "red")
lines(density(log(CPM_control_R1_24+1)),col="blue")
lines(density(log(CPM_protac_R1_6+1)), col = "green")
lines(density(log(CPM_control_R1_6+1)), col = "orange")
lines(density(log(CPM_protac_R2_6+1)), col = "pink")
lines(density(log(CPM_control_R2_6+1)), col = "yellow")
lines(density(log(CPM_protac_R2_6+1)), col = "purple")
lines(density(log(CPM_control_R2_6+1)), col = "black")
legend( x = "topright", legend=c("control_R1_24", "BI_protac_R1_24", "control_R1_6", "control_R1_6", "control_R2_24", "BI_protac_R2_24", "control_R2_6", "control_R2_6"), 
        col=c("blue", "red", "green", "orange", "pink", "yellow", "purple", "black"), pch=c("-","-", "-","-", "-","-", "-","-"))
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### `#?#` *Why do you think log-transforming is usually performed when looking at genomics data? What about adding 1 before log transforming? - 2 pt*

It is likely performed when looking at genomics because the unprocessed
data is very skewed (many 0s). 1 is added to avoid taking the log
transform of 0, which is common in the data.

### `#?#` *Some regions have very large CPMs. Inspect the peaks for which CPM>400. What do you notice about them? 3 pt*

``` r
largeCPM <- downsampledSeqData[cpm(downsampledSeqData$value) > 400, ]
```

They are all from times 6 hours or 24 hours. *Normally, we would remove
some of these regions before continuing (and would redo the above
steps). Since this is an assignment, we will continue with the data
as-is.*

*Often a good first step is to see if the data look good. One way to do
this is by seeing whether or not the signals in each sample correlate
with each other in ways you expect.*

### `#?#` *Calculate the pairwise correlations between log(CPM+1)s for the samples and plot them as a heatmap (samples x samples) - 3 pt*

``` r
df <- data.frame(CPM_control_R1_24, CPM_protac_R1_24, CPM_control_R1_6, CPM_protac_R1_6, CPM_control_R2_24, CPM_protac_R2_24, CPM_control_R2_6, CPM_protac_R2_6)

corr_mat <- cor(df)
melted_cormat <- melt(corr_mat)
```

    ## Warning in type.convert.default(X[[i]], ...): 'as.is' should be specified by the
    ## caller; using TRUE

    ## Warning in type.convert.default(X[[i]], ...): 'as.is' should be specified by the
    ## caller; using TRUE

``` r
ggplot(data = melted_cormat, aes(x=X1, y=X2, fill=value)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + geom_tile()
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

### `#?#` *What do you expect the correlations between replicates to look like? Is that what you see? - 2 pt*

We expected the correlations between replicates to be pretty high, which
is the case for most of the samples.

*It is common to exclude some regions from analysis. For instance, we
won’t be able to robustly identify those that are differential but have
low coverage even if they are truly differential, so there is no point
testing these. We will also remove mitochondrial regions, a common
contaminant of ATAC-seq data.*

### `#?#` *Filter your data, retaining only regions where the average counts per sample is greater than 10, and also remove mitochondrial regions - 3 pt*

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
    ## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
    ## ✓ readr   2.1.2     ✓ forcats 0.5.1
    ## ✓ purrr   0.3.4

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::collapse()        masks Biostrings::collapse(), IRanges::collapse()
    ## x dplyr::combine()         masks Biobase::combine(), BiocGenerics::combine()
    ## x purrr::compact()         masks XVector::compact()
    ## x dplyr::count()           masks matrixStats::count()
    ## x dplyr::desc()            masks IRanges::desc()
    ## x tidyr::expand()          masks S4Vectors::expand(), reshape::expand()
    ## x dplyr::filter()          masks stats::filter()
    ## x dplyr::first()           masks S4Vectors::first()
    ## x dplyr::lag()             masks stats::lag()
    ## x BiocGenerics::Position() masks ggplot2::Position(), base::Position()
    ## x purrr::reduce()          masks GenomicRanges::reduce(), IRanges::reduce()
    ## x dplyr::rename()          masks S4Vectors::rename(), reshape::rename()
    ## x dplyr::slice()           masks XVector::slice(), IRanges::slice()

``` r
filteredData <- downsampledSeqData[downsampledSeqData$value > 10,]
filteredData <- filteredData %>% filter(!str_detect(region, 'chrM'))
```

### `#?#` *How many peaks did you have before? How many do you have now? - 1 pt*

``` r
#assuming peak is a number in value, i might be wrong TODO: Alice please check 
print(sum(downsampledSeqData$value))
```

    ## [1] 26219406

``` r
sum(filteredData$value)
```

    ## [1] 22415384

Before we had 26219406 peaks but now we have 22415384 peaks.

# Part 3: Differential ATAC

*We want to know what regions are differentially accessible between
BI_protac and the control.*

*Today, we’re going to use edgeR, which is designed for RNA-seq, but
works well on ATAC-seq as well. The user guide is here:*
<https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>

### `#?#` *Make a count matrix called `countMatrix` for the BI_protac and control samples, including only the peaks we retained above - 2 pt*

``` r
library(reshape2)
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

    ## The following objects are masked from 'package:reshape':
    ## 
    ##     colsplit, melt, recast

``` r
library(dplyr)
countMatrix <- dcast(filteredData, region~variable)
rownames(countMatrix)<-(countMatrix$region)
countMatrix <- subset(countMatrix, select = -region)
countMatrix <- as.matrix(countMatrix)
#countMatrix[is.na(countMatrix)] <- 0
countMatrix <- na.omit(countMatrix)
```

*EdgeR is exceptionally versatile, with many different options for
analysis. Today, you’re going to use the GLM-quasi-likelihood approach
to calculate differential accessibility. We are providing the first
example analysis below, which you can modify in subsequent steps. You
will need to understand what the steps do, so read the appropriate
documentation. *

``` r
curSamples = samples[match(colnames(countMatrix), samples$ID),]; #some comments are added by Abhishek
y = DGEList(counts=countMatrix, group=curSamples$treatment)
y = calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes.
designPaired = model.matrix(~curSamples$treatment + curSamples$timeName)   #returning a dense or sparse object of class
# we are using timeName here to make sure that time is treated as a categorical variable. Had we more time points it might make sense to treat time as a value.
y = estimateDisp(y, designPaired) #Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fitPaired = glmQLFit(y, designPaired)
qlfPairedTime6vs24 = glmQLFTest(fitPaired, coef=3) 
qlfPairedTreatControlvsProtac = glmQLFTest(fitPaired, coef=2)

allDEStatsPairedTreatControlvsProtac = as.data.frame(topTags(qlfPairedTreatControlvsProtac,n=nrow(countMatrix)))
allDEStatsPairedTreatControlvsProtac$region=row.names(allDEStatsPairedTreatControlvsProtac)

allDEStatsPairedTime6vs24 = as.data.frame(topTags(qlfPairedTime6vs24,n=nrow(countMatrix)))
allDEStatsPairedTime6vs24$region=row.names(allDEStatsPairedTime6vs24)
```

*While the differential analysis has been done in this case, before we
look at the results, we are going to check if the data appear to be
normalized correctly. Also include a loess line of best fit, and the
line y=0.*

### `#?#` *Make an MA plot for allDEStatsPairedTreatControlvsProtac -2pt*

``` r
maPlot(logFC = allDEStatsPairedTreatControlvsProtac$logFC, logAbundance = allDEStatsPairedTreatControlvsProtac$logCPM)

ll1 <- loess(logFC ~ logCPM, data = allDEStatsPairedTreatControlvsProtac)
j <- order(allDEStatsPairedTreatControlvsProtac$logCPM)
lines(allDEStatsPairedTreatControlvsProtac$logCPM[j], ll1$fitted[j], col="red")
abline(h=0, col = "green")
legend( x = "topright", legend=c("y=0", "loess line of best fit"), 
        col=c("green", "red"), pch=c("-","-"))
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### `#?#` *Make an MA plot for allDEStatsPairedTime6vs24 - 1 pt*

``` r
maPlot(logFC = allDEStatsPairedTime6vs24$logFC, logAbundance = allDEStatsPairedTime6vs24$logCPM)

ll1 <- loess(logFC ~ logCPM, data = allDEStatsPairedTime6vs24)
j <- order(allDEStatsPairedTime6vs24$logCPM)
lines(allDEStatsPairedTime6vs24$logCPM[j], ll1$fitted[j], col="red")
abline(h=0, col = "green")
legend( x = "topright", legend=c("y=0", "loess line of best fit"), 
        col=c("green", "red"), pch=c("-","-"))
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

*Now we’re going to test loess normalization instead.*

### `#?#` *Perform the same differential peak analysis using loess regularization. - 1 pt*

``` r
#Note: the Bioconductor package csaw implements loess regularization in a way that is compatible with edgeR
## Tip: use the csaw library to implement the loess regularization
library(csaw)
curSamples = samples[match(colnames(countMatrix), samples$ID),]; #some comments are added by Abhishek
y = DGEList(counts=countMatrix, group=curSamples$treatment)
y = normOffsets(y) #Calculate normalization factors to scale the raw library sizes.
designPaired = model.matrix(~curSamples$treatment + curSamples$timeName)   #returning a dense or sparse object of class
# we are using timeName here to make sure that time is treated as a categorical variable. Had we more time points it might make sense to treat time as a value.
y = estimateDisp(y, designPaired) #Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fitPaired = glmQLFit(y, designPaired)
qlfPairedTime6vs24 = glmQLFTest(fitPaired, coef=3) 
qlfPairedTreatControlvsProtac = glmQLFTest(fitPaired, coef=2)

allDEStatsPairedTreatControlvsProtac = as.data.frame(topTags(qlfPairedTreatControlvsProtac,n=nrow(countMatrix)))
allDEStatsPairedTreatControlvsProtac$region=row.names(allDEStatsPairedTreatControlvsProtac)

allDEStatsPairedTime6vs24 = as.data.frame(topTags(qlfPairedTime6vs24,n=nrow(countMatrix)))
allDEStatsPairedTime6vs24$region=row.names(allDEStatsPairedTime6vs24)
```

### `#?#` *Make the same two MA plots as before, but this time using the loess normalized analysis - 1 pt*

``` r
maPlot(logFC = allDEStatsPairedTreatControlvsProtac$logFC, logAbundance = allDEStatsPairedTreatControlvsProtac$logCPM)

ll1 <- loess(logFC ~ logCPM, data = allDEStatsPairedTreatControlvsProtac)
j <- order(allDEStatsPairedTreatControlvsProtac$logCPM)
lines(allDEStatsPairedTreatControlvsProtac$logCPM[j], ll1$fitted[j], col="red")
abline(h=0, col = "green")
legend( x = "topright", legend=c("y=0", "loess line of best fit"), 
        col=c("green", "red"), pch=c("-","-"))
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
maPlot(logFC = allDEStatsPairedTime6vs24$logFC, logAbundance = allDEStatsPairedTime6vs24$logCPM)

ll1 <- loess(logFC ~ logCPM, data = allDEStatsPairedTime6vs24)
j <- order(allDEStatsPairedTime6vs24$logCPM)
lines(allDEStatsPairedTime6vs24$logCPM[j], ll1$fitted[j], col="red")
abline(h=0, col = "green")
legend( x = "topright", legend=c("y=0", "loess line of best fit"), 
        col=c("green", "red"), pch=c("-","-"))
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

### `#?#` *What was the first normalization method? What changed in the MA plots? Which analysis do you think is more reliable and why? - 4 pt*

The first normalization used the calcNormFactors method to calculate
normalization factors, which according to the documentation normalizes
by the TMM method, which trims the means of the M-values, which are
weighted using the inverse variances of the variables. Since we didn’t
set refColumn, the column with the CPM upper quartile closest to the
mean upper quartile was used as a reference library to scale the
observed library. Using the loess regulariazation, the plot generated
look similar, but in the upper and lower left hand corners the
previously linear lines have been evened out, and the data points are
more scattered.

# Part 4: GC bias

*Next, we will look at potential GC bias in the data. We will again use
bioconductor *

### `#?#` *Convert the region IDs to a GRanges object - 3 pt*

``` r
#note that the names of your peaks are of the format <chr>:<startPos>-<endPos>
## Tip: lookinto the GenomicRanges documentation 
library(GenomicRanges)
gr <- as(atacSeqData$region, "GRanges")
```

### `#?#` *Extract the genomic DNA sequences for each peak using hg38 - 3 pt*

*See for relevant documentation:
<https://bioconductor.org/packages/release/workflows/vignettes/sequencing/inst/doc/sequencing.html>
*

``` r
## Tip: Use the Biostring library 
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")# if not installed
library(BSgenome.Hsapiens.UCSC.hg38)
```

    ## Loading required package: BSgenome

    ## Loading required package: rtracklayer

``` r
genome <- BSgenome.Hsapiens.UCSC.hg38
j<-getSeq(genome, gr)
```

*Now we will see if there’s any relationship between peak CPM and GC
content for each of the samples.*

### `#?#` *Create scatter plots (one per sample, e.g. using facet_wrap), including lines of best fit (GAM), where each plot shows GC content (x axis) vs CPM (y axis) for each peak (points) -2pt*

``` r
#please limit the y axis to between 0 and 50
alf <- as.data.frame(alphabetFrequency(j))
#alf$Name <- gr
alf$Width <- width(j)
alf$gc <- (alf$C+alf$G) / (alf$C+alf$G+alf$A+alf$T)
alf$region <- atacSeqData$region

df2 <- data.frame(alf$gc, CPM_control_R1_24, CPM_protac_R1_24, CPM_control_R1_6, CPM_protac_R1_6, CPM_control_R2_24, CPM_protac_R2_24, CPM_control_R2_6, CPM_protac_R2_6)
meltedDF2 <- melt(df2, id.vars = "alf.gc")

ggplot(data = meltedDF2, aes(alf.gc, value)) +
  geom_line()+
  geom_point() + 
  labs(title = "CPM vs GC content",
       y = "CPM", x = "GC Content") +
      facet_wrap(~variable)
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

### `#?#` *Repeat the above, but this time showing only the lines of best fit and all on the same plot - 2 pt*

``` r
ggplot(data = meltedDF2, aes(alf.gc, value)) +
  geom_line()+
  labs(title = "CPM vs GC content",
       y = "CPM", x = "GC Content")
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### `#?#` *Given this result, predict whether we will see a significant relationship between GC content and logFC in our differential peak analysis (loess-normalized). Justify your prediction. Predicting “wrong” will not be penalized, as long as your justification is correct. Don’t retroactively change your answer. - 2 pt*

We don’t think there will be a significant relationship between GC
content and the logFC in our differential peak analysis. The CPM and GC
don’t seem to have an obvious relation from the graphs above.

### `#?#` *Plot the relationship between GC and logFC for the loess-normalized ControlvsProtac analysis. Also include a line of best fit (blue) and y=0 (red) - 2 pt*

### `#?#` *Now plot the same thing for the NON loess-normalized ControlvsProtac analysis. - 1 pt*

### `#?#` *Was your prediction correct? Do you think we should also account for GC normalization in our differential ATAC analysis? Why/why not? - 3 pt*

*We will leave GC normalization as an optional exercise, and will not
actually do it here.*

# Part 5: Differential analysis results

### `#?#` *Suppose we perform the analyses above, redoing the differential analysis once more with GC normalization, and also considering that we tested loess and the default normalization methods. Did we P-hack? Why or why not? - 2 pt*

*Going forward, we will only use the initial analysis (**not loess
normalized**)*

### `#?#` *Now considering the two comparisons (6 vs 24 hours, and protac vs control). EdgeR performed a correction for MHT, but if we want to analyze the results from both comparisons, do we need to re-adjust to account for the fact that we tested two different hypothesis sets (time and treatment)? Why/not? - 2 pt*

### `#?#` *How many differential peaks did you find (FDR\<0.01). - 1 pt*

``` r
pairedControlvsProtacFDR <- allDEStatsPairedTreatControlvsProtac[allDEStatsPairedTreatControlvsProtac$FDR < 0.01,]
```

### `#?#` *Make a volcano plot of the allDEStatsPairedTreatControlvsProtac, with -log10(p-value) on the y axis and logFC on the x. Colour points that are significant at an FDR\<0.01. - 2 pt*

``` r
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
```

    ## Loading required package: ggrepel

    ## Registered S3 methods overwritten by 'ggalt':
    ##   method                  from   
    ##   grid.draw.absoluteGrob  ggplot2
    ##   grobHeight.absoluteGrob ggplot2
    ##   grobWidth.absoluteGrob  ggplot2
    ##   grobX.absoluteGrob      ggplot2
    ##   grobY.absoluteGrob      ggplot2

``` r
allDEStatsPairedTreatControlvsProtac$logPvalue <- -log10(allDEStatsPairedTreatControlvsProtac$PValue)
keyvals.colour <- ifelse(rownames(allDEStatsPairedTreatControlvsProtac) %in% rownames(pairedControlvsProtacFDR), 'red','transparent')
names(keyvals.colour)[keyvals.colour == 'red'] <- 'significant'
EnhancedVolcano(toptable = allDEStatsPairedTreatControlvsProtac, x = 'logFC', y = 'PValue', lab= rownames(allDEStatsPairedTreatControlvsProtac), selectLab = rownames(pairedControlvsProtacFDR),colCustom = keyvals.colour  )
```

![](Assignment_6_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

### `#?#` *Plot the logCPM (x axis) by -log10(Pvalue) (y axis), again colouring by FDR\<0.01. - 2 pt*

``` r
EnhancedVolcano(toptable = allDEStatsPairedTreatControlvsProtac, x = 'logCPM', y = 'PValue', lab= rownames(allDEStatsPairedTreatControlvsProtac), selectLab = rownames(pairedControlvsProtacFDR),colCustom = keyvals.colour  )
```

    ## Warning: Removed 2 rows containing missing values (geom_vline).

![](Assignment_6_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

### `#?#` *Do you think our initial filtering on peaks with at least 10 reads on average per sample was a good choice? Why or why not?*

*At this point there are many other follow ups you can and would do for
a real differential analysis, but we leave these as optional exercises.
For example:* 1. Confirming that the differential peaks look correct
(e.g. CPM heatmap) 2. Confirming that peaks look differential on the
genome browser 3. Looking for motif enrichment 4. Performing a GREAT
analysis, including functional enrichment and assigning peaks to genes

*Knit your assignment as a github_document and submit the resulting .md
and this .Rmd to your github, and complete the assignment submission on
Canvas. Make sure to include the graphs with your submission. *

# Authors and contributions

Following completion of your assignment, please fill out this section
with the authors and their contributions to the assignment. If you
worked alone, only the author (e.g. your name and student ID) should be
included.

Authors: Name1 (studentID1) and Name2 (studentID2)

Contributions: (example) N1 and N2 worked together on the same computer
to complete the assignment. N1 typed for the first half and N2 typed for
the second half.
