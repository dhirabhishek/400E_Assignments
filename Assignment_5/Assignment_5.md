Assignment 5: ChIP analysis
================

-   [Assignment Overview](#assignment-overview)
    -   [0. Getting Started](#getting-started)
    -   [1. ChIP signal tracks](#chip-signal-tracks)
    -   [2. Narrow vs Broad peaks](#narrow-vs-broad-peaks)
    -   [3. Peak calling](#peak-calling)
        -   [a. Peak calling with macs2](#a.-peak-calling-with-macs2)
        -   [b. Understanding the peaks](#b.-understanding-the-peaks)
        -   [c. Peak enrichments](#c.-peak-enrichments)
-   [Assignment submission](#assignment-submission)
-   [Authors and contributions](#authors-and-contributions)

# Assignment Overview

By now you must have become vaguely familiar with ChIP-seq data but
might be a little confused about what to do after the alignment. Well,
this assignment’s aim is to walk you through a ChIP-seq pipeline
post-alignment. We will be analyzing 3 different histone modification
marks (H3K27me3, H3K4me3 and H3K27ac). In order to identify the
enrichments for each epigenetic mark, we also need to use the *input*
which represents the DNA content of the sheared chromatin sample prior
to immunoprecipitation. All the files can be found under the following
path: **/projects/bmeg/A5/** .

-   H3K27me3 (H3K27me3\_chr3\_subset.bam)

-   H3K4me3 (H3K4me3\_chr3\_subset.bam)

-   H3K27ac (H3K27ac\_chr3\_subset.bam)

-   input (input\_chr3\_subset.bam)

A couple of things to remember:

-   When uploading your completed assignment to your GitHub directory,
    remember to specify a **github\_document** instead of the default
    *html\_document* on your .Rmd file.

-   Double check that all the files have been uploaded to your
    repository and that you are able to visualize the html-like version
    on GitHub.

-   Be mindful of your space on the server! Delete ALL unnecessary files
    and keep a clean environment.

## 0. Getting Started

We will be using a couple of new tools this time. Before we move on to
the practical part, make sure you have them all installed.

-   Integrative Genomics Viewer (IGV): Interactive tool to visualize
    different data types of genetic information (e.g. bam, bed files).
    You will install this tool to your **local computer**. To visualize
    where the reads of our ChIP analysis mapped in the genome. To
    install it, follow the instructions on this website:
    *<https://software.broadinstitute.org/software/igv/home>*

-   Deeptools
    (<https://deeptools.readthedocs.io/en/develop/index.html>): Software
    to analyze high-throughput data that allows to create easy to
    visualize figures. This will be installed on the server as part of
    your conda environment.

-   macs2
    (<https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md>):
    Tool to capture enrichments of ChIP-seq analysis. This will be
    installed on the server as part of your conda environment.

``` bash
#?# Add macs2 and deeptools to your environment created on A1 - 1 pt
conda install -c bioconda macs2
conda install -c bioconda deeptools
## Install IGV to your local computer after downloading it from: https://software.broadinstitute.org/software/igv/home
```

## 1. ChIP signal tracks

ChIP-seq experiments require a control against which to compare the ChIP
signnal. Often this is the “input DNA”, sheared chromatin that has not
been immmunoprecipitated. For this assignment we are using an input and
three different epigenetic marks. These histone modifications mark
states of active (H3K27ac, H3K4me3) or inactive (H3K27me3) gene
transcription and have different coverage of the genomic region where
they are located. To better visualize the differences, we will create
bigWig files from previously aligned, filtered and indexed bam files.
BigWig files are indexed, compressed files that can be used to visualize
signals across the genome. Here, we will be using it to graph the
coverage across the genome.

``` bash
## Use the bamCoverage command (included in deepTools) to convert the bam files outlined below (located in the /projects/bmeg/A5/ directory) into bigwig files that represent the read coverage (e.g. density of reads across the genome). Each file has already been indexed using sambamba and hence, has it's bam index (.bai) file associated that you will need to run bamCoverage. 
## Tip: Remember to add the .bw extension to your specified output! 

#?# Type the commands you use for converting each of the 4 bam files to bigwig files - 2 pts
## input_chr3_subset.bam
bamCoverage -b /projects/bmeg/A5/input_chr3_subset.bam -o input_chr3.bw

## H3K4me3_chr3_subset.bam
bamCoverage -b /projects/bmeg/A5/H3K4me3_chr3_subset.bam -o H3K4me3_chr3.bw

## H3K27me3_chr3_subset.bam
bamCoverage -b /projects/bmeg/A5/H3K27me3_chr3_subset.bam -o H3K27me3_chr3.bw

## H3K27ac_chr3_subset.bam
bamCoverage -b /projects/bmeg/A5/H3K27ac_chr3_subset.bam -o H3K27ac_chr3.bw 


### Follow this steps: 

## 1. Copy the bigwig files from the previous steps to your computer 

## 2.Load all the files bigwig signal track files onto IGV on your local computer, select the "autoscale" option for each file on their individual tracks
## Tip: use the "File" tab to "load from file" option to choose the files from your computer directories

## 3. Change the visualization of the files to see the following position: ** chr3:93,470,124-93,471,058 **
#?# Include a screenshot of your IGV session right after this code chunk (using the Rmarkdown syntax)  - 2 pt

#?# Explore this region by zooming in and out and navigating around. What do you see? Is there something similar across all the different files that stands out on this region? Is there anything peculiar about the DNA sequence at this locus?- 3 pt
# All the files look the same to our eyes on IGV. It looks like they have very similar peaks in this region. But when autoscale is turned on, the peaks look identical to a naked eye, and the shape of it looks like a bell curve.
## Tip: sometimes track signal will be truncated to a pre-set maximum. If you right-click the track label (left) and select "autoscale", it will automatically adust the scales of the track to the range of your data. Try playing around with the track settings to see their effect.

## 4. This file (/projects/bmeg/A5/hg38_blacklisted_regions.bed) contains the hg38 blacklisted regions. Load it into IGV along with all the other files. 

## 5. Look at the following region again (chr3:93,470,124-93,471,058). 
#?# What can you say now about the similarities between the files at this region? In your own words, explain what a blacklisted region is and if you think this region should be excluded a ChIP-seq analysis. - 1.5 pt
# A blacklisted region is where the genome assembly results in erroneous signal. These signals need to be filtered out befor
#ChIP-seq to ensure that there are not biases in the results.
```

**Add screenshot of your IGV session here: **. (Use Rmarkdown syntax)
![](igv_snapshot_autoscale.png)

## 2. Narrow vs Broad peaks

While exploring the bigwig files of the epigenetic marks on IGV, you
probably noticed that they can look very different from each other and
some of them reassemble the input more closely than others. Different
epigenetic marks can have very different signatures based on their
distribution patterns across the genome. When calling peaks, we often
lump these into two different categories of reads: broad and narrow.
Active transcription marks (H3K4me3 and H3K27ac) tend to form a sharper
coverage peaks at transcription start sites (H3K27ac is also at
enhancers), while repression marks cover a broader area (H3K27me3).
Here, we’re going to inspect their distributions relative to genes.

``` bash

## Here, we have created three bigWig track files, one for each epigenetic mark, which show the read coverage normalized using the input. They are found here: /projects/bmeg/A5/
# H3K4me3_norm.bw
# H3K27me3_norm.bw
# H3K27ac_norm.bw


## The deepTools ** computeMatrix reference-point ** command calculates scores to represent the reads mapping to specified regions of the genome across different files. 
## Use computeMatrix to compute a matrix for the signal tracks for each histone modification outlined above (which we will use to create a plot in the following step), with the following criteria: 

## - We will use the regions in reference_genes.bed located under the /projects/bmeg/A5/ directory as the basis for the plot.
## - Include the surrounding 1kb
## - Use all 3 input-normalized bigWig files (H3K4me3, H3K27ac, H3K27me3) as signal tracks
#?# Write the command you used to run it below: - 1.5 pt

#computeMatrix reference-point -S /projects/bmeg/A5/H3K4me3_norm.bw /projects/bmeg/A5/H3K27me3_norm.bw /projects/bmeg/A5/H3K27ac_norm.bw -R /projects/bmeg/A5/reference_genes.bed -o reference.tab -a 1000 -b 1000

## Now that the scores matrix has been computed, we can use it to create a heatmap to provide a better visual representation of the reads distrubution across our reference genes (provided in the reference_genes.bed file)
## Use the deepTools ** plotHeatmap ** function to create a heatmap following this criteria: 
## - Use the matrix from the previous point
## - Use the Blues colormap
## - Create 3 clusters within the heatmap according to the patterns of the reads distrubution across the files using heirarchical clustering
#?# Type the command you used to run it below: - 1.5
#?# Add a screenshot of the plot right after this code chunk using Rmarkdown syntaxis - 1 pt 

#plotHeatmap -m reference.tab --colorMap Blues --hclust 3 -o heatmap.svg


#?# Explain what you are looking at (Axes, colours, curves). Where are the marks located? What are the differences between the clusters? - 3 pts
# The x axis of the heatmap shows the gene distance in base pairs and the y axis shows which cluster each row corresponds to. The colours are shown in the legend on the right, with darker colors representing larger values. The curves show the data points that are being represented by the heatmap, which correspond to scores we calculated earlier. The marks are located in the darkest areas of the heatmap. The first cluster and second cluster have higher scores in general than the third cluster. Within each cluster the values are relatively similar to each other, which makes sense since we used hierarchical clustering. 
```

**Add screenshot here: ** (Use Rmarkdown syntax)

![](heatmap.svg)

``` bash
## Now the above heatmap was made with the ratio of ChIP to input. Repeat the process above, but this time using the raw bigwig files (not input-normalized). 
#?# Type the commands you used for this below - 1 pt
#computeMatrix reference-point -S H3K4me3_chr3.bw H3K27me3_chr3.bw H3K27ac_chr3.bw -R /projects/bmeg/A5/reference_genes.bed -o reference_notnorm.tab -a 1000 -b 1000
#plotHeatmap -m reference_notnorm.tab --colorMap Blues --hclust 3 -o heatmap_notnorm.svg
#?# Include a screenshot of this analysis, below this code block. - 1 pt
#?# How does this compare to the input-normalized data? Why do you think this is? - 1 pt
# The normalized data has more area covered in higher relative values (shown by the darker areas). This is because normalizing the data changes it to a common scale so we can better compare the data. So the relative ranges for the normalized data are the same as the non-normalized, but scaled to allow comparisons. 
```

**Add screenshot here: ** (Use Rmarkdown syntax)
![](heatmap_notnorm.svg)

## 3. Peak calling

Now we want to identify enriched regions of the genome for each of our
three histone marks. In order to get the enrichments, we will run the
**macs2** program to call the peaks for each epigenetic mark.

### a. Peak calling with macs2

``` bash
## Tip: Make sure to read the documentation (using the -h flag) for the *masc2 callpeak* command
## Run the callpeak command of the macs2 tool, once for each of H3K27ac, H3K27Me3, H3K4Me3
## Each time, use the input as the control 
## For H3K27Me3, you should call "broad" peaks because these signals tend to form longer domains and not acute peaks.
#?# Type the commands you used below: - 1.5 pt
# macs2 callpeak -t /projects/bmeg/A5/H3K4me3_chr3_subset.bam -c /projects/bmeg/A5/input_chr3_subset.bam -n H3K4me3 & macs2 callpeak -t /projects/bmeg/A5/H3K27me3_chr3_subset.bam -c /projects/bmeg/A5/input_chr3_subset.bam -n H3K27me3 --broad & macs2 callpeak -t /projects/bmeg/A5/H3K27ac_chr3_subset.bam -c /projects/bmeg/A5/input_chr3_subset.bam -n H3K27ac 
```

### b. Understanding the peaks

Macs2 calls the peaks by analyzing the enrichments of sequences in the
genome compared to the input. The algorithm uses different measures to
ensure the enrichments are statistically significant and make biological
sense based on the sequence length, the expected peaks (narrow or broad)
and other parameters. We will focus on the H3K4me3 mark to visualize how
are the peaks are identified.


    ## 1. Copy the H3K4me3 .narrowPeak file to your local computer
    ## 2. Open IGV with and load the hg38 reference 

    ## 3. Load the following files:
    ### a. H3K4me3 bigwig file that you created in section 1 
    ### b. Input  bigwig file that you created in section 1 
    ## Note: remember to autoscale each file track!

    ## 4. Go to this position: * chr3:44,608,952-44,670,849 *

    #?# Compare the input and H3K4me3 signal trakcs (a and b), are there regions that you consider to have an enriched signal in the H3K4me3 track compared to the input? If yes, how many?- 0.5 pt 
    # Yes, there are two sections.

    #?# Why do you consider those regions to have an enrichment? Would you need to confirm that the enrichment is statistically significant? If so, what do you think would be a way to test for the significance (you can mention the tools used on this assignment or explain conceptually what would you do to confirm an enrichment)? - 1 pt
    # These regions have an enrichment because they have peaks in the H3K4me3 track. The other regions are all flat. 

    ## 5. Load into the same session the following files:
    ### c. H3K4me3 normalized using the input bigWig file that was pre-computed and you used on section 2
    ### d. H3K4me3 .narrowPeak file that you created in section 3a
    ## Note: remember to autoscale each file track!

    ## 6. Go to the same position as in step 4

    #?# Does the region/s you thought had an enrichment show differences between the input and the H3K4me3 mark in the bamCompare (file c)? - 0.5 pt
    # Yes, it does.

    #?# Are these visually enriched regions you selected, found in file d (narrowPeak)? What does that mean? - 1 pt
    # Yes, it does. It means that there are peaks that represent signal enrichment in the selected regions. This means that the peaks identified in the H3K4me3 track were correct. 

    #?# Add a screenshot of your IGV session right after this chunk, using Rmarkdown syntax - 1pt

![](peaks_summary.png)

### c. Peak enrichments

For this assignment, we are working with 3 different epigenetic marks:
H3K4me3, H3K27me3 and H3K27ac. Each of these marks a different
transcription state (i.e., activation or repression). Thus, to better
visualize the differences in the called peaks for each of these
epigenetic marks you will create a heatmap plot (like the one you
created on part 2) using their called peaks files (.narrowPeak or
.broadPeak).

``` bash

### Create 3 heatmaps following the specifications you used on part 2 (one for the peaks called for each epigenetic mark, both containing data from all three tracks). Use the peak files of each of the epigenetic marks as reference files. Use ONLY the non-input normalized files: H3K27ac_norm.bw H3K27me3_norm.bw H3K4me3_norm.bw
#?# Write the commands you used to compute the matrices and create the heatmaps below: - 3 pt


#computeMatrix reference-point -S /projects/bmeg/A5/H3K4me3_norm.bw /projects/bmeg/A5/H3K27me3_norm.bw /projects/bmeg/A5/H3K27ac_norm.bw -R H3K27ac_peaks.narrowPeak -o h3k27ac_peaks_refenence.tab -a 1000 -b 1000
#plotHeatmap -m h3k27ac_peaks_refenence.tab --colorMap Blues --hclust 3 -o heatmap_h3k27ac_norm.svg

#computeMatrix reference-point -S /projects/bmeg/A5/H3K4me3_norm.bw /projects/bmeg/A5/H3K27me3_norm.bw /projects/bmeg/A5/H3K27ac_norm.bw -R H3K27me3_peaks.broadPeak -o h3k27me3_peaks_refenence.tab -a 1000 -b 1000
#plotHeatmap -m h3k27me3_peaks_refenence.tab --colorMap Blues --hclust 3 -o heatmap_h3k27me3_norm.svg

#computeMatrix reference-point -S /projects/bmeg/A5/H3K4me3_norm.bw /projects/bmeg/A5/H3K27me3_norm.bw /projects/bmeg/A5/H3K27ac_norm.bw -R H3K4me3_peaks.narrowPeak -o h3k4me3_peaks_refenence.tab -a 1000 -b 1000
#plotHeatmap -m  h3k4me3_peaks_refenence.tab --colorMap Blues --hclust 3 -o heatmap_h3k4me3_norm.svg
#?# Add screenshots of the 3 heatmaps you got using the epigenetic marks' peak files as reference files. Add them after this code chunk (using Rmarkdown syntax) in the following order: H3K4me3, H3K27ac, H3K27me3 - 1.5 pt
#Screenshots are added in order as listed above:
```

![](heatmap_h3k4me3_norm.svg) ![](heatmap_h3k27ac_norm.svg)
![](heatmap_h3k27me3_norm.svg)


    #?# Do you see an overlap between the peaks of different epigenetic marks? Which epigenetic marks? - 1 pt
    # Yes, there is some overlap between the peaks of H3K27ac and H3K27me3. This is because both epigenetic marks are associated with higher transcription activation and and enhancer activity. 

    #?# Why do you think these epigenetic marks overlap? - 1 pt
    # It is common for parts of epigenetic marks to overlap. It can be useful to compress information into genetic sequences, and is also a way to regulate gene expression since polypeptides with similar functions can be coupled together. 

    #?# Explain the pattern you see for the peaks of H3K27me3, do they look similar or different to the other two marks? Why do you think this happens? - 1pt
    # The graph is fairly flat until TSS, where a hump forms. All three clusters follow this pattern. This looks different to the other marks (H3K27 is more flat, and H3K4me3 is somewhat flat with larger humps). This could be because H3K27me3 is involved in transcriptional silencing. 

    #?# Why do you think the borders of the elements have such clearly-defined borders with respect to the ChIP signal? -1 pt
    # The borders of elements have such clearly-defined borders with respect to the ChIP signal is due to the mechanism ChIP-seq analysis uses. It allows us to get clearly defined borders.

# Assignment submission

Please knit your Rmd file to *github\_document* and include the link to
the .md file in your Canvas submission.

# Authors and contributions

Following completion of your assignment, please fill out this section
with the authors and their contributions to the assignment. If you
worked alone, only the author (e.g. your name and student ID) should be
included.

Authors: Alice Zhang (52721552) and Abhishek Dhir (87603866)

Contributions: Alice and Abhishek worked together asynchronously on the
assignment together and then pushed the combined changes in the end,
after a synchronous working session.
