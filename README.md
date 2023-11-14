
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEAnaly

<!-- badges: start -->
<!-- https://www.codefactor.io/repository/github/ReaganGen/GEAnaly/issues -->

[![GitHub
issues](https://img.shields.io/github/issues/ReaganGen/GEAnaly)](https://github.com/ReaganGen/GEAnaly/issues)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
![GitHub language
count](https://img.shields.io/github/languages/count/ReaganGen/GEAnaly)
![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/y/ReaganGen/GEAnaly/master)

<!-- https://shields.io/category/license -->
<!-- badges: end -->

## Description

`GEAnaly` is an R package that implements a pipeline, which integrates
gene differential expression analysis, gene enrichment analysis and gene
correlation analysis together for the analysis of gene expression matrix
data and visualization of the results of each analysis in the pipeline.
The whole analysis pipeline can be executed in a single command
efficiently.

The pipeline of the package is basically: Gene expression matrix data →
Gene differential expression analysis → Select genes that are
differentially expressed → Enrichment analysis for those differentially
expressed genes, correlation analysis would also be applied on Gene
expression matrix data. Apart from that, for each analysis,
visualizations of results are generated.

The biological data analyzed by GEAnaly would be the gene expression
data (usually generated from RNA-sequencing experiments), which are the
expression matrix of genes in different samples, and the table that
records the information about the experimental samples (like sample
names and if they are controlled sample or treated sample). In the
expression matrix file, the names of rows are gene names, and the names
of columns are sample names. The expression matrix file records the
estimated count of sequencing reads of each gene (un-normalized counts).

Currently, performing a complete bioinformatics pipeline to analyze gene
expression matrix data requires programming skills and experience on
processing data. However, it is hard for biologists to learn programming
in a short period of time. Therefore, GEAnaly aims to make the analysis
of gene expression data as easy as possible without adjusting file
formats between different analysis. Also, this package provides a very
simplified data analysis experience for biologists with no programming
experiences, eliminating the manual construction of the workflow. In
addition, users can use very simple commands to visualize the analysis
results, once again simplifying the analysis process, and improving
analysis efficiency. The `GEAnaly` package was developed using
`R version 4.3.1 (2023-06-16)`,
`Platform: aarch64-apple-darwin20 (64-bit)` and
`Running under: macOS Sonoma 14.1`.

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("ReaganGen/GEAnaly", build_vignettes = TRUE)
library("GEAnaly")
```

To run the Shiny app:

``` r
Under construction
```

## Overview

``` r
# To list all available functions in the package, use command:
ls("package:GEAnaly")

# To list all example datasets available in the package:
data(package = "GEAnaly") 

# To see the tutorial of GEAnaly
browseVignettes("GEAnaly")
```

`GEAnaly` contains 11 functions:

3.  ***diffExpressionAnalysis***: Perform differential expression
    analysis on gene expression matrix data.

4.  ***corrAnalysis***: Perform gene expression correlation analysis on
    gene expression matrix data.

5.  ***extractSignificantGene***: Filter out the genes that have
    significantly different expression levels in different samples based
    on the threshold provided and the results of differential expression
    analysis.

6.  ***labelGenes***: Label genes with “UP”, “DOWN” and “NOCHANGE”
    according to the pValue and foldChange thresholds provided.

7.  ***enrichAnalysis***: Perform genes’ functional enrichment analysis,
    which is the procedure of identifying functions that are over- or
    under-represented among a set of genes.

8.  ***visDeAnaly***: Generate a volcano plot for the differential
    expression analysis result.

9.  ***visCorrelationAnaly***: Generate a heatmap plot for the gene
    correlation analysis result.

10. ***visEnrichAnaly***: Generate a Manhattan plot for the enrichment
    analysis result.

11. ***visEnrichAnalyLollipop***: Generate a Lollipop plot for the
    enrichment analysis result (show the top 10 biological pathways with
    top 10 highest gene ratio values).

The package also contains 8 datasets.

``` r
# use the following command to list all example datasets available in the package:
data(package = "GEAnaly") 

# Use ? to check details of the data.
?diffExpressionResult
```

Refer to package vignettes for more details. An overview of the package
is illustrated below:
