% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterDeGene.R
\name{labelGenes}
\alias{labelGenes}
\title{Label genes with "UP", "DOWN" and "NOCHANGE" according to the pValue and
foldChange thresholds.}
\usage{
labelGenes(
  diffExpressionResult = NULL,
  filePath = NULL,
  save = FALSE,
  pValue = 0.05,
  foldChange = 2
)
}
\arguments{
\item{diffExpressionResult}{A DESeqResults object, a simple subclass of
DataFrame defined by DESeq2. It contains the columns: baseMean,
log2FoldChange, lfcSE, stat, pvalue and padj, and also includes
metadata columns of variable information. This input is the output
of function diffExpressionAnalysis. Default value is NULL.}

\item{filePath}{A character string path to the directory that the user want
to store the labelled genes based on the p value and fold change threshold.
It is recommended to create a new directory to store the output .csv file.
Default value is NULL. Should in the format: "/Path/to/the/directory".}

\item{save}{A boolean value to indicate if store the result into a file
at the file path. If set to TRUE, the output would be stored at the
filePath. Default value is FALSE.}

\item{pValue}{A positive double, which is the threshold set by the users used
to decide if the gene UP or DOWN expressed in different samples. Only
genes that has smaller p values than the threshold would be regarded as
significant genes (significantly up expressed or significantly down
expressed). Default value is 0.05. p Value is used to measure if a
difference is significant. The pValue should be a value between 0 and 1.}

\item{foldChange}{A positive integer, which is the threshold set by the users
used to decide how the expression levels of a gene change in different
samples. Only genes that has larger absolute values of fold change than
the threshold would be regarded as genes who expressed differently.
Default value is 2.
Fold change is a measure describing how much the expression levels of
genes change between different samples. The value of foldChange should
be a positive integer.}
}
\value{
Returns a data frame, which has an additional column called
"group" compared to the input data frame "diffExpressionResult",
It also contains the columns: baseMean, log2FoldChange,
lfcSE, stat, pvalue and padj, and also includes metadata columns of
variable information. Genes are labelled as up-expressed or down-expressed
or no significant change on expression levles.
}
\description{
This is a function that create a new column called "group" to the original
input "diffExpressionResult" and save the new data frame to a csv. In the new
column, each gene is labelled with "UP", "DOWN" and "NOCHANGE". Basically,
genes that shows significant higher expression levels would be labelled with
"UP", genes that shows significant lower expression levels would be labelled
with "DOWN". Otherwise, genes would be labelled with "NOCHANGE". The function
regards a gene as "UP" expressed when the gene has higher foldChange and
lower pValue than thresholds set by users. The function regards a gene as
"DOWN" expressed when the gene has lower -foldChange and lower pValue than
thresholds set by users. The rest of genes are labelled as "NOCHANGE".
}
\examples{
# Example 1:
# Using example differential expression analysis result
# (diffExpressionResult) available with the package
\dontrun{
dim(diffExpressionResult) # 18309 rows, 6 columns

# label genes based on the input threshold
labelGenesE <- labelGenes(diffExpressionResult,
                          pValue = 0.05,
                          foldChange = 2)
labelGenesE
}

}
\references{
Geistlinger L, Csaba G, Zimmer R (2016). “Bioconductor's EnrichmentBrowser:
seamless navigation through combined results of set- & network-based
enrichment analysis.” BMC Bioinformatics, 17, 45. doi:10.1186/s12859-016-0884-1.

R Core Team (2023). R: A Language and Environment for Statistical Computing. R Foundation
for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Wickham H., François R., Henry L., Müller K., Vaughan D. (2023).
dplyr: A Grammar of Data Manipulation. R package version 1.1.3,
https://CRAN.R-project.org/package=dplyr
}
