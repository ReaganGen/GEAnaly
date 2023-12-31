% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{labelledGenes}
\alias{labelledGenes}
\title{A data frame that has the output of differential expression analysis and
label each gene with "UP", "DOWN" and "NOCHANGE"}
\format{
A matrix with 18309 rows and 7 columns, row names are gene names:
\describe{
\item{baseMean}{The mean of normalized counts of all samples, normalizing
for sequencing depth}
\item{log2FoldChange}{The effect size estimate, which tells about how much
the gene's expression seems to have changed in
different samples}
\item{lfcSE}{Standard error of the log2FoldChange.}
\item{stat}{log2FoldChange divided by lfcSE}
\item{pvalue}{Use stat to compare to a standard Normal distribution to
generate a two-tailed pvalue, which tells about whether the
difference is significant}
\item{padj}{Adjusted p value}
\item{group}{Up-expressed genes are labelled with "UP", down-expressed genes
are labelled with "DOWN". Otherwise, genes are labelled
with "NOCHANGE".}
}
}
\source{
Raw data obtained from Gene Expression Omnibus, generated by
GEAnaly::diffExpressionAnalysis and GEAnaly::labelGenes.
\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103}{link}
}
\usage{
labelledGenes
}
\description{
Basically similar to the diffExpressionResult, with an additional column
to label each gene. Genes that shows significant higher expression levels
are labelled with "UP", genes that shows significant lower expression
levels are labelled with "DOWN". Otherwise, genes are labelled
with "NOCHANGE".
}
\examples{
\dontrun{
 labelledGenes
}
}
\keyword{datasets}
