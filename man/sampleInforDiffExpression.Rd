% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sampleInforDiffExpression}
\alias{sampleInforDiffExpression}
\title{A comparison table for sample IDs and sample conditions, supplementary
information for the gene count matrix used in differential expression analysis.}
\format{
A data frame with 126 rows and 2 columns, row names are sample names:
\describe{
\item{sample}{List all sample names}
\item{condition}{Indicate the condition for each sample, with Covid or not}
}
}
\source{
Gene Expression Omnibus.
\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103}{link}
}
\usage{
sampleInforDiffExpression
}
\description{
A data frame containing the comparison table between
sample IDs and sample conditions, which are Covid and Noncovid.
The first column (sample) is the IDs or names of the samples, the second
column (condition) is the group that the sample belongs to.
}
\examples{
\dontrun{
 sampleInforDiffExpression
}
}
\keyword{datasets}
