#' Perform gene differential expression analysis based on the DESeq2 package
#'
#' This is a function that perform the automatic gene differential expression
#' analysis given gene expression matrix and sample information about the
#' samples found in the gene expression matrix. This function perform
#' differential expression to identify difference in gene expressions in
#' different samples.
#'
#' @param geneCounts A data frame containing the matrix of un-normalized read
#'    counts of gene expression data, normally obtained from RNA-seq or other
#'    sequencing experiments. The count number in the matrix should be
#'    non-negative integers. The row names of the matrix should be the genes,
#'    and the column names should be the sample IDs. Default value is NULL.
#' @param sampleInfor A data frame containing the comparison table between
#'    sample IDs and sample group, e.g. control group or treated group.
#'    The first column is the IDs or names of the samples, the second column
#'    is the group that the sample belongs to. First column name shoud be
#'    "sample", second column name should be "condition". Default value is NULL.
#'
#' @return Returns a DESeqResults object, a simple subclass of DataFrame
#'    defined by DESeq2.This object contains the results columns: baseMean,
#'    log2FoldChange, lfcSE, stat, pvalue and padj, and also includes
#'    metadata columns of variable information. The lfcSE gives the standard
#'    error of the log2FoldChange. stat is the log2FoldChange divided by lfcSE,
#'    which is compared to a standard Normal distribution to generate
#'    a two-tailed pvalue.
#'
#' @examples
#' # Example 1:
#' # Using GeneCounts dataset available with package
#' dim(GeneCounts) # a n = 30 by d = 3 dataset
#'
#' # Example 2:
#'
#' # Example 3:
#'
#' @references
#' Love, M. I., Huber, W., &amp; Anders, S. (2014). Moderated estimation of fold
#' change and dispersion for RNA-seq data with deseq2. \emph{Genome Biology}, 15(12).
#' \href{https://doi.org/10.1186/s13059-014-0550-8}{link}
#'
#' Geistlinger L, Csaba G, Zimmer R (2016). “Bioconductor's EnrichmentBrowser:
#' seamless navigation through combined results of set- & network-based
#' enrichment analysis.” BMC Bioinformatics, 17, 45. doi:10.1186/s12859-016-0884-1.
#'
#' @export
#' @import DESeq2

diffExpressionAnalysis <- function(geneCounts = NULL,
                                   sampleInfor = NULL) {

  # Performing checks of user input
  if (is.null(geneCounts) == TRUE) {
    stop("Please input a gene count matrix as the input.")
  } else {
    ;
  }

  if (is.null(sampleInfor) == TRUE) {
    stop("Please input a sample information table as the input.")
  } else {
    ;
  }

  if (identical(colnames(geneCounts), sampleInfor$sample) != TRUE) {
    stop("The content and order of the sample names in Gene Counts table
         and Sample Infotmation table should be the same.")
  } else {
    ;
  }

  if (all((geneCounts - floor(geneCounts)) == 0) == FALSE) {
    stop("The count in your Gene Count table should be integers.")
  } else {
    ;
  }

  # Start analysis
  # Provide information about the analysis to the users
  message("Performing gene differential expression analysis")

  # Change the sampleInfor to factor type
  sampleInfor$sample <- factor(sampleInfor$sample)
  sampleInfor$condition <- factor(sampleInfor$condition)

  # Remove genes that has no expression in all samples
  geneCounts <- geneCounts[rowSums(geneCounts) != 0, ]

  # Perform differential expression analysis
  deseqDataset <- DESeq2::DESeqDataSetFromMatrix(countData = geneCounts,
                                         colData = sampleInfor,
                                         design = ~condition)

  # Access the analysis result and return
  deseqDataset <- DESeq2::DESeq(deseqDataset)
  result <- DESeq2::results(deseqDataset)
  return(result)
}

# [END]
