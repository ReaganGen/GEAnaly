#' Perform gene correlation analysis on the gene expression matrix data, using
#' method indicated by users
#'
#' This is a function that can conduct genetic correlation analysis to identify
#' pairs of genes that have a relationship in their expression. It calculates
#' the correlation coefficients between each pair of genes and store them into a
#' matrix. Different ways for the calculation of the correlation coefficients
#' are provided.
#'
#' @param geneCounts A data frame containing the matrix of gene expression data,
#'    normally obtained from RNA-seq or other sequencing experiments.
#'    The row names of the matrix should be the genes, and the column names
#'    should be the sample IDs. Default value is NULL.
#' @param method A character string indicating which correlation coefficient is
#'    to be computed. One of "pearson" (default), "kendall", or "spearman". If
#'    you are not sure which coefficient to use, just use "pearson".
#'    \href{https://www.R-project.org/}{Link for more detail about correlation
#'    coefficient}
#' @param filePath A character string path to the directory that the user want
#'    to store the correlation analysis result. It is recommended to create a
#'    new directory to store the output .csv file. Default value is NULL.
#'    Should in the format: "/Path/to/the/directory/".
#'
#' @return A correlation coefficient matrix that has n * n dimensions (n is the
#'    number of genes). Each cell contains the correlation coefficient between
#'    corresponding genes. The function would save the matrix into a .csv file.
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
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' \href{https://www.R-project.org/}{link}.
#'
#' Geistlinger L, Csaba G, Zimmer R (2016). “Bioconductor's EnrichmentBrowser:
#' seamless navigation through combined results of set- & network-based
#' enrichment analysis.” BMC Bioinformatics, 17, 45. doi:10.1186/s12859-016-0884-1.
#'
#' @export
#' @import utils
#' @importFrom stats cor

corrAnalysis <- function(geneCounts = NULL,
                         method = "pearson",
                         filePath = NULL) {

  # Performing checks of user input
  if (is.null(geneCounts) == TRUE) {
    stop("Please input a gene count matrix as the input.")
  } else {
    ;
  }

  if (is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.")
  } else {
    ;
  }

  # Remove genes that has no expression in all samples
  geneCountsRemoved <- geneCounts[rowSums(geneCounts[])>0,]

  # Transpose the matrix and calculate the log value of gene counts and calculate
  # the correlation between genes
  message("Performing gene correlation analysis")
  geneCountsRemoved <- t(geneCountsRemoved)
  geneCountsLog <- log(geneCountsRemoved + 1)
  geneCor <- stats::cor(geneCountsLog,
                        method = method)

  # Set the diagonal of the correlation matrix to 0
  diag(geneCor) <- 0

  # Save the output as a file or return the output
  message("Saving correlation analysis output to ", filePath)
  utils::write.csv(geneCor,
                   file = file.path(filePath, "correlation_analysis_result.csv"),
                   quote = FALSE)

  return(geneCor)
}

# [END]
