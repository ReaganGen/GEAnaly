#' Filter out the genes that have significantly different expression levels by
#' the user-customized threshold
#'
#' This is a function that select genes that have significantly different
#' expression levels in different samples according to the p value and
#' fold change provided by the users. As long as the gene has smaller p value
#' than the p value threshold and has larger fold change value than the
#' fold change threshold inputted by the users, these genes are regarded as
#' genes with significantly different expression levels in different samples.
#' The function save these significant genes into a .csv file.
#'
#' @param diffExpressionResult A DESeqResults object, a simple subclass of
#'    DataFrame defined by DESeq2. It contains the columns: baseMean,
#'    log2FoldChange, lfcSE, stat, pvalue and padj, and also includes
#'    metadata columns of variable information. This input is the output
#'    of function diffExpressionAnalysis. Default value is NULL.
#' @param filePath A character string path to the directory that the user want
#'    to store the significant genes selected by the p value and fold change
#'    threshold. It is recommended to create a new directory to store the output
#'    .csv file. Default value is NULL. Should in the format:
#'    "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#' @param pValue A positive double, which is the threshold set by the
#'    users used to filter genes that are differentially expressed in different
#'    samples. Only genes that has smaller p values than the threshold would be
#'    regarded as significant genes. Default value is 0.05. pValue is used to
#'    measure if a difference is significant. The pValue should be a value
#'    between 0 and 1.
#' @param foldChange A positive integer, which is the threshold set by the users
#'    used to filter genes that are differentially expressed in different
#'    samples. Only genes that has larger absolute values of fold change than
#'    the threshold would be regarded as genes express differently. Default
#'    value is 2. Fold change is a measure describing how much the expression
#'    levels of genes change between different samples.
#'
#' @return Returns a data frame, which is a subset of the input dataframe
#'    "diffExpressionResult", containing the genes that have significantly
#'    different expression level in different samples (based on the input
#'    thresholds). It also contains the columns: baseMean, log2FoldChange,
#'    lfcSE, stat, pvalue and padj, and also includes
#'    metadata columns of variable information.
#'
#' @examples
#' # Example 1:
#' # Using example differential expression analysis result
#' # (diffExpressionResult) available with the package
#' \dontrun{
#' dim(diffExpressionResult) # 18309 rows, 6 columns
#'
#' # Filter out genes with significantly different expression levels
#' # This creates a significant_de_genes.csv file stores the result
#' # in the current working directory
#' significantGenesE <- extractSignificantGene(diffExpressionResult,
#'                                             filePath = getwd(),
#'                                             save = TRUE,
#'                                             pValue = 0.05,
#'                                             foldChange = 2)
#' significantGenesE
#' }
#'
#' @references
#' Geistlinger L, Csaba G, Zimmer R (2016). “Bioconductor's EnrichmentBrowser:
#' seamless navigation through combined results of set- & network-based
#' enrichment analysis.” BMC Bioinformatics, 17, 45. doi:10.1186/s12859-016-0884-1.
#'
#' R Core Team (2023). R: A Language and Environment for Statistical Computing. R Foundation
#' for Statistical Computing, Vienna, Austria. https://www.R-project.org/
#'
#' Wickham H., François R., Henry L., Müller K., Vaughan D. (2023).
#' dplyr: A Grammar of Data Manipulation. R package version 1.1.3,
#' https://CRAN.R-project.org/package=dplyr
#'
#' @export
#' @import dplyr
#' @import utils

extractSignificantGene <- function(diffExpressionResult = NULL,
                                   filePath = NULL,
                                   save = FALSE,
                                   pValue = 0.05,
                                   foldChange = 2) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.")
  } else {
    ;
  }

  if (is.null(diffExpressionResult) == TRUE) {
    stop("Please input a differential expression analysis result as the input,
         which should be the output of function \"diffExpressionAnalysis\"")
  } else {
    ;
  }

  if (typeof(pValue) != "double") {
    stop("Please input a p value which is a double type number.")
  } else {
    ;
  }

  if (typeof(foldChange) != "double") {
    stop("Please input a foldChange which is a double type number.")
  } else {
    ;
  }

  if (save == TRUE & typeof(filePath) != "character") {
    stop("Please input a character string as the path.e.g./Path/to/the/directory")
  } else {
    ;
  }

  if (foldChange <= 0) {
    stop("Please input a positive foldChange.")
  } else {
    ;
  }

  if (pValue >= 1 | pValue <= 0) {
    stop("Please input a p value which is between 0 and 1.")
  } else {
    ;
  }

  message("Extracting sinificant genes")

  # Convert the differential expression input as a dataframe
  diffExpressionResult <- data.frame(diffExpressionResult)

  # Obtain the significant genes that are differentially expressed
  # by threshold provided
  significantGenes <- diffExpressionResult %>%
    dplyr::filter(padj < pValue & abs(log2FoldChange) > foldChange)

  # Save the output as a file and return the output
  if (save == TRUE) {
    utils::write.csv(significantGenes,
                     file = file.path(filePath, "significant_de_genes.csv"),
                     quote = FALSE)
    return(significantGenes)
  } else {
    return(significantGenes)
  }
  message("Extraction finished!")
}

#' Label genes with "UP", "DOWN" and "NOCHANGE" according to the pValue and
#' foldChange thresholds.
#'
#' This is a function that create a new column called "group" to the original
#' input "diffExpressionResult" and save the new data frame to a csv. In the new
#' column, each gene is labelled with "UP", "DOWN" and "NOCHANGE". Basically,
#' genes that shows significant higher expression levels would be labelled with
#' "UP", genes that shows significant lower expression levels would be labelled
#' with "DOWN". Otherwise, genes would be labelled with "NOCHANGE". The function
#' regards a gene as "UP" expressed when the gene has higher foldChange and
#' lower pValue than thresholds set by users. The function regards a gene as
#' "DOWN" expressed when the gene has lower -foldChange and lower pValue than
#' thresholds set by users. The rest of genes are labelled as "NOCHANGE".
#'
#' @param diffExpressionResult A DESeqResults object, a simple subclass of
#'    DataFrame defined by DESeq2. It contains the columns: baseMean,
#'    log2FoldChange, lfcSE, stat, pvalue and padj, and also includes
#'    metadata columns of variable information. This input is the output
#'    of function diffExpressionAnalysis. Default value is NULL.
#' @param filePath A character string path to the directory that the user want
#'    to store the labelled genes based on the p value and fold change threshold.
#'    It is recommended to create a new directory to store the output .csv file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#' @param pValue A positive double, which is the threshold set by the users used
#'    to decide if the gene UP or DOWN expressed in different samples. Only
#'    genes that has smaller p values than the threshold would be regarded as
#'    significant genes (significantly up expressed or significantly down
#'    expressed). Default value is 0.05. p Value is used to measure if a
#'    difference is significant. The pValue should be a value between 0 and 1.
#' @param foldChange A positive integer, which is the threshold set by the users
#'    used to decide how the expression levels of a gene change in different
#'    samples. Only genes that has larger absolute values of fold change than
#'    the threshold would be regarded as genes who expressed differently.
#'    Default value is 2.
#'    Fold change is a measure describing how much the expression levels of
#'    genes change between different samples. The value of foldChange should
#'    be a positive integer.
#'
#' @return Returns a data frame, which has an additional column called
#'    "group" compared to the input data frame "diffExpressionResult",
#'    It also contains the columns: baseMean, log2FoldChange,
#'    lfcSE, stat, pvalue and padj, and also includes metadata columns of
#'    variable information. Genes are labelled as up-expressed or down-expressed
#'    or no significant change on expression levles.
#'
#' @examples
#' # Example 1:
#' # Using example differential expression analysis result
#' # (diffExpressionResult) available with the package
#' \dontrun{
#' dim(diffExpressionResult) # 18309 rows, 6 columns
#'
#' # label genes based on the input threshold
#' # This creates a de_genes_with_label.csv file stores the result
#' # in the current working directory
#' labelGenesE <- labelGenes(diffExpressionResult,
#'                           filePath = getwd(),
#'                           save = TRUE,
#'                           pValue = 0.05,
#'                           foldChange = 2)
#' labelGenesE
#' }
#'
#' @references
#' Geistlinger L, Csaba G, Zimmer R (2016). “Bioconductor's EnrichmentBrowser:
#' seamless navigation through combined results of set- & network-based
#' enrichment analysis.” BMC Bioinformatics, 17, 45. doi:10.1186/s12859-016-0884-1.
#'
#' R Core Team (2023). R: A Language and Environment for Statistical Computing. R Foundation
#' for Statistical Computing, Vienna, Austria. https://www.R-project.org/
#'
#' Wickham H., François R., Henry L., Müller K., Vaughan D. (2023).
#' dplyr: A Grammar of Data Manipulation. R package version 1.1.3,
#' https://CRAN.R-project.org/package=dplyr
#'
#' @export
#' @import dplyr
#' @import utils
#'
labelGenes <- function(diffExpressionResult = NULL,
                       filePath = NULL,
                       save = FALSE,
                       pValue = 0.05,
                       foldChange = 2) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.")
  } else {
    ;
  }

  if (is.null(diffExpressionResult) == TRUE) {
    stop("Please input a differential expression analysis result as the input,
         which should be the output of function \"diffExpressionAnalysis\"")
  } else {
    ;
  }

  if (typeof(pValue) != "double") {
    stop("Please input a p value which is a double type number.")
  } else {
    ;
  }

  if (typeof(foldChange) != "double") {
    stop("Please input a foldChange which is a double type number.")
  } else {
    ;
  }

  if (save == TRUE & typeof(filePath) != "character") {
    stop("Please input a character string as the path.e.g./Path/to/the/directory")
  } else {
    ;
  }

  if (foldChange <= 0) {
    stop("Please input a positive foldChange.")
  } else {
    ;
  }

  if (pValue >= 1 | pValue <= 0) {
    stop("Please input a p value which is between 0 and 1.")
  } else {
    ;
  }

  message("Labeling genes")

  # Convert the differential expression output to a data frame
  diffExpressionResult <- data.frame(diffExpressionResult)

  # Generate labels for each gene according to the fold change and p-value
  differentialExpressionResultlabel <- diffExpressionResult %>%
    dplyr::mutate(group = case_when(
      log2FoldChange >= foldChange & padj <= pValue ~ "UP",
      log2FoldChange <= -foldChange & padj <= pValue ~ "DOWN",
      TRUE ~ "NOCHANGE"
    ))

  # Save the output as a file or return the output
  if (save == TRUE) {
    utils::write.csv(differentialExpressionResultlabel,
                     file = file.path(filePath, "de_genes_with_label.csv"),
                     quote = FALSE)
  } else {
    ;
  }

  message("Labeling finished!")

  return(differentialExpressionResultlabel)

}

# [END]
