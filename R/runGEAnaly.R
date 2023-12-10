#' Execute the differential expression analysis and enrichment analysis and
#' their visualizations in a single call of a function.
#'
#' This is a function that can execute the gene expression analysis
#' pipeline. It accept a gene counts matrix and sample information, and perform
#' differential expression analysis, after that, genes with significantly
#' different expression levels are selected and labelled. Then a volcano plot is
#' used for the visualization of the differential expression analysis. Subsequently,
#' those selected significant genes are sent for an enrichment analysis. Then
#' a Manhattan plot and a Lollipop plot are generated for the visualization of
#' enrichment analysis. All steps are run by default settings.
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
#' @param filePath A character string path to the directory that the user want
#'    to store the output files. It is recommended to create a new directory.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#'
#' @return A NULL, save the output files to "filePath"
#'
#' @examples
#' # Example 1:
#' # Using gene counts matrix (geneCountsDiffExpression) and the sample
#' # information (sampleInforDiffExpression) available with the package
#' \dontrun{
#' dim(geneCountsDiffExpression) # 11000 rows, 126 columns
#' dim(sampleInforDiffExpression) # 126 rows, 2 columns
#'
#' # Perform the whole gene expression analysis pipeline on the input gene
#' # counts matrix
#' runGEAnaly(geneCountsDiffExpression,
#'            sampleInforDiffExpression,
#'            filePath = getwd(),
#'            save = TRUE)
#' # You should see all visualizations and all results from each stage of the
#' # pipeline in the current working directory, which is the path from getwd()
#' }
#'
#' @export

runGEAnaly <- function(geneCounts = NULL,
                       sampleInfor = NULL,
                       filePath = NULL,
                       save = TRUE) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.")
  } else {
    ;
  }

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

  if (save == TRUE & typeof(filePath) != "character") {
    stop("Please input a character string as the path.e.g./Path/to/the/directory")
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

  # Perform differential expression analysis and save the output
  deResult <- diffExpressionAnalysis(geneCounts, sampleInfor)
  message("Saving differential expression analysis output to ", filePath)
  if (save == TRUE) {
    write.csv(deResult,
              file = file.path(filePath, "diffExpression_analysis_result.csv"),
              quote = FALSE)
  } else {
    ;
  }


  # Filter out the genes that have significantly different expression levels
  sigGenes <- extractSignificantGene(deResult, filePath = filePath, save = save)

  # Label genes with "UP", "DOWN" and "NOCHANGE"
  labelGenes <- labelGenes(deResult, filePath = filePath, save = save)

  # Visualize the differential expression analysis
  visDeAnaly(labelGenes, filePath = filePath, save = save)

  # Perform enrichment analysis
  enrichOutputListE <- enrichAnalysis(sigGenes, filePath = filePath, save = save)

  # Visualize enrichment analysis result as Manhattan plot and Lollipop plot
  visEnrichAnaly(enrichOutputListE, filePath = filePath, save = save)
  visEnrichAnalyLollipop(enrichOutputListE, filePath = filePath, save = save)

  message("The whole pipeline has been finished!")

  return(invisible(NULL))
}

#' Execute the gene correlation analysis and its visualization in a single
#' call of a function
#'
#' This is a function that can perform gene correlation analysis and visualize
#' the result. It accepts a gene counts matrix, and perform gene correlation
#' analysis, after that, the result is visualized by a heatmap. All steps are
#' dun by default settings.
#'
#' @param geneCounts A data frame containing the matrix of gene expression data,
#'    normally obtained from RNA-seq or other sequencing experiments.
#'    The row names of the matrix should be the genes, and the column names
#'    should be the sample IDs. Default value is NULL.
#' @param filePath A character string path to the directory that the user want
#'    to store the output files. It is recommended to create a new directory.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#'
#' @return The correlation analysis result as a heatmap, if the save set as TRUE,
#'    save the output files to "filePath"
#'
#' @examples
#' # Example 1:
#' # Using gene counts matrix (geneCountsCorrelation) available with the package
#' \dontrun{
#' dim(geneCountsCorrelation) # 105 rows, 41 columns
#'
#' # Perform the gene correlation analysis pipeline on the input gene
#' # counts matrix
#' runGEAnalyCor(geneCountsCorrelation,
#'               filePath = getwd(),
#'               save = TRUE)
#' # You should see visualizations and results from each stage of the
#' # pipeline in the current working directory, which is the path from getwd()
#' }
#'
#' @export

runGEAnalyCor <- function(geneCounts = NULL,
                          filePath = NULL,
                          save = TRUE) {

  # Performing checks of user input
  if (is.null(geneCounts) == TRUE) {
    stop("Please input a gene count matrix as the input.")
  } else {
    ;
  }

  if (save == TRUE & is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.")
  } else {
    ;
  }

  if (save == TRUE & typeof(filePath) != "character") {
    stop("Please input a character string as the path.e.g./Path/to/the/directory")
  } else {
    ;
  }

  # Perform gene correlation analysis
  geneCor <- corrAnalysis(geneCounts,
                                filePath = filePath, save = save)

  # Visualize the correlation analysis result as a heatmap
  corPlot <- visCorrelationAnaly(geneCor, filePath = filePath, save = save)

  message("The correlation analysis has been finished!")

  return(corPlot)
}

# [END]
