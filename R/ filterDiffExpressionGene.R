

library(dplyr)
extractSignificantGene <- function(diffExpressionResult,
                                   filePath = "./significant_de_genes.csv",
                                   pValue = 0.1,
                                   foldChange = 1) {

  # Obtain the genes that are significantly differentially expressed
  # by threshold provided
  diffExpressionResult <- data.frame(diffExpressionResult)
  significantGenes <- diffExpressionResult %>%
    dplyr::filter(padj < pValue & log2FoldChange > foldChange)

  # Save the output as a file and return the output
    write.csv(significantGenes, file = filePath, quote = FALSE)
    return(significantGenes)
}


labelGenes <- function(diffExpressionResult,
                       filePath = "./de_genes_with_label.csv",
                       pValue = 0.05,
                       foldChange = 2) {

  # Convert the differential expression output to a dataframe
  diffExpressionResult <- data.frame(diffExpressionResult)

  # Generate labels for each gene according to the foldchange and p-value
  differentialExpressionResultlabel <- diffExpressionResult %>%
    dplyr::mutate(group = case_when(
      log2FoldChange >= foldChange & padj <= pValue ~ "UP",
      log2FoldChange <= -foldChange & padj <= pValue ~ "DOWN",
      TRUE ~ "NOCHANGE"
    ))

  # Save the output as a file or return the output

    write.csv(differentialExpressionResultlabel, file = filePath, quote = FALSE)
    return(differentialExpressionResultlabel)

}

