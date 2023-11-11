runGEAnaly <- function(geneCounts = NULL,
                       sampleInfor = NULL,
                       filePath = ".") {



  deResult <- diffExpressionAnalysis(geneCounts, sampleInfor)

  message("Saving differential expression analysis output to ", filePath)
  write.csv(deResult,
            file = file.path(filePath, "diffExpression_analysis_result.csv"),
            quote = FALSE)

  message("Filtering and labeling genes")
  significantGenes <- extractSignificantGene(deResult, filePath = filePath)
  labelledGenes <- labelGenes(deResult, filePath = filePath)

  message("Plotting and saving DE analysis to ", filePath)
  visDeAnaly(labelledGenes, filePath = filePath)

  message("Performing enrichment analysis on significant genes")
  enrichOutputList <- enrichAnalysis(significantGenes, filePath = filePath)

  message("Plotting enrichment analysis result")
  visEnrichAnaly(enrichOutputList, filePath = filePath)
  visEnrichAnalyLollipop(enrichOutputList, filePath = filePath)

  message("The whole pipeline has been finished!")
}

runGEAnalyCor <- function(geneCounts = NULL,
                          filePath = ".") {
  message("Performing gene correlation analysis")
  geneCorResult <- corrAnalysis(geneCounts, method = "pearson")

  message("Saving correlation analysis output to ", filePath)
  write.csv(geneCorResult,
            file = file.path(filePath, "correlation_analysis_result.csv"),
            quote = FALSE)

  message("Plotting and saving gene correlation analysis to ", filePath)
  visCorrelationAnaly(geneCorResult, filePath = filePath)
}
