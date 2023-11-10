library(stats)
corrAnalysis <- function(geneCounts,
                         filePath = "./correlation_analysis.csv",
                         method = "pearson") {

  geneCountsRemoved <- geneCounts[rowSums(geneCounts[])>0,]
  geneCountsRemoved <- t(geneCountsRemoved)
  geneCountsLog <- log(geneCountsRemoved + 1)
  geneCor <- stats::cor(geneCountsLog,
                        method = method)
  diag(geneCor) <- 0
  write.csv(geneCor, file = filePath, quote = FALSE)
  return(geneCor)
}
