library(stats)
corrAnalysis <- function(geneCounts,
                         method = "pearson") {

  geneCountsRemoved <- geneCounts[rowSums(geneCounts[])>0,]
  geneCountsRemoved <- t(geneCountsRemoved)
  geneCountsLog <- log(geneCountsRemoved + 1)
  geneCor <- stats::cor(geneCountsLog,
                        method = method)
  diag(geneCor) <- 0
  return(geneCor)
}
