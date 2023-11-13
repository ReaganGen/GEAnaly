library(GEAnaly)

test_that("The function perform differential expression analysis properly", {

  # Use data available with the package
  deResult <- diffExpressionAnalysis(geneCountsDiffExpression,
                                     sampleInforDiffExpression)

  expect_s4_class(deResult, "DESeqResults")
  expect_length(deResult, 6)
  expect_identical(colnames(deResult), c("baseMean",
                                         "log2FoldChange",
                                         "lfcSE",
                                         "stat",
                                         "pvalue",
                                         "padj"))
  expect_equal(length(rownames(deResult)), 10393)
})

test_that("Checking for invalid input", {

  invalidDEMatrix1 <- geneCountsDiffExpression[ , 1:10]

  invalidDEMatrix2 <- geneCountsDiffExpression
  invalidDEMatrix2[1, 1] <- 2.5

  # Information in the comparison table not match expression data
  expect_error(diffExpressionAnalysis(invalidDEMatrix1,
                                      sampleInforDiffExpression))

  # Counts in the expression matrix are not integers
  expect_error(diffExpressionAnalysis(invalidDEMatrix2,
                                      sampleInforDiffExpression))

  # Missing input data
  expect_error(diffExpressionAnalysis(sampleInfor = sampleInforDiffExpression))
  expect_error(diffExpressionAnalysis(geneCounts = geneCountsDiffExpression))
})

# [END]
