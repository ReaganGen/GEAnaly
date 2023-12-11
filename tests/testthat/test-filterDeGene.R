library(GEAnaly)

test_that("Check if extractSignificantGene function perform filtering properly", {

  significantGenesE <- extractSignificantGene(diffExpressionResult,
                                              filePath = getwd(),
                                              pValue = 0.05,
                                              foldChange = 2)

  expect_type(significantGenesE, "list")
  expect_equal(length(row.names(significantGenesE)), 229)
  expect_equal(length(colnames(significantGenesE)), 6)
})

test_that("Checking for invalid input", {

  expect_message(extractSignificantGene(filePath = getwd(),
                                      pValue = 0.05,
                                      foldChange = 2))
  expect_message(extractSignificantGene(diffExpressionResult,
                                      filePath = 1,
                                      save = TRUE,
                                      pValue = 0.05,
                                      foldChange = 2))
  expect_message(extractSignificantGene(diffExpressionResult,
                                      filePath = getwd(),
                                      pValue = "0.05",
                                      foldChange = 2))
  expect_message(extractSignificantGene(diffExpressionResult,
                                      filePath = getwd(),
                                      pValue = 12,
                                      foldChange = 2))
  expect_message(extractSignificantGene(diffExpressionResult,
                                      filePath = getwd(),
                                      pValue = 0.05,
                                      foldChange = "A"))
  expect_message(extractSignificantGene(diffExpressionResult,
                                      filePath = getwd(),
                                      pValue = 0.05,
                                      foldChange = -2))
})

test_that("Check if labelGenes function labelling properly", {

  labelGenesE <- labelGenes(diffExpressionResult,
                                  filePath = getwd(),
                                  pValue = 0.05,
                                  foldChange = 2)

  expect_type(labelGenesE, "list")
  expect_equal(length(row.names(labelGenesE)), 18309)
  expect_equal(length(colnames(labelGenesE)), 7)
})

test_that("Checking for invalid input", {

  expect_message(labelGenes(diffExpressionResult,
                          pValue = 0.05,
                          foldChange = 2))
  expect_message(labelGenes(filePath = getwd(),
                          pValue = 0.05,
                          foldChange = 2))
  expect_message(labelGenes(diffExpressionResult,
                          filePath = 1,
                          pValue = 0.05,
                          foldChange = 2))
  expect_message(labelGenes(diffExpressionResult,
                          filePath = getwd(),
                          pValue = "0.05",
                          foldChange = 2))
  expect_error(labelGenes(diffExpressionResult,
                          filePath = getwd(),
                          pValue = 12,
                          foldChange = 2))
  expect_message(labelGenes(diffExpressionResult,
                          filePath = getwd(),
                          pValue = 0.05,
                          foldChange = "A"))
  expect_message(labelGenes(diffExpressionResult,
                          filePath = getwd(),
                          pValue = 0.05,
                          foldChange = -2))
})

# [END]
