library(GEAnaly)

test_that("The function perform enrichment analysis properly with g_SCS", {

  enrichOutputListE <- enrichAnalysis(significantGenes,
                                     pvalueCutoff = 0.05,
                                     correctionMethod = "g_SCS",
                                     filePath = getwd())

  expect_type(enrichOutputListE, "list")
  expect_type(enrichOutputListE$enrichmentVis, "list")
  expect_type(enrichOutputListE$gProfilerResult, "list")
  expect_equal(length(enrichOutputListE), 2)
  expect_equal(length(enrichOutputListE$enrichmentVis), 7)
  expect_equal(length(enrichOutputListE$gProfilerResult), 2)
})

test_that("The function perform enrichment analysis properly with fdr", {

  enrichOutputListE <- enrichAnalysis(significantGenes,
                                     pvalueCutoff = 0.1,
                                     correctionMethod = "fdr",
                                     filePath = getwd())

  expect_type(enrichOutputListE, "list")
  expect_type(enrichOutputListE$enrichmentVis, "list")
  expect_type(enrichOutputListE$gProfilerResult, "list")
  expect_equal(length(enrichOutputListE), 2)
  expect_equal(length(enrichOutputListE$enrichmentVis), 7)
  expect_equal(length(enrichOutputListE$gProfilerResult), 2)
})

test_that("The function perform enrichment analysis properly with bonferroni", {

  enrichOutputListE <- enrichAnalysis(significantGenes,
                                     pvalueCutoff = 0.1,
                                     correctionMethod = "bonferroni",
                                     filePath = getwd())

  expect_type(enrichOutputListE, "list")
  expect_type(enrichOutputListE$enrichmentVis, "list")
  expect_type(enrichOutputListE$gProfilerResult, "list")
  expect_equal(length(enrichOutputListE), 2)
  expect_equal(length(enrichOutputListE$enrichmentVis), 7)
  expect_equal(length(enrichOutputListE$gProfilerResult), 2)
})

test_that("Checking for invalid input", {

  expect_error(enrichAnalysis(significantGenes,
                              pvalueCutoff = 0.05,
                              correctionMethod = "g_SCS"))
  expect_error(enrichAnalysis(pvalueCutoff = 0.1,
                              correctionMethod = "bonferroni",
                              filePath = getwd()))
  expect_error(enrichAnalysis(significantGenes,
                              pvalueCutoff = "a",
                              correctionMethod = "bonferroni",
                              filePath = getwd()))
  expect_error(enrichAnalysis(significantGenes,
                              pvalueCutoff = 0.05,
                              correctionMethod = 3,
                              filePath = getwd()))
  expect_error(enrichAnalysis(significantGenes,
                              pvalueCutoff = 10,
                              correctionMethod = "g_SCS",
                              filePath = getwd()))
})

# [END]
