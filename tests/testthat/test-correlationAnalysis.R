library(GEAnaly)

test_that("Check the function can perform correlation analysis properly with
          pearson correlation coefficient", {

  corrResult <- corrAnalysis(geneCountsCorrelation,
                             filePath = getwd(),
                             method = "pearson")

  expect_type(corrResult, "double")
  expect_equal(length(row.names(corrResult)), 103)
  expect_equal(length(colnames(corrResult)), 103)
  expect_identical(colnames(corrResult), row.names(corrResult))
})

test_that("Check the function can perform correlation analysis properly with
          kendall correlation coefficient", {

  corrResult <- corrAnalysis(geneCountsCorrelation,
                             filePath = getwd(),
                             method = "kendall")

  expect_type(corrResult, "double")
  expect_equal(length(row.names(corrResult)), 103)
  expect_equal(length(colnames(corrResult)), 103)
  expect_identical(colnames(corrResult), row.names(corrResult))
})

test_that("Check the function can perform correlation analysis properly with
          spearman correlation coefficient", {

  corrResult <- corrAnalysis(geneCountsCorrelation,
                             filePath = getwd(),
                             method = "spearman")

  expect_type(corrResult, "double")
  expect_equal(length(row.names(corrResult)), 103)
  expect_equal(length(colnames(corrResult)), 103)
  expect_identical(colnames(corrResult), row.names(corrResult))
})

test_that("Checking for invalid input", {

  # Missing input data
  expect_error(corrAnalysis(filePath = getwd()))
  expect_error(corrAnalysis(geneCounts = geneCountsCorrelation))
  expect_error(corrAnalysis(geneCountsCorrelation,
                            filePath = 2,
                            method = "spearman"))
})

# [END]
