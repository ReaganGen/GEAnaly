library(GEAnaly)

test_that("The function works properly", {

  pipelineResult <- runGEAnaly(geneCountsDiffExpression,
                               sampleInforDiffExpression,
                               filePath = getwd())

  expect_identical(pipelineResult, NULL)
})

test_that("Checking for invalid input", {

  invalidDEMatrix1 <- geneCountsDiffExpression[ , 1:10]

  invalidDEMatrix2 <- geneCountsDiffExpression
  invalidDEMatrix2[1, 1] <- 2.5

  expect_error(runGEAnaly(invalidDEMatrix1,
                          sampleInforDiffExpression,
                          filePath = getwd()))

  expect_error(runGEAnaly(invalidDEMatrix2,
                          sampleInforDiffExpression,
                          filePath = getwd()))

  expect_error(runGEAnaly(geneCountsDiffExpression,
                          sampleInforDiffExpression,
                          filePath = 1))

  expect_error(runGEAnaly(sampleInforDiffExpression,
                          filePath = getwd()))
  expect_error(runGEAnaly(geneCountsDiffExpression,
                          filePath = getwd()))
  expect_error(runGEAnaly(geneCountsDiffExpression,
                          sampleInforDiffExpression))
})

test_that("The function works properly", {

  pipelineResultC <- runGEAnalyCor(geneCountsCorrelation,
                                   filePath = getwd())

  expect_identical(pipelineResultC, NULL)
})

test_that("Checking for invalid input", {

  expect_error(runGEAnalyCor(geneCountsCorrelation,
                             filePath = 1))

  expect_error(runGEAnalyCor(filePath = getwd()))
  expect_error(runGEAnalyCor(geneCountsCorrelation))
})

# [END]
