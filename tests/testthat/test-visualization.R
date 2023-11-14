library(GEAnaly)

test_that("Check visDeAnaly works properly", {

  visResult <- visDeAnaly(labelledGenes, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "DiffExpresion_volcano_plot.png")), TRUE)
})

test_that("Checking for invalid input", {

  expect_error(visDeAnaly(filePath = getwd()))
  expect_error(visDeAnaly(labelledGenes))
})

test_that("Check visCorrelationAnaly works properly", {

  visResult <- visCorrelationAnaly(geneCorResult, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "correlation_analysis_vis.png")), TRUE)
})

test_that("Checking for invalid input", {

  expect_error(visCorrelationAnaly(geneCorResult))
  expect_error(visCorrelationAnaly(filePath = getwd()))
})

test_that("Check visEnrichAnaly works properly", {

  visResult <- visEnrichAnaly(enrichOutputList, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "enrich_analysis_vis.html")), TRUE)

})

test_that("Checking for invalid input", {

  expect_error(visEnrichAnaly(enrichOutputList))
  expect_error(visEnrichAnaly(filePath = getwd()))
})

test_that("Check visEnrichAnalyLollipop works properly", {

  visResult <- visEnrichAnalyLollipop(enrichOutputList, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "enrich_analysis_vis_Lollipop.png")), TRUE)
})

test_that("Checking for invalid input", {

  expect_error(visEnrichAnalyLollipop(enrichOutputList))
  expect_error(visEnrichAnalyLollipop(filePath = getwd()))
})

# [END]
