library(GEAnaly)

test_that("Check visDeAnaly works properly", {

  visResult <- visDeAnaly(labelledGenes, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "DiffExpresion_volcano_plot.png")), TRUE)
})

test_that("Checking for invalid input", {

  expect_message(visDeAnaly(filePath = getwd()))
  expect_message(visDeAnaly(labelledGenes, save = TRUE))
})

test_that("Check visCorrelationAnaly works properly", {

  visResult <- visCorrelationAnaly(geneCorResult, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "correlation_analysis_vis.png")), TRUE)
})

test_that("Checking for invalid input", {

  expect_message(visCorrelationAnaly(geneCorResult, save = TRUE))
  expect_message(visCorrelationAnaly(filePath = getwd()))
})

test_that("Check visEnrichAnaly works properly", {

  visResult <- visEnrichAnaly(enrichOutputList, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "enrich_analysis_vis.html")), TRUE)

})

test_that("Checking for invalid input", {

  expect_message(visEnrichAnaly(enrichOutputList, save = TRUE))
  expect_message(visEnrichAnaly(filePath = getwd()))
})

test_that("Check visEnrichAnalyLollipop works properly", {

  visResult <- visEnrichAnalyLollipop(enrichOutputList, filePath = getwd())

  expect_type(visResult, "list")
  expect_equal(file.exists(file.path(getwd(),
                                     "enrich_analysis_vis_Lollipop.png")), TRUE)
})

test_that("Checking for invalid input", {

  expect_message(visEnrichAnalyLollipop(enrichOutputList, save = TRUE))
  expect_message(visEnrichAnalyLollipop(filePath = getwd()))
})

# [END]
