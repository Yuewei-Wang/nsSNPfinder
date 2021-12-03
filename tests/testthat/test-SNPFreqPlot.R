context("Generate SNP position versus number of different nucleotides plot
         within the gene range")
library(nsSNPfinder)

test_that("invalid input1: wrong chromosome name", {
  chrName = 25
  startPosition = 232
  endPosition = 272

  expect_error(SNPplot <- SNPFreqPlot(chrName = chrName,
                                      startPosition = startPosition,
                                      endPosition = endPosition))
})


test_that("invalid input2: invalid gene range", {
  chrName = 22
  startPosition = -1
  endPosition = 272

  expect_error(SNPplot <- SNPFreqPlot(chrName = chrName,
                                      startPosition = startPosition,
                                      endPosition = endPosition))
})


test_that("invalid input3: gene range too long", {
  chrName = 1
  startPosition = 2321253
  endPosition = 2391707

  expect_error(SNPplot <- SNPFreqPlot(chrName = chrName,
                                      startPosition = startPosition,
                                      endPosition = endPosition))
})


test_that("invalid input4: no encoding protein within the range", {
  chrName = "X"
  startPosition = 2794535
  endPosition = 2808377

  expect_error(testResult <- SNPFreqPlot(chrName = chrName,
                                         startPosition = startPosition,
                                         endPosition = endPosition))
})

# [END]
