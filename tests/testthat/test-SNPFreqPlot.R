context("Generate SNP position versus number of different nucleotides plot
         within the gene range")
library(nsSNPfinder)

test_that("invalid input1: wrong chromosome name", {
  chrName = 25
  start_position = 232
  end_position = 272

  expect_error(SNPplot <- SNPFreqPlot(chrName = chrName,
                                      start_position = start_position,
                                      end_position = end_position))
})


test_that("invalid input2: invalid gene range", {
  chrName = 22
  start_position = -1
  end_position = 272

  expect_error(SNPplot <- SNPFreqPlot(chrName = chrName,
                                      start_position = start_position,
                                      end_position = end_position))
})


test_that("invalid input3: gene range too long", {
  chrName = 1
  start_position = 2321253
  end_position = 2391707

  expect_error(SNPplot <- SNPFreqPlot(chrName = chrName,
                                      start_position = start_position,
                                      end_position = end_position))
})


test_that("invalid input4: no encoding protein within the range", {
  chrName = "X"
  start_position = 2794535
  end_position = 2808377

  expect_error(testResult <- SNPFreqPlot(chrName = chrName,
                                         start_position = start_position,
                                         end_position = end_position))
})

# [END]
