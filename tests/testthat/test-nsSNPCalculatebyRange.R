context("Calculate nsSNP percentage within the gene range")
library(nsSNPfinder)

test_that("valid input: approper arguments in chromosome 1", {
  chrName = 1
  startPosition = 2321253
  endPosition = 2391707

  testResult <- nsSNPCalculatebyRange(chrName = chrName,
                                      startPosition = startPosition,
                                      endPosition = endPosition)

  expect_type(testResult, "list")
  expect_length(testResult, 4)
  expect_s3_class(testResult, "data.frame")
  expect_equal(testResult$geneName[1], "MORN1")
  expect_identical(trunc(testResult$lengths[1]), 70302)
  expect_identical(trunc(testResult$nsSNPs[1]), 3072)
  expect_equal(as.numeric(testResult$percent[1]), 0.0437)
  expect_equal(testResult$geneName[2], "MORN1")
  expect_identical(trunc(testResult$lengths[2]), 36533)
  expect_identical(trunc(testResult$nsSNPs[2]), 1521)
  expect_equal(as.numeric(testResult$percent[2]), 0.0416)
  expect_equal(testResult$geneName[3], "MORN1")
  expect_identical(trunc(testResult$lengths[3]), 13485)
  expect_identical(trunc(testResult$nsSNPs[3]), 686)
  expect_equal(as.numeric(testResult$percent[3]), 0.0509)
  expect_equal(testResult$geneName[4], "MORN1")
  expect_identical(trunc(testResult$lengths[4]), 1993)
  expect_identical(trunc(testResult$nsSNPs[4]), 119)
  expect_equal(as.numeric(testResult$percent[4]), 0.0597)
})


test_that("invalid input1: wrong chromosome name", {
  chrName = 'a'
  startPosition = 2321253
  endPosition = 2391707

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    startPosition = startPosition,
    endPosition = endPosition))

  chrName = 25

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    startPosition = startPosition,
    endPosition = endPosition))
})


test_that("invalid input2: wrong position input", {
  chrName = 1
  startPosition = -2
  endPosition = 2391604

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    startPosition = startPosition,
    endPosition = endPosition))
})


test_that("invalid input3: no encoding protein", {
  chrName = "X"
  startPosition = 2794535
  endPosition = 2808377

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    startPosition = startPosition,
    endPosition = endPosition))
})

# [END]
