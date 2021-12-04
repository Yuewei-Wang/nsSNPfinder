context("Calculate nsSNP percentage within the gene range")
library(nsSNPfinder)

test_that("valid input: approper arguments in chromosome 1", {
  chrName = 3
  startPosition = 49359145
  endPosition = 49411645

  testResult <- nsSNPCalculatebyRange(chrName = chrName,
                                      startPosition = startPosition,
                                      endPosition = endPosition)

  expect_type(testResult, "list")
  expect_length(nrow(testResult), 1)
  expect_s3_class(testResult, "data.frame")
  expect_equal(testResult$geneName[1], "RHOA")
  expect_identical(trunc(testResult$lengths[1]), 52501)
  expect_identical(trunc(testResult$nsSNPs[1]), 1614)
  expect_equal(as.numeric(testResult$percent[1]), 0.0307)
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
  startPosition = 2816501
  endPosition = 2816600

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    startPosition = startPosition,
    endPosition = endPosition))
})

# [END]
