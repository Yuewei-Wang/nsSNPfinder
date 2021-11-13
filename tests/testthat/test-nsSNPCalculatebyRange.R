context("Calculate nsSNP percentage within the gene range")
library(nsSNPfinder)

test_that("valid input: approper arguments in chromosome 1", {
  chrName = 1
  start_position = 2321253
  end_position = 2391707

  testResult <- nsSNPCalculatebyRange(chrName = chrName,
                                      start_position = start_position,
                                      end_position = end_position)

  expect_type(testResult, "list")
  expect_length(testResult, 4)
  expect_s3_class(testResult, "data.frame")
  expect_equal(testResult$gene_Name[1], "MORN1")
  expect_identical(trunc(testResult$lengths[1]), 70302)
  expect_identical(trunc(testResult$all_nssnps[1]), 3072)
  expect_equal(as.numeric(testResult$percent[1]), 0.0437)
  expect_equal(testResult$gene_Name[2], "MORN1")
  expect_identical(trunc(testResult$lengths[2]), 36533)
  expect_identical(trunc(testResult$all_nssnps[2]), 1521)
  expect_equal(as.numeric(testResult$percent[2]), 0.0416)
  expect_equal(testResult$gene_Name[3], "MORN1")
  expect_identical(trunc(testResult$lengths[3]), 13485)
  expect_identical(trunc(testResult$all_nssnps[3]), 686)
  expect_equal(as.numeric(testResult$percent[3]), 0.0509)
  expect_equal(testResult$gene_Name[4], "MORN1")
  expect_identical(trunc(testResult$lengths[4]), 1993)
  expect_identical(trunc(testResult$all_nssnps[4]), 119)
  expect_equal(as.numeric(testResult$percent[4]), 0.0597)
})


test_that("invalid input1: wrong chromosome name", {
  chrName = 'a'
  start_position = 2321253
  end_position = 2391707

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    start_position = start_position,
    end_position = end_position))

  chrName = 25

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    start_position = start_position,
    end_position = end_position))
})


test_that("invalid input2: wrong position input", {
  chrName = 1
  start_position = -2
  end_position = 2391604

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    start_position = start_position,
    end_position = end_position))
})


test_that("invalid input3: no encoding protein", {
  chrName = "X"
  start_position = 2794535
  end_position = 2808377

  expect_error(testResult <- nsSNPCalculatebyRange(
    chrName = chrName,
    start_position = start_position,
    end_position = end_position))
})

# [END]
