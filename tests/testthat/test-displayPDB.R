context("Display 3D protein strucutre for input gene and point the input
        nsSNP position.")
library(nsSNPfinder)

test_that("invalid input1: incorrect chromosome name", {
  chrName = 'r'
  geneName = 'RHOA'
  nsPos = 49395585

  expect_error(testresult <- displayPDB(chrName = chrName,
                                   geneName = geneName,
                                   nsPos = nsPos))
})

test_that("invalid input2: incorrect gene name", {
  chrName = 'X'
  geneName = 'afalsename'
  nsPos = 49395585

  expect_error(testresult <- displayPDB(chrName = chrName,
                                   geneName = geneName,
                                   nsPos = nsPos))
})


test_that("invalid input3: incorrect choice for pdb ID", {
  chrName = 3
  geneName = 'RHOA'
  nsPos = 49395585
  choice = 3

  expect_error(testresult <- displayPDB(chrName = chrName,
                                   geneName = geneName,
                                   nsPos = nsPos,
                                   choice = choice))
})


test_that("invalid input4: incorrect nsSNP position", {
  chrName = 3
  geneName = 'RHOA'
  nsPos = -200

  expect_error(testresult <- displayPDB(chrName = chrName,
                                   geneName = geneName,
                                   nsPos = nsPos,
                                   choice = choice))
})


test_that("invalid input5: no encoding protein in PDB", {
  chrName = 1
  geneName = 'MORN1'
  nsPos = 2321283

  expect_error(testresult <- displayPDB(chrName = chrName,
                                   geneName = geneName,
                                   nsPos = nsPos,
                                   choice = choice))
})

# [END]
