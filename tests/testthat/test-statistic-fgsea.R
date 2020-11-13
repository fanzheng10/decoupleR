library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "fgsea", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
  readRDS()

# test for run_fgsea() ----------------------------------------------------

test_that("test run_fgsea with dorothea gene sets", {
  res_1 <- run_fgsea(emat, dorothea_genesets)
  exp_1 <- file.path(expected_dir, "output-fgsea_dorothea_default.rds") %>%
    readRDS()

  res_2 <- run_fgsea(emat, dorothea_genesets, tf, target)
  exp_2 <- file.path(expected_dir, "output-fgsea_dorothea_tidy-evaluation.rds") %>%
    readRDS()

  res_3 <- run_fgsea(emat, dorothea_genesets, options = list(minSize = 4))
  exp_3 <- file.path(expected_dir, "output-fgsea_dorothea_minsize4.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1)
  expect_equal(res_2, exp_2)
  expect_equal(res_3, exp_3)
})
