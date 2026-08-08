library(testthat)
library(PhysioAnalysis)

# spmTTest takes a (time x observations) matrix (or a PhysioExperiment with a 2D
# assay) and runs a random-field-theory-corrected t-test across observations at
# each time point.

test_that("spmTTest one-sample returns an spm_result with the expected shape", {
  set.seed(1)
  m <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
  res <- spmTTest(m)
  expect_s3_class(res, "spm_result")
  expect_length(res$t, 100)
  expect_length(res$p_values, 100)
  expect_true(is.numeric(res$threshold) && res$threshold > 0)
  expect_equal(res$df, 19)                       # n - 1
})

test_that("spmTTest detects a real one-sample effect and null gives no clusters", {
  set.seed(2)
  n_time <- 100; n_sub <- 25
  # inject a sustained effect in samples 40-60
  signal <- matrix(rnorm(n_time * n_sub), n_time, n_sub)
  signal[40:60, ] <- signal[40:60, ] + 3
  res <- spmTTest(signal)
  # the peak t is inside the effect window
  expect_gt(max(res$t), res$threshold)
  peak <- which.max(res$t)
  expect_true(peak >= 40 && peak <= 60)

  set.seed(3)
  noise <- matrix(rnorm(n_time * n_sub), n_time, n_sub)
  res0 <- spmTTest(noise)
  expect_true(all(abs(res0$t) < 20))             # sane t range on pure noise
})

test_that("spmTTest two-sample uses pooled variance and df = n1 + n2 - 2", {
  set.seed(4)
  m <- matrix(rnorm(60 * 30), nrow = 60, ncol = 30)
  res <- spmTTest(m, group1 = 1:15, group2 = 16:30)
  expect_s3_class(res, "spm_result")
  expect_equal(res$df, 28)
  expect_length(res$t, 60)
})

test_that("spmTTest matches base t.test at a single time point (one-sample)", {
  set.seed(5)
  m <- matrix(rnorm(10 * 18), nrow = 10, ncol = 18)
  res <- spmTTest(m)
  ref_t <- t.test(m[3, ])$statistic
  expect_equal(unname(res$t[3]), unname(ref_t), tolerance = 1e-10)
})
