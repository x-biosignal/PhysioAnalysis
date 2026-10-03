test_that("spmRegression node-wise t equals the lm() slope t exactly", {
  set.seed(42)
  n_time <- 60; n_obs <- 18
  pred <- rnorm(n_obs)
  data <- matrix(rnorm(n_time * n_obs), nrow = n_time)
  data[20:35, ] <- data[20:35, ] + outer(rep(1.2, 16), pred)
  r <- spmRegression(data, pred)
  t_lm <- vapply(seq_len(n_time), function(i)
    summary(lm(data[i, ] ~ pred))$coefficients["pred", "t value"], numeric(1))
  expect_equal(r$t, t_lm, tolerance = 1e-10)
  expect_equal(r$df, n_obs - 2L)
  # slope + intercept fields match lm()
  expect_equal(unname(r$slope[25]), unname(coef(lm(data[25, ] ~ pred))["pred"]),
               tolerance = 1e-10)
  # uncorrected pointwise p at the peak matches lm()
  pk <- which.max(abs(r$t))
  expect_equal(r$p_values[pk],
               summary(lm(data[pk, ] ~ pred))$coefficients["pred", "Pr(>|t|)"],
               tolerance = 1e-10)
})

test_that("spmRegression detects an injected slope window and is class spm_result", {
  set.seed(7)
  pred <- rnorm(20)
  data <- matrix(rnorm(80 * 20), nrow = 80)
  data[30:50, ] <- data[30:50, ] + outer(rep(2, 21), pred)  # strong slope
  r <- spmRegression(data, pred)
  expect_s3_class(r, "spm_result")
  expect_equal(r$test_type, "regression")
  expect_gt(length(r$clusters), 0)
  # a detected cluster overlaps the injected window
  starts <- vapply(r$clusters, `[[`, numeric(1), "start")
  ends <- vapply(r$clusters, `[[`, numeric(1), "end")
  expect_true(any(starts <= 50 & ends >= 30))
})

test_that("spmRegression validates its input", {
  data <- matrix(rnorm(10 * 6), nrow = 10)
  expect_error(spmRegression(data, rnorm(5)), "one value per observation")
  expect_error(spmRegression(data, rep(1, 6)), "zero variance")
  expect_error(spmRegression(matrix(rnorm(10 * 2), 10), c(1, 2)),
               "at least 3 observations")
  expect_error(spmRegression("not a matrix", 1:3), "PhysioExperiment or matrix")
})
