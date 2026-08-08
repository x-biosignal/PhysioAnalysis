library(testthat)
library(PhysioAnalysis)

# internal accessors
rc_t <- function(a, df, rc, tt = TRUE) PhysioAnalysis:::.rftCriticalT(a, df, rc, tt)
rc_f <- function(a, d1, d2, rc) PhysioAnalysis:::.rftCriticalF(a, d1, d2, rc)
ecT <- function(u, df) PhysioAnalysis:::.ecDensityT(u, df)

# Gaussian-smoothed columns (edge-padded), for smooth random fields
.sm <- function(m, k = 9) {
  w <- dnorm(seq(-3, 3, length.out = k)); w <- w / sum(w)
  apply(m, 2, function(col) {
    padded <- c(rep(col[1], k), col, rep(col[length(col)], k))
    as.numeric(stats::filter(padded, w, sides = 2))[(k + 1):(k + length(col))]
  })
}

test_that("RFT critical-t solves the expected-EC equation at the target level", {
  for (df in c(6, 10, 20, 40)) for (rc in c(2, 5, 10, 20)) {
    u <- rc_t(0.05, df, rc, TRUE)
    ec <- ecT(u, df)
    expect_equal(unname(ec["rho0"] + rc * ec["rho1"]), 0.025, tolerance = 1e-6)
    expect_gt(u, stats::qt(0.975, df))              # above single comparison
  }
})

test_that("RFT thresholds are monotone in resel count (up) and df (down)", {
  u_rc <- vapply(c(1, 5, 10, 20, 40), function(rc) rc_t(0.05, 10, rc), numeric(1))
  expect_true(all(diff(u_rc) > 0))
  u_df <- vapply(c(6, 10, 20, 40, 100), function(df) rc_t(0.05, df, 10), numeric(1))
  expect_true(all(diff(u_df) < 0))
  uf <- vapply(c(2, 5, 10, 20), function(rc) rc_f(0.05, 3, 40, rc), numeric(1))
  expect_true(all(diff(uf) > 0))
})

test_that("as smoothness grows the threshold approaches the single comparison", {
  expect_equal(rc_t(0.05, 20, 1e-6, TRUE), stats::qt(0.975, 20), tolerance = 1e-3)
  expect_equal(rc_f(0.05, 3, 40, 1e-6), stats::qf(0.95, 3, 40), tolerance = 1e-3)
  expect_lt(rc_t(0.05, 20, (101 - 1) / 80), stats::qt(0.975, 20) + 0.6)
})

test_that("critical thresholds reproduce the stored RFT reference fixture", {
  ref <- readRDS(test_path("fixtures", "rft-reference.rds"))
  t_now <- mapply(function(a, df, rc) rc_t(a, df, rc, TRUE),
                  ref$t$alpha, ref$t$df, ref$t$resel)
  expect_equal(t_now, ref$t$crit_t, tolerance = 1e-6)
  f_now <- mapply(function(a, d1, d2, rc) rc_f(a, d1, d2, rc),
                  ref$f$alpha, ref$f$df1, ref$f$df2, ref$f$resel)
  expect_equal(f_now, ref$f$crit_f, tolerance = 1e-6)
})

test_that("FWHM is estimated from residuals, differing from the statistic value", {
  # With a real effect the t-statistic field has structure (a sharp supra-
  # threshold bump) that does NOT reflect the residual smoothness, so the old
  # statistic-based estimator and the correct residual-based estimator diverge.
  set.seed(11)
  n_time <- 101; n <- 12
  g1 <- .sm(matrix(rnorm(n_time * n), n_time, n))
  g2 <- .sm(matrix(rnorm(n_time * n), n_time, n)); g2[40:60, ] <- g2[40:60, ] + 1.8
  res <- spmTTest(cbind(g1, g2), group1 = 1:n, group2 = (n + 1):(2 * n))
  fwhm_stat <- PhysioAnalysis:::.estimateFWHM(res$t)   # legacy: from the statistic
  expect_gt(abs(res$fwhm - fwhm_stat), 0.5)            # guards the specific bug
  expect_true(res$fwhm > 1 && res$fwhm < n_time)
})

test_that("a smoothed two-sample effect is recovered as the correct cluster", {
  set.seed(21)
  n_time <- 101; n <- 12
  g1 <- .sm(matrix(rnorm(n_time * n), n_time, n))
  g2 <- .sm(matrix(rnorm(n_time * n), n_time, n)); g2[40:60, ] <- g2[40:60, ] + 1.6
  res <- spmTTest(cbind(g1, g2), group1 = 1:n, group2 = (n + 1):(2 * n))
  expect_gte(length(res$clusters), 1)
  cl <- res$clusters[[which.max(vapply(res$clusters, `[[`, numeric(1), "extent"))]]
  expect_true(cl$start >= 30 && cl$end <= 72)       # near the effect
  expect_true(cl$start <= 45 && cl$end >= 55)       # covers the effect core
  expect_equal(cl$direction, "negative")            # g1 - g2 < 0
  expect_lt(cl$p_cluster, 0.05)
  set.seed(22)
  n0 <- cbind(.sm(matrix(rnorm(n_time * n), n_time, n)),
              .sm(matrix(rnorm(n_time * n), n_time, n)))
  expect_equal(length(spmTTest(n0, group1 = 1:n, group2 = (n + 1):(2 * n))$clusters), 0)
})

test_that("ANOVA uses the residual FWHM and the F critical threshold", {
  set.seed(31)
  n_time <- 101; per <- 8
  d <- cbind(.sm(matrix(rnorm(n_time * per), n_time, per)),
             .sm(matrix(rnorm(n_time * per), n_time, per)),
             .sm(matrix(rnorm(n_time * per), n_time, per)))
  d[40:60, (2 * per + 1):(3 * per)] <- d[40:60, (2 * per + 1):(3 * per)] + 1.9
  grp <- factor(rep(c("A", "B", "C"), each = per))
  res <- spmAnova(d, grp)
  expect_equal(res$df1, 2); expect_equal(res$df2, 3 * per - 3)
  expect_true(res$fwhm > 1 && res$fwhm < n_time)
  expect_gt(res$threshold, stats::qf(0.95, res$df1, res$df2))  # RFT-corrected
  expect_gte(length(res$clusters), 1)
})

# --- regression tests for adversarial-review findings (WS8-01) -----------------

test_that("degenerate low-df fields yield an infinite threshold, not a crash", {
  # df = 1 t-field: the EC density does not decay, so no finite RFT threshold
  expect_equal(rc_t(0.05, 1, 10, TRUE), Inf)
  expect_true(is.finite(rc_t(0.05, 1, 0.05, TRUE)))   # tiny resel -> root exists
  # F-field: only df2 = 1 has a non-decaying EC density; df2 > 1 is finite
  expect_equal(rc_f(0.05, 3, 1, 10), Inf)
  expect_true(is.finite(rc_f(0.05, 2, 2, 10)))
  # a 3-observation two-sample test (df = 1) runs end-to-end without error
  set.seed(1)
  d3 <- matrix(rnorm(50 * 3), 50, 3)
  res <- spmTTest(d3, group1 = 1:2, group2 = 3)
  expect_s3_class(res, "spm_result")
  expect_equal(length(res$clusters), 0)               # Inf threshold -> none
})

test_that("spmAnova handles unused factor levels and rejects df2 < 1", {
  set.seed(2)
  d <- matrix(rnorm(50 * 12), 50, 12)
  grp <- factor(rep(c("A", "B", "C"), each = 4), levels = c("A", "B", "C", "D"))
  res <- spmAnova(d, grp)                              # level D is empty
  expect_equal(res$df1, 2L)                            # dropped to 3 groups
  expect_equal(res$df2, 12L - 3L)
  # one observation per group -> residual df 0 -> clean error, not a NaN crash
  expect_error(spmAnova(matrix(rnorm(50 * 3), 50, 3), factor(c("A", "B", "C"))),
               "residual df")
})

test_that("print.spm_result reports the resel count", {
  res <- spmTTest(matrix(rnorm(50 * 20), 50, 20), group1 = 1:10, group2 = 11:20)
  expect_output(print(res), "Resel count")
})

test_that("the 1D F-field EC density obeys the F(1, df2) = t(df2)^2 identity", {
  # rho1 of an F(1, df2) field equals twice the one-sided t(df2) density
  for (df2 in c(10, 20, 40)) for (u in c(2, 4, 8, 15)) {
    fF <- PhysioAnalysis:::.ecDensityF(u, 1, df2)["rho1"]
    fT <- 2 * PhysioAnalysis:::.ecDensityT(sqrt(u), df2)["rho1"]
    expect_equal(unname(fF), unname(fT), tolerance = 1e-9)
  }
  # and the RFT critical F threshold equals the squared two-tailed t threshold
  for (df2 in c(10, 20, 40)) for (rc in c(5, 10, 20)) {
    uF <- rc_f(0.05, 1, df2, rc)
    uT <- rc_t(0.05, df2, rc, TRUE)
    expect_equal(uF, uT^2, tolerance = 1e-4)
  }
  # the density is non-negative on the search range (no spurious 2D bracket)
  expect_gt(PhysioAnalysis:::.ecDensityF(0.05, 3, 16)["rho1"], 0)
})
