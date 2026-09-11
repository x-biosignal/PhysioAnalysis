# Functional regression, verified against known coefficient functions.

test_that("function-on-scalar recovers a known coefficient curve", {
  set.seed(1)
  P <- 101; t <- seq(0, 1, length.out = P); N <- 60
  beta_shape <- sin(2 * pi * t)                       # true effect of speed on the curve
  speed <- runif(N, -1, 1)
  curves <- outer(speed, beta_shape) + matrix(rnorm(N * P, 0, 0.1), N, P)
  fit <- functionalRegression(curves, speed, n_perm = 300, seed = 42)
  expect_s3_class(fit, "fosr_result")
  # slope coefficient curve (row 2) correlates with the true shape
  expect_gt(cor(fit$coefficients[2, ], beta_shape), 0.95)
  expect_lt(fit$p_global[["x1"]], 0.05)              # speed effect significant
  expect_equal(dim(fit$coefficients), c(2, P))
})

test_that("function-on-scalar reports a null effect as non-significant", {
  set.seed(2)
  P <- 60; N <- 50
  curves <- matrix(rnorm(N * P), N, P)               # no predictor relationship
  z <- rnorm(N)
  fit <- functionalRegression(curves, z, n_perm = 500, seed = 7)
  expect_gt(fit$p_global[["x1"]], 0.05)
})

test_that("function-on-scalar handles multiple predictors", {
  set.seed(3)
  P <- 40; t <- seq(0, 1, length.out = P); N <- 80
  x1 <- rnorm(N); x2 <- rnorm(N)
  curves <- outer(x1, sin(2 * pi * t)) + outer(x2, cos(2 * pi * t)) +
    matrix(rnorm(N * P, 0, 0.1), N, P)
  fit <- functionalRegression(curves, cbind(x1 = x1, x2 = x2), n_perm = 300, seed = 1)
  expect_equal(nrow(fit$coefficients), 3)            # intercept + 2
  expect_lt(fit$p_global[["x1"]], 0.05)
  expect_lt(fit$p_global[["x2"]], 0.05)
  expect_output(print(fit), "Function-on-scalar")
})

test_that("scalar-on-function predicts from a curve and localises the effect", {
  set.seed(4)
  P <- 80; t <- seq(0, 1, length.out = P); N <- 120
  X <- matrix(rnorm(N * P), N, P) +
    outer(rnorm(N), sin(2 * pi * t))                 # add smooth structure
  beta_true <- dnorm(t, 0.3, 0.05); beta_true <- beta_true / max(beta_true)
  y <- as.numeric(X %*% beta_true) + rnorm(N, 0, 0.5)
  m <- scalarOnFunctionRegression(y, X, npc = 8)
  expect_s3_class(m, "sofr_result")
  expect_gt(m$r_squared, 0.6)
  # recovered beta(t) peaks near the true peak (t = 0.3)
  expect_lt(abs(t[which.max(m$beta)] - 0.3), 0.12)
  expect_length(m$beta, P)
})

test_that("scalar-on-function auto-selects npc by variance explained", {
  set.seed(5)
  P <- 50; N <- 60
  X <- matrix(rnorm(N * P), N, P)
  y <- rnorm(N)
  m <- scalarOnFunctionRegression(y, X, ve = 0.9)
  expect_true(m$npc >= 1 && m$npc <= min(N, P))
  expect_gte(m$pve, 0.9 - 1e-8)
})
