test_that("spmSnPM runs for t, F, and T2 and returns a valid permutation field", {
  set.seed(1)
  g1 <- matrix(rnorm(50 * 10), 50); g2 <- matrix(rnorm(50 * 10), 50)
  g2[20:30, ] <- g2[20:30, ] + 1.5
  st <- spmSnPM(cbind(g1, g2), group1 = 1:10, group2 = 11:20,
                n_permutations = 300, seed = 1)
  expect_s3_class(st, "spm_result")
  expect_equal(st$test_type, "snpm")
  expect_length(st$perm_max, 300)
  expect_true(all(st$p_values >= 0 & st$p_values <= 1))
  expect_gt(length(st$clusters), 0)

  sf <- spmSnPM(matrix(rnorm(50 * 24), 50), groups = rep(1:3, each = 8),
                statistic = "F", n_permutations = 200, seed = 2)
  expect_equal(sf$statistic, "F")
  expect_true(is.finite(sf$threshold))

  arr <- array(rnorm(50 * 16 * 3), c(50, 16, 3))
  arr[20:30, 1:8, ] <- arr[20:30, 1:8, ] + 1.2
  s2 <- spmSnPM(arr, groups = rep(c("A", "B"), each = 8), statistic = "T2",
                n_permutations = 200, seed = 3)
  expect_equal(s2$statistic, "T2")
  expect_gt(length(s2$clusters), 0)
})

test_that("the permutation threshold is the (1 - alpha) field-max quantile", {
  set.seed(11)
  data <- matrix(rnorm(40 * 20), 40)
  s <- spmSnPM(data, group1 = 1:10, group2 = 11:20, n_permutations = 500,
               alpha = 0.05, seed = 11)
  expect_equal(s$threshold,
               as.numeric(quantile(s$perm_max, 0.95, type = 1, names = FALSE)))
  # the observed field-max is the first permutation entry (identity)
  expect_equal(s$perm_max[1], max(abs(s$t)))
})

test_that("spmSnPM is reproducible with a seed", {
  data <- matrix(rnorm(30 * 16), 30)
  a <- spmSnPM(data, group1 = 1:8, group2 = 9:16, n_permutations = 200, seed = 5)
  b <- spmSnPM(data, group1 = 1:8, group2 = 9:16, n_permutations = 200, seed = 5)
  expect_equal(a$threshold, b$threshold)
  expect_equal(a$perm_max, b$perm_max)
})

test_that("the permutation threshold converges to the parametric RFT threshold", {
  # smooth null field: kernel-smoothed white noise (FWHM ~ 15 nodes)
  set.seed(3)
  fwhm <- 15; sd <- fwhm / sqrt(8 * log(2)); k <- ceiling(3 * sd)
  w <- exp(-(-k:k)^2 / (2 * sd^2)); w <- w / sqrt(sum(w^2))
  smooth1 <- function(n) as.numeric(stats::filter(rnorm(n + 2 * k), w))[(k + 1):(k + n)]
  data <- sapply(1:24, function(i) smooth1(120))
  sn <- spmSnPM(data, group1 = 1:12, group2 = 13:24, n_permutations = 1000,
                seed = 3)
  pr <- spmTTest(data, group1 = 1:12, group2 = 13:24)
  expect_lt(abs(sn$threshold - pr$threshold) / pr$threshold, 0.05)  # within 5%
})

test_that("spmSnPM handles n_permutations = 1 without corrupting the identity entry", {
  set.seed(2)
  data <- matrix(rnorm(20 * 8), 20)
  s <- spmSnPM(data, group1 = 1:4, group2 = 5:8, n_permutations = 1, seed = 2)
  expect_length(s$perm_max, 1L)                          # not 2 (the 2:1 antipattern)
  expect_equal(s$perm_max[1], max(abs(s$t)))             # still the observed max
  expect_error(spmSnPM(data, group1 = 1:4, group2 = 5:8, n_permutations = 0),
               ">= 1")
})

test_that("spmMANOVA rejects component matrices of mismatched dimensions", {
  bad <- list(matrix(rnorm(20 * 12), 20), matrix(rnorm(18 * 12), 18))
  expect_error(spmMANOVA(bad, groups = rep(c("A", "B"), each = 6)),
               "same dimensions")
})

test_that("spmSnPM validates the T2 group count", {
  arr <- array(rnorm(30 * 24 * 2), c(30, 24, 2))
  expect_error(spmSnPM(arr, groups = rep(1:3, each = 8), statistic = "T2",
                       n_permutations = 50), "two groups")
})
