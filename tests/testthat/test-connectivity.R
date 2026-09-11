# Tests for connectivity analysis functions

test_that("coherence computes channel-pair coherence", {
  set.seed(123)
  # Create data with two correlated channels
  n_samples <- 500
  sr <- 100
  t <- seq(0, n_samples/sr, length.out = n_samples)

  # Create signals with shared oscillation
  shared <- sin(2 * pi * 10 * t)  # 10 Hz shared component
  ch1 <- shared + 0.5 * rnorm(n_samples)
  ch2 <- shared + 0.5 * rnorm(n_samples)
  ch3 <- rnorm(n_samples)  # Uncorrelated channel

  data <- cbind(ch1, ch2, ch3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = c("Ch1", "Ch2", "Ch3")),
    samplingRate = sr
  )

  result <- coherence(pe)

  expect_type(result, "list")
  expect_true("coherence" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("channel_names" %in% names(result))

  # Coherence should be 3D array: freq x channel x channel
  expect_equal(length(dim(result$coherence)), 3)
  expect_equal(dim(result$coherence)[2], 3)  # 3 channels
  expect_equal(dim(result$coherence)[3], 3)  # 3 channels
})

test_that("coherence works with frequency range filter", {
  set.seed(123)
  data <- matrix(rnorm(500 * 4), nrow = 500, ncol = 4)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- coherence(pe, freq_range = c(8, 12))

  expect_true(all(result$frequencies >= 8))
  expect_true(all(result$frequencies <= 12))
})

test_that("coherence works with specific channels", {
  set.seed(123)
  data <- matrix(rnorm(500 * 5), nrow = 500, ncol = 5)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = paste0("Ch", 1:5)),
    samplingRate = 100
  )

  result <- coherence(pe, channels = c(1, 2, 3))

  # Should have 3 channels
  expect_equal(dim(result$coherence)[2], 3)
  expect_equal(dim(result$coherence)[3], 3)
})

test_that("crossSpectrum computes cross-spectral density", {
  set.seed(123)
  data <- matrix(rnorm(500 * 3), nrow = 500, ncol = 3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- crossSpectrum(pe)

  expect_type(result, "list")
  expect_true("csd" %in% names(result))
  expect_true("frequencies" %in% names(result))
  expect_true("channel_names" %in% names(result))
})

test_that("plv computes phase locking value", {
  set.seed(123)
  # Create phase-locked signals
  n_samples <- 1000
  sr <- 256
  t <- seq(0, n_samples/sr, length.out = n_samples)

  # Two channels with similar phase in alpha band
  phase <- 2 * pi * 10 * t
  ch1 <- sin(phase) + 0.3 * rnorm(n_samples)
  ch2 <- sin(phase + 0.1) + 0.3 * rnorm(n_samples)  # Small phase difference
  ch3 <- sin(phase + runif(n_samples, 0, 2*pi))  # Random phase

  data <- cbind(ch1, ch2, ch3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = c("Ch1", "Ch2", "Ch3")),
    samplingRate = sr
  )

  result <- plv(pe, freq_band = c(8, 12))

  # plv returns a matrix (channel x channel)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))

  # PLV values should be between 0 and 1
  expect_true(all(result >= 0 & result <= 1))

  # Diagonal should be 1
  expect_equal(unname(diag(result)), c(1, 1, 1))
})

test_that("pli computes phase lag index", {
  set.seed(123)
  data <- matrix(rnorm(500 * 3), nrow = 500, ncol = 3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- pli(pe, freq_band = c(8, 12))

  # pli returns a matrix (channel x channel)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))

  # PLI values should be between 0 and 1
  expect_true(all(result >= 0 & result <= 1))
})

test_that("wPLI computes weighted phase lag index", {
  set.seed(123)
  data <- matrix(rnorm(500 * 3), nrow = 500, ncol = 3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- wPLI(pe, freq_band = c(8, 12))

  # wPLI returns a matrix (channel x channel)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))

  # wPLI values should be between 0 and 1
  expect_true(all(result >= 0 & result <= 1))
})

test_that("correlationMatrix computes pairwise correlations", {
  set.seed(123)
  # Create data with known correlation structure
  n_samples <- 200
  base <- rnorm(n_samples)
  ch1 <- base + 0.1 * rnorm(n_samples)
  ch2 <- base + 0.1 * rnorm(n_samples)
  ch3 <- rnorm(n_samples)

  data <- cbind(ch1, ch2, ch3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = c("Ch1", "Ch2", "Ch3")),
    samplingRate = 100
  )

  result <- correlationMatrix(pe)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))

  # Diagonal should be 1
  expect_equal(unname(diag(result)), c(1, 1, 1))

  # Matrix should be symmetric
  expect_equal(result, t(result))

  # Ch1 and Ch2 should be highly correlated
  expect_true(result[1, 2] > 0.8)
})

test_that("correlationMatrix works with different methods", {
  set.seed(123)
  data <- matrix(rnorm(200 * 4), nrow = 200, ncol = 4)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  pearson <- correlationMatrix(pe, method = "pearson")
  spearman <- correlationMatrix(pe, method = "spearman")

  expect_true(is.matrix(pearson))
  expect_true(is.matrix(spearman))

  # Results should differ
  expect_false(identical(pearson, spearman))
})

test_that("connectivityMatrix works with coherence method", {
  set.seed(123)
  data <- matrix(rnorm(500 * 4), nrow = 500, ncol = 4)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- connectivityMatrix(pe, method = "coherence", freq_band = c(8, 12))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 4))

  # Matrix should be symmetric
  expect_equal(result, t(result))

  # Diagonal should be 1 (coherence with self)
  expect_equal(diag(result), rep(1, 4))
})

test_that("connectivityMatrix works with plv method", {
  set.seed(123)
  data <- matrix(rnorm(500 * 3), nrow = 500, ncol = 3)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- connectivityMatrix(pe, method = "plv", freq_band = c(8, 12))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
})

test_that("connectivityMatrix works with correlation method", {
  set.seed(123)
  data <- matrix(rnorm(200 * 4), nrow = 200, ncol = 4)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  result <- connectivityMatrix(pe, method = "correlation")

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 4))
})

test_that("connectivity functions error on insufficient data", {
  data <- matrix(rnorm(10 * 2), nrow = 10, ncol = 2)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    samplingRate = 100
  )

  # Should work but may warn about short signals
  expect_type(correlationMatrix(pe), "double")
})

test_that("connectivity functions work with epoched data", {
  set.seed(123)
  epochs <- array(rnorm(100 * 4 * 10), dim = c(100, 4, 10))
  pe <- PhysioExperiment(
    assays = list(epoched = epochs),
    samplingRate = 100
  )

  # For epoched data, correlationMatrix should work on first epoch
  result <- correlationMatrix(pe)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 4))
})
