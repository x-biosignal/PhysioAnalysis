# Golden regression tests (WSF-08).
#
# Each test rebuilds the SAME deterministic input used in data-raw/golden.R,
# runs the PACKAGE function, and compares against a golden captured from an
# INDEPENDENT reference (igraph / base R / analytic one-sided FFT PSD).
# Goldens live in tests/testthat/_golden/ and skip cleanly when absent.

library(testthat)
library(PhysioAnalysis)

# --- 1. Network: clustering coefficient vs igraph ---------------------------
test_that("clusteringCoefficient matches igraph transitivity (golden)", {
  set.seed(101)
  n_net <- 12
  A <- matrix(0, n_net, n_net)
  ut <- upper.tri(A)
  A[ut] <- rbinom(sum(ut), 1, 0.35)
  A <- A + t(A)
  diag(A) <- 0

  cc <- clusteringCoefficient(A)
  expect_equal_golden(cc, "net_clustering_igraph", tol = 1e-8)
})

# --- 1b. Network: betweenness centrality vs igraph --------------------------
test_that("betweennessCentrality matches igraph betweenness (golden)", {
  set.seed(101)
  n_net <- 12
  A <- matrix(0, n_net, n_net)
  ut <- upper.tri(A)
  A[ut] <- rbinom(sum(ut), 1, 0.35)
  A <- A + t(A)
  diag(A) <- 0

  bc <- betweennessCentrality(A, normalized = TRUE)
  expect_equal_golden(bc, "net_betweenness_igraph", tol = 1e-8)
})

# --- 2. Correlation connectivity vs base cor() ------------------------------
test_that("correlationMatrix (pearson) matches base cor (golden)", {
  set.seed(202)
  Xc <- matrix(rnorm(500 * 5), 500, 5)
  colnames(Xc) <- paste0("Ch", 1:5)
  pe <- PhysioExperiment(
    assays = list(raw = Xc),
    colData = S4Vectors::DataFrame(label = paste0("Ch", 1:5)),
    samplingRate = 100
  )
  cm <- correlationMatrix(pe, method = "pearson")
  expect_equal_golden(cm, "conn_correlation_pearson_base", tol = 1e-10)
})

test_that("correlationMatrix (spearman) matches base cor (golden)", {
  set.seed(202)
  Xc <- matrix(rnorm(500 * 5), 500, 5)
  colnames(Xc) <- paste0("Ch", 1:5)
  pe <- PhysioExperiment(
    assays = list(raw = Xc),
    colData = S4Vectors::DataFrame(label = paste0("Ch", 1:5)),
    samplingRate = 100
  )
  cm <- correlationMatrix(pe, method = "spearman")
  expect_equal_golden(cm, "conn_correlation_spearman_base", tol = 1e-10)
})

# --- 3. Spectrogram vs hand-computed one-sided FFT PSD ----------------------
test_that("spectrogram (single window) matches hand FFT PSD (golden)", {
  set.seed(404)
  sig_s <- rnorm(256)
  Xs <- matrix(sig_s, nrow = 256, ncol = 1)
  pe <- PhysioExperiment(assays = list(raw = Xs), samplingRate = 256)

  spec <- spectrogram(pe, window_size = 256L, overlap = 0.5,
                      window_type = "hanning", channel = 1L)
  # Single 256-sample window -> one column.
  expect_equal(ncol(spec$power), 1L)
  expect_equal_golden(spec$power[, 1], "spectrogram_onesided_psd_handfft",
                      tol = 1e-9)
})

# --- 4. bandPower (Welch) vs independent from-scratch Welch PSD -------------
test_that("bandPower (welch) matches independent Welch PSD (golden)", {
  set.seed(505)
  sig_b <- rnorm(512)
  Xb <- matrix(sig_b, nrow = 512, ncol = 1)
  pe <- PhysioExperiment(
    assays = list(raw = Xb),
    colData = S4Vectors::DataFrame(label = "Ch1"),
    samplingRate = 256
  )
  bands <- list(theta = c(4, 8), alpha = c(8, 13), beta = c(13, 30))
  bp <- bandPower(pe, bands = bands, method = "welch")

  actual <- c(theta = bp$theta[1], alpha = bp$alpha[1], beta = bp$beta[1])
  expect_equal_golden(actual, "bandpower_welch_handfft", tol = 1e-9)
})

# --- 5. SPM t-map vs base t.test per time point -----------------------------
test_that("spmTTest two-sample t-map matches base t.test (golden)", {
  set.seed(303)
  nt <- 60L; n1 <- 8L; n2 <- 9L
  g1 <- matrix(rnorm(nt * n1), nt, n1)
  g2 <- matrix(rnorm(nt * n2), nt, n2) + 0.3
  dat_two <- cbind(g1, g2)
  pe <- PhysioExperiment(assays = list(values = dat_two), samplingRate = 100)

  res <- spmTTest(pe, group1 = 1:8, group2 = 9:17, two_tailed = TRUE)
  expect_equal_golden(as.numeric(res$t), "spm_ttest_twosample_base", tol = 1e-8)
})

test_that("spmTTest one-sample t-map matches base t.test (golden)", {
  set.seed(303)
  nt <- 60L; n1 <- 8L; n2 <- 9L
  g1 <- matrix(rnorm(nt * n1), nt, n1)
  g2 <- matrix(rnorm(nt * n2), nt, n2) + 0.3
  dat_two <- cbind(g1, g2)

  res <- spmTTest(dat_two[, 1:8], two_tailed = TRUE)
  expect_equal_golden(as.numeric(res$t), "spm_ttest_onesample_base", tol = 1e-8)
})
