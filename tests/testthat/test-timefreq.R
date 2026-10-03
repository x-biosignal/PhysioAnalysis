library(testthat)
library(PhysioAnalysis)

make_test_signal <- function(sr = 500, duration = 2) {
  n <- as.integer(sr * duration)
  t <- seq(0, duration - 1/sr, length.out = n)

  # Create signal with known frequency components
  # 10 Hz and 30 Hz sine waves
  signal <- sin(2 * pi * 10 * t) + 0.5 * sin(2 * pi * 30 * t)

  data <- matrix(signal, ncol = 1)
  # rowData must have n rows (matching dim[1] = time points)
  row_data <- S4Vectors::DataFrame(time_idx = seq_len(n))
  # colData must have 1 row (matching dim[2] = channels)
  col_data <- S4Vectors::DataFrame(label = "Ch1", type = "EEG")

  # Convert to 3D array (time x channels x samples)
  assays <- S4Vectors::SimpleList(raw = array(data, dim = c(n, 1, 1)))
  PhysioExperiment(assays, rowData = row_data, colData = col_data, samplingRate = sr)
}

test_that("spectrogram computes STFT", {
  x <- make_test_signal(sr = 500, duration = 2)

  spec <- spectrogram(x, window_size = 256, overlap = 0.5)

  expect_type(spec, "list")
  expect_true("power" %in% names(spec))
  expect_true("frequencies" %in% names(spec))
  expect_true("times" %in% names(spec))

  # Check dimensions
  expect_equal(nrow(spec$power), floor(256 / 2) + 1)
  expect_true(length(spec$times) > 0)
  expect_equal(length(spec$frequencies), nrow(spec$power))
})

test_that("spectrogram identifies frequency peaks", {
  x <- make_test_signal(sr = 500, duration = 2)

  spec <- spectrogram(x, window_size = 256, overlap = 0.5)

  # Find frequency with maximum average power
  avg_power <- rowMeans(spec$power)

  # 10 Hz should have significant power
  idx_10hz <- which.min(abs(spec$frequencies - 10))
  idx_30hz <- which.min(abs(spec$frequencies - 30))

  # Power at 10 Hz should be higher than at 50 Hz (no signal there)
  idx_50hz <- which.min(abs(spec$frequencies - 50))
  expect_true(avg_power[idx_10hz] > avg_power[idx_50hz])
})

test_that("spectrogram window types work", {
  x <- make_test_signal()

  for (wtype in c("hanning", "hamming", "blackman", "rectangular")) {
    spec <- spectrogram(x, window_size = 128, window_type = wtype)
    expect_type(spec, "list")
    expect_true(all(!is.na(spec$power)))
  }
})

test_that("waveletTransform computes CWT", {
  x <- make_test_signal(sr = 500, duration = 2)

  wt <- waveletTransform(x, frequencies = seq(5, 40, by = 5), n_cycles = 5)

  expect_type(wt, "list")
  expect_true("power" %in% names(wt))
  expect_true("phase" %in% names(wt))
  expect_true("frequencies" %in% names(wt))
  expect_true("times" %in% names(wt))

  # Check dimensions
  expect_equal(nrow(wt$power), length(seq(5, 40, by = 5)))
  expect_equal(ncol(wt$power), 1000)  # sr * duration
})

test_that("waveletTransform identifies frequency components", {
  x <- make_test_signal(sr = 500, duration = 2)

  wt <- waveletTransform(x, frequencies = seq(5, 50, by = 5))

  # Average power at each frequency
  avg_power <- rowMeans(wt$power)

  # 10 Hz should have peak power
  idx_10hz <- which(wt$frequencies == 10)
  idx_30hz <- which(wt$frequencies == 30)
  idx_45hz <- which(wt$frequencies == 45)

  expect_true(avg_power[idx_10hz] > avg_power[idx_45hz])
})

test_that("bandPower extracts power in frequency bands", {
  x <- make_test_signal(sr = 500, duration = 2)

  bp <- bandPower(x, method = "welch")

  expect_s3_class(bp, "data.frame")
  expect_true("channel" %in% names(bp))
  expect_true("alpha" %in% names(bp))  # 8-13 Hz
  expect_true("beta" %in% names(bp))   # 13-30 Hz

  # Our 10 Hz signal should show up in alpha band
  expect_true(bp$alpha[1] > 0)
})

test_that("bandPower with custom bands", {
  x <- make_test_signal(sr = 500, duration = 2)

  custom_bands <- list(
    low = c(5, 15),
    high = c(25, 35)
  )

  bp <- bandPower(x, bands = custom_bands, method = "welch")

  expect_true("low" %in% names(bp))
  expect_true("high" %in% names(bp))

  # 10 Hz signal in "low" band, 30 Hz in "high" band
  expect_true(bp$low[1] > 0)
  expect_true(bp$high[1] > 0)
})

test_that("bandPower relative option works", {
  x <- make_test_signal()

  bp_abs <- bandPower(x, relative = FALSE)
  bp_rel <- bandPower(x, relative = TRUE)

  # Relative power should sum to 1
  total_rel <- sum(bp_rel[1, -1])  # Exclude channel column
  expect_equal(total_rel, 1, tolerance = 0.01)
})

test_that("hilbertTransform computes analytic signal", {
  x <- make_test_signal()

  y <- hilbertTransform(x, output_assay = "analytic")

  # Check that analytic assay was created
  expect_true("analytic" %in% SummarizedExperiment::assayNames(y))

  analytic <- SummarizedExperiment::assay(y, "analytic")
  expect_true(is.complex(analytic))
})

test_that("instantaneousAmplitude extracts envelope", {
  x <- make_test_signal()

  x <- hilbertTransform(x, output_assay = "analytic")
  y <- instantaneousAmplitude(x, assay_name = "analytic", output_assay = "amplitude")

  expect_true("amplitude" %in% SummarizedExperiment::assayNames(y))

  amp <- SummarizedExperiment::assay(y, "amplitude")
  expect_true(is.numeric(amp))
  expect_true(all(amp >= 0))  # Amplitude should be non-negative
})

test_that("instantaneousPhase extracts phase", {
  x <- make_test_signal()

  x <- hilbertTransform(x, output_assay = "analytic")
  y <- instantaneousPhase(x, assay_name = "analytic", output_assay = "phase")

  expect_true("phase" %in% SummarizedExperiment::assayNames(y))

  phase <- SummarizedExperiment::assay(y, "phase")
  expect_true(is.numeric(phase))
  # Phase should be in [-pi, pi]
  expect_true(all(phase >= -pi & phase <= pi))
})

test_that("spectrogram handles missing sampling rate",
  {
  x <- make_test_signal()
  # Use NA_real_ instead of NA (logical) to match numeric slot type
  samplingRate(x) <- NA_real_

  expect_error(spectrogram(x), "sampling rate")
})

test_that("waveletTransform handles 3D data", {
  # Create 3D data (samples x channels x epochs)
  x <- make_test_signal(sr = 250, duration = 1)
  data <- SummarizedExperiment::assay(x)
  data3d <- array(data, dim = c(nrow(data), 1, 2))
  SummarizedExperiment::assay(x, "raw") <- data3d

  wt <- waveletTransform(x, frequencies = c(5, 10, 20))

  expect_type(wt, "list")
  expect_equal(nrow(wt$power), 3)  # 3 frequencies
})

test_that("oneWayAnova reproduces base R oneway.test(var.equal=TRUE) bit-for-bit", {
  set.seed(42)
  groups <- list(a = rnorm(30), b = rnorm(40, 0.6), c = rnorm(25, 1.2), d = rnorm(35, 0.3))
  op <- oneWayAnova(groups)

  y  <- unlist(groups, use.names = FALSE)
  gf <- factor(rep(names(groups), vapply(groups, length, integer(1))), levels = names(groups))
  ref <- stats::oneway.test(y ~ gf, var.equal = TRUE)

  # bit-for-bit against the base R equal-variance one-way ANOVA
  expect_equal(op$statistic, unname(ref$statistic), tolerance = 0)
  expect_equal(op$p_value,   ref$p.value,           tolerance = 0)
  # degrees of freedom and bookkeeping
  expect_identical(op$df1, 3L)                       # k - 1
  expect_identical(op$df2, length(y) - 4L)           # N - k
  expect_identical(op$k, 4L)
  expect_identical(op$n, length(y))

  # single-vector-plus-dots calling form equals the list form
  op2 <- oneWayAnova(groups$a, groups$b, groups$c, groups$d)
  expect_equal(op2$statistic, op$statistic, tolerance = 0)

  # input guards
  expect_error(oneWayAnova(list(1:5)), "at least 2 groups")
  expect_error(oneWayAnova(list(1:5, 3)), "at least 2 observations")
})

test_that("friedmanTest reproduces base R friedman.test bit-for-bit", {
  set.seed(11)
  # 12 blocks x 4 treatments with a shift + ties
  M <- round(matrix(rnorm(48), nrow = 12) + rep(c(0, 0.5, 1, 1.6), each = 12), 1)
  op <- friedmanTest(M)
  ref <- stats::friedman.test(M)

  expect_equal(op$statistic, unname(ref$statistic), tolerance = 0)   # tie-corrected chi-square, bit-for-bit
  expect_equal(op$p_value,   ref$p.value,           tolerance = 0)
  expect_identical(op$df, 3L)                                        # k - 1
  expect_identical(op$k, 4L)
  expect_identical(op$n, 12L)

  # list-of-columns form equals the matrix form
  op2 <- friedmanTest(lapply(seq_len(ncol(M)), function(j) M[, j]))
  expect_equal(op2$statistic, op$statistic, tolerance = 0)

  # rows with a missing value are dropped (complete blocks only)
  M2 <- M; M2[1, 2] <- NA
  expect_identical(friedmanTest(M2)$n, 11L)

  # input guards
  expect_error(friedmanTest(matrix(1:4, nrow = 4)), "at least 2 treatments")
})

test_that("kruskalTest reproduces base R kruskal.test bit-for-bit", {
  set.seed(23)
  groups <- list(a = round(rnorm(30), 1), b = round(rnorm(40, 0.6), 1),
                 c = round(rnorm(25, 1.2), 1), d = round(rnorm(35, 0.3), 1))  # ties
  op <- kruskalTest(groups)
  y  <- unlist(groups, use.names = FALSE)
  gf <- factor(rep(names(groups), vapply(groups, length, integer(1))), levels = names(groups))
  ref <- stats::kruskal.test(y ~ gf)

  expect_equal(op$statistic, unname(ref$statistic), tolerance = 0)   # tie-corrected H, bit-for-bit
  expect_equal(op$p_value,   ref$p.value,           tolerance = 0)
  expect_identical(op$df, 3L)
  expect_identical(op$k, 4L)
  expect_identical(op$n, length(y))

  op2 <- kruskalTest(groups$a, groups$b, groups$c, groups$d)
  expect_equal(op2$statistic, op$statistic, tolerance = 0)
  expect_error(kruskalTest(list(1:5)), "at least 2 groups")
})

test_that("linearRegression reproduces base R lm and scipy-style OLS", {
  set.seed(31); x <- rnorm(80); y <- 1.5 - 0.7 * x + rnorm(80, 0, 0.9)
  op <- linearRegression(x, y)
  co <- summary(stats::lm(y ~ x))$coefficients
  expect_equal(op$slope,     unname(co[2, 1]), tolerance = 1e-9)   # slope vs lm (QR) to machine precision
  expect_equal(op$intercept, unname(co[1, 1]), tolerance = 1e-9)
  expect_equal(op$statistic, unname(co[2, 3]), tolerance = 1e-7)   # slope t
  expect_equal(op$p_value,   unname(co[2, 4]), tolerance = 1e-7)
  expect_equal(op$r_squared, summary(stats::lm(y ~ x))$r.squared, tolerance = 1e-9)
  expect_identical(op$df, 78L); expect_identical(op$n, 80L)
  # two-column-x calling form equals the (x, y) form
  op2 <- linearRegression(cbind(x, y))
  expect_equal(op2$slope, op$slope, tolerance = 0)
  expect_error(linearRegression(1:2, 1:2), "at least 3")
})

test_that("kendallTau reproduces base R cor.test(method='kendall') bit-for-bit", {
  set.seed(41)
  x <- round(rnorm(60), 1); y <- round(0.5 * x + rnorm(60), 1)   # ties in both
  op <- kendallTau(x, y)
  ref <- suppressWarnings(stats::cor.test(x, y, method = "kendall"))
  expect_equal(op$statistic, unname(ref$estimate),  tolerance = 0)   # tau-b, bit-for-bit
  expect_equal(op$p_value,   ref$p.value,           tolerance = 0)   # normal-approx p, bit-for-bit
  expect_equal(op$z,         unname(ref$statistic), tolerance = 1e-9)
  expect_identical(op$n, 60L)
  # two-column-x calling form equals the (x, y) form
  op2 <- kendallTau(cbind(x, y))
  expect_equal(op2$statistic, op$statistic, tolerance = 0)
  expect_error(kendallTau(1:2, 1:2), "at least 3")
})

test_that("partialCorrelation matches the residual method and closed-form formula", {
  set.seed(51); z <- rnorm(80); x <- 0.6 * z + rnorm(80); y <- 0.5 * z + rnorm(80)
  op <- partialCorrelation(x, y, z, method = "pearson")
  # residual method (base R lm): correlate residuals of x~z and y~z
  ref <- stats::cor(stats::resid(stats::lm(x ~ z)), stats::resid(stats::lm(y ~ z)))
  expect_equal(op$statistic, ref, tolerance = 1e-9)
  # closed-form single-covariate formula
  rxy <- cor(x, y); rxz <- cor(x, z); ryz <- cor(y, z)
  cf <- (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
  expect_equal(op$statistic, cf, tolerance = 1e-9)
  expect_identical(op$df, 77L); expect_identical(op$k, 1L); expect_identical(op$n, 80L)
  # matrix (x, y, z) calling form equals the (x, y, z) form
  op2 <- partialCorrelation(cbind(x, y, z), method = "pearson")
  expect_equal(op2$statistic, op$statistic, tolerance = 0)
  # spearman branch runs and returns a finite coefficient
  expect_true(is.finite(partialCorrelation(x, y, z, method = "spearman")$statistic))
})

test_that("tukeyHSD reproduces base R TukeyHSD", {
  set.seed(61)
  g <- list(a = rnorm(30), b = rnorm(35, 0.8), c = rnorm(40, 1.6), d = rnorm(25, 0.4))
  op <- tukeyHSD(g); rownames(op) <- op$comparison
  y <- unlist(g, use.names = FALSE)
  grp <- factor(rep(names(g), vapply(g, length, integer(1))), levels = names(g))
  th <- stats::TukeyHSD(stats::aov(y ~ grp))$grp
  expect_equal(max(abs(op[rownames(th), "diff"] - th[, "diff"])),  0, tolerance = 1e-9)
  expect_equal(max(abs(op[rownames(th), "p_adj"] - th[, "p adj"])), 0, tolerance = 1e-9)
  expect_equal(max(abs(op[rownames(th), "lwr"] - th[, "lwr"])),    0, tolerance = 1e-9)
  expect_identical(nrow(op), 6L)             # 4 choose 2
  expect_identical(attr(op, "df"), 126L)     # N - k
  expect_error(tukeyHSD(list(rnorm(5))), "at least 2 groups")
})

test_that("dunnTest reproduces the Dunn z and Bonferroni p (base R recompute)", {
  set.seed(71)
  g <- list(a = rnorm(30), b = rnorm(35, 0.7), c = rnorm(40, 1.5), d = rnorm(25, 0.3))
  op <- dunnTest(g, p_adjust = "bonferroni"); rownames(op) <- op$comparison
  # independent recompute of one z (b-a): (Rbar_b - Rbar_a)/(sigma*sqrt(1/nb+1/na))
  y <- unlist(g, use.names = FALSE); grp <- rep(names(g), vapply(g, length, integer(1)))
  N <- length(y); r <- rank(y); tt <- table(y)
  sigma <- sqrt(N * (N + 1) / 12 - sum(tt^3 - tt) / (12 * (N - 1)))
  Rb <- tapply(r, grp, mean); nn <- tapply(r, grp, length)
  z_ba <- (Rb["b"] - Rb["a"]) / (sigma * sqrt(1 / nn["b"] + 1 / nn["a"]))
  expect_equal(op$z[op$comparison == "b-a"], unname(z_ba), tolerance = 1e-12)
  # Bonferroni = min(1, raw * n_pairs)
  expect_equal(op$p_adj, pmin(1, op$p_value * nrow(op)), tolerance = 0)
  expect_identical(nrow(op), 6L)             # 4 choose 2
  expect_identical(attr(op, "n"), N)
  expect_error(dunnTest(list(rnorm(5))), "at least 2 groups")
})

test_that("nemenyiTest reproduces the Nemenyi (Friedman) statistic and p", {
  set.seed(81)
  M <- matrix(rnorm(60), nrow = 15) + rep(c(0, 0.4, 0.9, 1.5), each = 15)
  colnames(M) <- c("a", "b", "c", "d")
  op <- nemenyiTest(M); rownames(op) <- op$comparison
  N <- nrow(M); k <- ncol(M); Rbar <- colMeans(t(apply(M, 1, rank))); se <- sqrt(k * (k + 1) / (6 * N))
  q_dc <- abs(Rbar["d"] - Rbar["c"]) / se
  expect_equal(op$statistic[op$comparison == "d vs c"], unname(q_dc), tolerance = 1e-12)
  # p_adj is ptukey of the op's own statistic (self-consistent mapping)
  expect_equal(op$p_adj[op$comparison == "d vs c"],
               stats::ptukey(op$statistic[op$comparison == "d vs c"] * sqrt(2), k, Inf, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_identical(nrow(op), 6L)             # 4 choose 2
  expect_identical(attr(op, "n"), N)
  expect_error(nemenyiTest(matrix(1:4, ncol = 1)), "at least 2 treatments")
})

test_that("multipleRegression reproduces base R lm", {
  set.seed(91); n <- 80; a <- rnorm(n); b <- rnorm(n); y <- 1 + 2 * a - 1.3 * b + rnorm(n)
  M <- cbind(y = y, a = a, b = b)
  op <- multipleRegression(M); rownames(op) <- op$term
  s <- summary(stats::lm(y ~ a + b))
  expect_equal(op$estimate, unname(s$coefficients[, 1]), tolerance = 1e-9)      # coef bit-for-bit
  expect_equal(op$std_error, unname(s$coefficients[, 2]), tolerance = 1e-9)
  expect_equal(op$p_value, unname(s$coefficients[, 4]), tolerance = 1e-9)
  expect_equal(attr(op, "r_squared"), s$r.squared, tolerance = 1e-9)
  expect_equal(attr(op, "f_statistic"), unname(s$fstatistic[1]), tolerance = 1e-7)
  expect_identical(attr(op, "df2"), 77L); expect_identical(attr(op, "n"), 80L)
  expect_error(multipleRegression(cbind(y, a)), "predictors")
})

test_that("cronbachAlpha reproduces the raw-alpha formula", {
  set.seed(101); k <- 6; n <- 120
  latent <- rnorm(n); X <- sapply(1:k, function(j) round(latent + rnorm(n, 0, 0.8)))
  op <- cronbachAlpha(X)
  vi <- apply(X, 2, var); vt <- var(rowSums(X))
  alpha_ref <- (k / (k - 1)) * (1 - sum(vi) / vt)     # independent recompute
  expect_equal(op$alpha, alpha_ref, tolerance = 1e-12)
  R <- cor(X); expect_equal(op$average_r, mean(R[upper.tri(R)]), tolerance = 1e-12)
  expect_identical(op$k, 6L); expect_identical(op$n, 120L)
  expect_error(cronbachAlpha(matrix(rnorm(10), ncol = 1)), "at least 2 items")
})

test_that("intraclassCorrelation reproduces the McGraw-Wong ICC formulas", {
  set.seed(111); n <- 30; k <- 3; subj <- rnorm(n, 0, 2)
  M <- sapply(1:k, function(j) subj + rnorm(n, 0, 0.6) + (j - 1) * 0.15)
  op <- intraclassCorrelation(M); rownames(op) <- op$form
  gm <- mean(M); SST <- sum((M - gm)^2); SSR <- k * sum((rowMeans(M) - gm)^2); SSC <- n * sum((colMeans(M) - gm)^2)
  MSR <- SSR / (n - 1); MSC <- SSC / (k - 1); MSE <- (SST - SSR - SSC) / ((n - 1) * (k - 1)); MSW <- (SST - SSR) / (n * (k - 1))
  expect_equal(op$ICC[op$form == "ICC2"], (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n), tolerance = 1e-12)
  expect_equal(op$ICC[op$form == "ICC3"], (MSR - MSE) / (MSR + (k - 1) * MSE), tolerance = 1e-12)
  expect_equal(op$ICC[op$form == "ICC2k"], (MSR - MSE) / (MSR + (MSC - MSE) / n), tolerance = 1e-12)
  expect_identical(nrow(op), 6L); expect_identical(attr(op, "n"), 30L); expect_identical(attr(op, "k"), 3L)
  expect_error(intraclassCorrelation(matrix(rnorm(10), ncol = 1)), "at least 2 raters")
})

test_that("leveneTest (Brown-Forsythe) reproduces the ANOVA-on-deviations formula", {
  set.seed(121)
  g <- list(a = rnorm(30, 0, 1), b = rnorm(35, 0, 2), c = rnorm(40, 0, 0.6), d = rnorm(25, 0, 1.3))
  op <- PhysioAnalysis::leveneTest(g, center = "median")
  # independent recompute: one-way ANOVA F on |x - median|
  z <- lapply(g, function(v) abs(v - median(v)))
  y <- unlist(z); grp <- factor(rep(seq_along(z), vapply(z, length, integer(1))))
  ref <- summary(stats::aov(y ~ grp))[[1]][["F value"]][1]
  expect_equal(op$statistic, ref, tolerance = 1e-9)
  expect_identical(op$df1, 3L); expect_identical(op$df2, length(y) - 4L)
  expect_identical(op$center, "median")
  # mean-centred classic Levene also runs
  expect_true(is.finite(PhysioAnalysis::leveneTest(g, center = "mean")$statistic))
  expect_error(PhysioAnalysis::leveneTest(list(1:5)), "at least 2 groups")
})
