library(testthat)
library(PhysioAnalysis)

.gauss <- function(f, c, h, w) h * exp(-(f - c)^2 / (2 * w^2))

# Synthetic PSD = 10^(aperiodic + Gaussian peaks).
syn_psd <- function(off = 1, expo = 1.5,
                    peaks = list(c(10, 0.7, 1.2), c(22, 0.4, 2.0)),
                    fr = seq(1, 45, by = 0.5), mode = "fixed", knee = NULL) {
  lp <- if (mode == "fixed") off - expo * log10(fr) else
    off - log10(knee + fr^expo)
  for (p in peaks) lp <- lp + .gauss(fr, p[1], p[2], p[3])
  list(freqs = fr, power = 10^lp)
}

fit1 <- function(s, mode = "fixed", ...) {
  PhysioAnalysis:::.specparam_fit_one(
    s$freqs, s$power, freq_range = range(s$freqs),
    peak_width_limits = c(1, 12), max_n_peaks = 6L,
    aperiodic_mode = mode, min_peak_height = 0.05, peak_threshold = 2, ...)
}

test_that("recovers 1/f exponent within 0.1 and peak CF within 0.5 Hz", {
  s <- syn_psd(off = 1, expo = 1.5)
  f <- fit1(s, mode = "fixed")
  expect_lt(abs(f$aperiodic[["exponent"]] - 1.5), 0.1)
  expect_gt(f$r_squared, 0.95)
  cfs <- sort(f$peaks$CF)
  expect_equal(length(cfs), 2L)
  expect_lt(abs(cfs[1] - 10), 0.5)
  expect_lt(abs(cfs[2] - 22), 0.5)
})

test_that("R-squared > 0.95 and graceful handling when no peaks present", {
  set.seed(1)
  fr <- seq(1, 45, by = 0.5)
  lp <- 1 - 1.2 * log10(fr) + stats::rnorm(length(fr), sd = 0.008)  # smooth 1/f
  f <- fit1(list(freqs = fr, power = 10^lp), mode = "fixed")
  expect_equal(nrow(f$peaks), 0L)
  expect_gt(f$r_squared, 0.95)
  expect_lt(abs(f$aperiodic[["exponent"]] - 1.2), 0.1)
})

test_that("knee mode recovers the knee frequency on a synthetic knee spectrum", {
  s <- syn_psd(off = 1, expo = 2.0, knee = 100, mode = "knee",
               peaks = list(c(10, 0.5, 1.5)))
  f <- fit1(s, mode = "knee")
  expect_gt(f$r_squared, 0.95)
  expect_lt(abs(f$aperiodic[["exponent"]] - 2.0), 0.2)
  # knee recovered within 25% of the true value (100)
  expect_lt(abs(f$aperiodic[["knee"]] - 100) / 100, 0.25)
})

# ---- full wrapper on a 1/f signal with an alpha oscillation ----
make_1f_pe <- function(expo = 1.2, sr = 250, n = 8000, alpha_amp = 8,
                       n_channels = 3, seed = 1) {
  set.seed(seed)
  t <- (0:(n - 1)) / sr
  freq <- (0:(n - 1)) * sr / n; freq[1] <- freq[2]
  one_over_f <- function() {
    amp <- freq^(-expo / 2)
    spec <- amp * exp(1i * stats::runif(n, 0, 2 * pi))
    sig <- Re(stats::fft(spec, inverse = TRUE)) / n
    sig / stats::sd(sig)
  }
  mat <- sapply(seq_len(n_channels), function(j)
    20 * one_over_f() + alpha_amp * sin(2 * pi * 10 * t + stats::runif(1, 0, 2 * pi)))
  PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(label = paste0("C", seq_len(n_channels)),
                                   type = rep("EEG", n_channels)),
    samplingRate = sr)
}

test_that("specparam() fits a PhysioExperiment and detects the alpha peak", {
  pe <- make_1f_pe(expo = 1.2)
  sp <- specparam(pe, freq_range = c(2, 40))
  expect_s3_class(sp, "specparam_result")
  expect_equal(nrow(sp$aperiodic), 3L)
  expect_true(all(sp$aperiodic$exponent > 0.6))          # recovered ~1.2
  expect_true(all(sp$fit$r_squared > 0.9))
  # an alpha peak near 10 Hz on at least one channel
  alpha_peaks <- sp$peaks[sp$peaks$CF >= 8 & sp$peaks$CF <= 12, ]
  expect_gt(nrow(alpha_peaks), 0L)
})

test_that("plotSpecparam returns a ggplot and specparamBiomarker a PhysioBiomarker", {
  pe <- make_1f_pe()
  sp <- specparam(pe, freq_range = c(2, 40))
  p <- plotSpecparam(sp, channel = 1)
  expect_s3_class(p, "ggplot")

  bm <- specparamBiomarker(sp, channel = 1, param = "exponent")
  expect_true(is.PhysioBiomarker(bm))
  expect_equal(bm@name, "aperiodic_exponent")
  expect_equal(bm@provenance_info$method, "specparam:fixed")
  expect_equal(biomarkerValue(bm), sp$aperiodic$exponent[1])
})

test_that("print method and knee biomarker guard work", {
  pe <- make_1f_pe()
  sp <- specparam(pe, freq_range = c(2, 40))
  expect_output(print(sp), "specparam_result")
  expect_error(specparamBiomarker(sp, param = "knee"), "fixed")
})

test_that("specparam validates inputs", {
  pe <- make_1f_pe()
  expect_error(specparam(pe, freq_range = c(40, 2)), "freq_range")
  expect_error(specparam(list()), "PhysioExperiment")
})
