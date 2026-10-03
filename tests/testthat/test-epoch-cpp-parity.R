library(testthat)
library(PhysioAnalysis)

# Numerical parity of the 0.3.5 compiled fixed-length epoch extractor against
# the pure-R gather reference (.epochFixedR), plus edge cases.

make_pe <- function(n_time = 5000, n_ch = 8, sr = 256, n_events = 40) {
  set.seed(31)
  X <- matrix(rnorm(n_time * n_ch), ncol = n_ch)
  colnames(X) <- paste0("ch", seq_len(n_ch))
  pe <- PhysioExperiment(assays = list(raw = X), samplingRate = sr)
  on <- seq(1, (n_time / sr) - 1, length.out = n_events)
  addEvents(pe, onset = on, type = "stim")
}

test_that("epochData fixed-length output matches the pure-R gather", {
  sr <- 256
  pe <- make_pe(sr = sr)
  X <- SummarizedExperiment::assay(pe, "raw")
  ep <- epochData(pe, tmin = -0.2, tmax = 0.8, event_type = "stim")
  got <- SummarizedExperiment::assay(ep, SummarizedExperiment::assayNames(ep)[1])

  onsets <- getEvents(pe)@events$onset
  es <- timeToSamples(pe, onsets)
  pre <- as.integer(round(0.2 * sr)); post <- as.integer(round(0.8 * sr))
  vi <- which(es - pre >= 1 & es + post <= nrow(X))
  d3 <- array(X, dim = c(nrow(X), ncol(X), 1))
  ref <- PhysioAnalysis:::.epochFixedR(d3, es[vi], pre, post)

  expect_equal(dim(got), dim(ref))
  expect_equal(dim(got), c(pre + post + 1L, ncol(X), length(vi), 1L))
  expect_equal(max(abs(got - ref)), 0)       # exact values (dimnames differ)
})

test_that("epochData matches gather for 3D (multi-sample) input", {
  set.seed(32)
  sr <- 200
  A <- array(rnorm(3000 * 4 * 3), dim = c(3000, 4, 3))
  pe <- PhysioExperiment(assays = list(raw = A), samplingRate = sr)
  pe <- addEvents(pe, onset = seq(1, 13, length.out = 20), type = "stim")
  ep <- epochData(pe, tmin = -0.1, tmax = 0.4, event_type = "stim")
  got <- SummarizedExperiment::assay(ep, SummarizedExperiment::assayNames(ep)[1])

  es <- timeToSamples(pe, getEvents(pe)@events$onset)
  pre <- as.integer(round(0.1 * sr)); post <- as.integer(round(0.4 * sr))
  vi <- which(es - pre >= 1 & es + post <= dim(A)[1])
  ref <- PhysioAnalysis:::.epochFixedR(A, es[vi], pre, post)
  expect_equal(max(abs(got - ref)), 0)
})

test_that("cpp_epoch_fixed copies the correct windows (direct check)", {
  set.seed(33)
  n_time <- 1000; n_ch <- 3
  X <- matrix(seq_len(n_time * n_ch) + 0.0, nrow = n_time)  # distinct values
  centers <- c(100L, 500L, 900L)
  pre <- 10L; post <- 20L
  out <- PhysioAnalysis:::cpp_epoch_fixed(X, n_time, n_ch, 1L, centers, pre, post)
  dim(out) <- c(pre + post + 1L, n_ch, length(centers), 1L)
  # epoch 2 (center 500), channel 3, spans source rows 490..520
  expect_equal(out[, 3, 2, 1], X[490:520, 3])
  expect_equal(out[, 1, 1, 1], X[90:120, 1])
})

test_that("cpp_epoch_fixed guards out-of-bounds and NA centers", {
  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  expect_error(PhysioAnalysis:::cpp_epoch_fixed(X, 1000L, 3L, 1L, 1e8L, 10L, 20L),
               "outside the signal bounds")
  expect_error(PhysioAnalysis:::cpp_epoch_fixed(X, 1000L, 3L, 1L, NA_integer_, 10L, 20L),
               "outside the signal bounds")
  expect_error(PhysioAnalysis:::cpp_epoch_fixed(X, 1000L, 3L, 1L, 5L, 10L, 20L),
               "outside the signal bounds")       # center-pre < 1
  # valid centers succeed
  out <- PhysioAnalysis:::cpp_epoch_fixed(X, 1000L, 3L, 1L, c(100L, 500L), 10L, 20L)
  expect_length(out, 31L * 3L * 2L * 1L)
})

test_that("epochData still applies baseline correction on the fast path", {
  sr <- 256
  pe <- make_pe(sr = sr)
  ep <- epochData(pe, tmin = -0.2, tmax = 0.8, event_type = "stim",
                  baseline = c(-0.2, 0))
  got <- SummarizedExperiment::assay(ep, SummarizedExperiment::assayNames(ep)[1])
  # baseline window mean per epoch/channel should be ~0 after correction
  bl_end <- as.integer(round((0 - (-0.2)) * sr)) + 1L
  bl_means <- apply(got[1:bl_end, , , , drop = FALSE], c(2, 3, 4), mean)
  expect_lt(max(abs(bl_means)), 1e-9)
})
