library(testthat)
library(PhysioAnalysis)

# Numerical parity of the 0.3.4 compiled/batched spectral paths against the
# pure-R reference implementations they replaced.

test_that("fftSignals (mvfft batch) matches the per-column apply reference", {
  set.seed(21)
  X <- matrix(rnorm(2048 * 8), ncol = 8)
  colnames(X) <- paste0("ch", 1:8)
  pe <- PhysioExperiment(assays = list(raw = X), samplingRate = 256)

  got <- SummarizedExperiment::assay(fftSignals(pe), "fft")
  want <- apply(X, 2, function(v) Mod(stats::fft(v)))
  expect_equal(dim(got), dim(X))
  expect_lt(max(abs(got - want)), 1e-10)
  expect_identical(dimnames(got), dimnames(X))
})

test_that("cpp_stft_power matches the pure-R per-window loop (power-of-two)", {
  set.seed(22)
  sr <- 256
  sig <- rnorm(4096)
  for (ws in c(64L, 128L, 256L, 512L)) {
    for (ov in c(0, 0.5, 0.75)) {
      step <- as.integer(ws * (1 - ov))
      nf <- ws %/% 2L + 1L
      nw <- floor((length(sig) - ws) / step) + 1L
      wf <- PhysioAnalysis:::.getWindow("hanning", ws)
      wp <- sum(wf^2)
      got <- PhysioAnalysis:::cpp_stft_power(sig, wf, step, sr, wp)
      want <- PhysioAnalysis:::.spectrogramLoop(sig, wf, step, sr, wp, ws, nf, nw)
      expect_equal(dim(got), c(nf, nw))
      expect_lt(max(abs(got - want)), 1e-9)
    }
  }
})

test_that("spectrogram() dispatches to matching results for both window types", {
  set.seed(23)
  X <- matrix(rnorm(3000 * 3), ncol = 3)
  pe <- PhysioExperiment(assays = list(raw = X), samplingRate = 250)

  # power-of-two window -> compiled kernel
  sp2 <- spectrogram(pe, window_size = 256L, overlap = 0.5, channel = 2L)
  wf2 <- PhysioAnalysis:::.getWindow("hanning", 256L)
  nw2 <- ncol(sp2$power)
  ref2 <- PhysioAnalysis:::.spectrogramLoop(X[, 2], wf2, 128L, 250,
                                            sum(wf2^2), 256L, 129L, nw2)
  expect_lt(max(abs(sp2$power - ref2)), 1e-9)

  # non-power-of-two window -> mvfft fallback
  sp3 <- spectrogram(pe, window_size = 200L, overlap = 0.5, channel = 2L)
  wf3 <- PhysioAnalysis:::.getWindow("hanning", 200L)
  nw3 <- ncol(sp3$power)
  ref3 <- PhysioAnalysis:::.spectrogramLoop(X[, 2], wf3, 100L, 250,
                                            sum(wf3^2), 200L, 101L, nw3)
  expect_lt(max(abs(sp3$power - ref3)), 1e-9)
})

test_that("spectrogram rejects too-high overlap instead of crashing (C1)", {
  set.seed(25)
  X <- matrix(rnorm(2000), ncol = 1)
  pe <- PhysioExperiment(assays = list(raw = X), samplingRate = 256)
  # step_size = as.integer(256 * (1 - overlap)) < 1 -> clean error, not SIGFPE
  expect_error(spectrogram(pe, window_size = 256L, overlap = 1),
               "too high")
  expect_error(spectrogram(pe, window_size = 256L, overlap = 0.999),
               "too high")
  # kernel guards step and short signals directly
  expect_error(PhysioAnalysis:::cpp_stft_power(rnorm(512), rep(1, 128), 0L,
                                               256, 1), "step")
  empty <- PhysioAnalysis:::cpp_stft_power(rnorm(100), rep(1, 128), 64L, 256, 1)
  expect_equal(dim(empty), c(65L, 0L))
})

test_that("fftSignals handles 3D single-channel data (C2)", {
  set.seed(26)
  a1 <- array(rnorm(256 * 1 * 3), dim = c(256, 1, 3))
  pe1 <- PhysioExperiment(assays = list(raw = a1), samplingRate = 256)
  got <- SummarizedExperiment::assay(fftSignals(pe1), "fft")
  want <- apply(a1, c(2, 3), function(v) Mod(stats::fft(v)))
  expect_equal(dim(got), dim(a1))
  expect_lt(max(abs(got - want)), 1e-10)
  # multi-channel 3D still correct
  a2 <- array(rnorm(256 * 4 * 2), dim = c(256, 4, 2))
  pe2 <- PhysioExperiment(assays = list(raw = a2), samplingRate = 256)
  g2 <- SummarizedExperiment::assay(fftSignals(pe2), "fft")
  w2 <- apply(a2, c(2, 3), function(v) Mod(stats::fft(v)))
  expect_lt(max(abs(g2 - w2)), 1e-10)
})

test_that(".spectrogramLoop reference is correct for tiny windows (C4)", {
  # ws=2 (nf=2): neither DC nor Nyquist should be doubled. The reference must
  # not double them (old code did) and must not error.
  set.seed(27)
  sig <- rnorm(64)
  wf <- PhysioAnalysis:::.getWindow("hamming", 2L)   # c(0.08, 0.08), nonzero
  nw <- floor((64 - 2) / 1) + 1
  ref <- PhysioAnalysis:::.spectrogramLoop(sig, wf, 1L, 100, sum(wf^2),
                                           2L, 2L, nw)
  got <- PhysioAnalysis:::cpp_stft_power(sig, wf, 1L, 100, sum(wf^2))
  expect_equal(dim(ref), c(2L, nw))
  expect_lt(max(abs(ref - got)), 1e-12)              # reference == kernel now
})

test_that("spectrogram Parseval / frequency structure is preserved", {
  set.seed(24)
  sr <- 256
  t <- seq_len(2048) / sr
  sig <- sin(2 * pi * 30 * t)                 # pure 30 Hz tone
  X <- matrix(sig, ncol = 1)
  pe <- PhysioExperiment(assays = list(raw = X), samplingRate = sr)
  sp <- spectrogram(pe, window_size = 256L, overlap = 0.5, channel = 1L)
  # peak frequency bin should sit at 30 Hz
  peak_hz <- sp$frequencies[which.max(rowMeans(sp$power))]
  expect_lt(abs(peak_hz - 30), sr / 256)
})
