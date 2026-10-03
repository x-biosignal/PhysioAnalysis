#' Time-Frequency Analysis for PhysioExperiment
#'
#' Functions for time-frequency analysis including wavelet transforms,
#' spectrograms, and band power extraction.

#' Compute spectrogram (Short-Time Fourier Transform)
#'
#' Computes the spectrogram using STFT with proper power spectral density
#' normalization. The returned power values are one-sided PSD estimates in
#' V^2/Hz (non-DC/Nyquist bins are doubled).
#'
#' @param x A PhysioExperiment object.
#' @param window_size Window size in samples.
#' @param overlap Overlap between windows (0-1).
#' @param window_type Window function: "hanning", "hamming", "blackman", or "rectangular".
#' @param channel Channel index to analyze.
#' @param sample Sample index (for 3D data).
#' @return A list with the following components:
#'   \item{power}{Power spectrogram matrix (frequency x time)}
#'   \item{frequencies}{Numeric vector of frequencies in Hz}
#'   \item{times}{Numeric vector of time points in seconds}
#'   \item{sampling_rate}{The sampling rate used}
#'   \item{window_size}{Window size in samples}
#'   \item{overlap}{Overlap fraction}
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [waveletTransform()] for wavelet-based time-frequency analysis,
#'   [plotSpectrogram()] for visualization, [bandPower()] for band power
#'   extraction, [fftSignals()] for simple FFT.
#' @export
#' @examples
#' # Create example data
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
#'   samplingRate = 250
#' )
#'
#' # Compute spectrogram for channel 1
#' spec <- spectrogram(pe, channel = 1)
#'
#' # Plot spectrogram
#' plotSpectrogram(spec, freq_range = c(1, 40))
spectrogram <- function(x, window_size = 256L, overlap = 0.5,
                        window_type = c("hanning", "hamming", "blackman", "rectangular"),
                        channel = 1L, sample = 1L) {
  stopifnot(inherits(x, "PhysioExperiment"))
  window_type <- match.arg(window_type)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  # Extract signal
  if (length(dims) == 2) {
    signal <- data[, channel]
  } else if (length(dims) == 3) {
    signal <- data[, channel, sample]
  } else {
    stop("Data must be 2D or 3D", call. = FALSE)
  }

  n <- length(signal)
  step_size <- as.integer(window_size * (1 - overlap))
  if (step_size < 1L) {
    stop(sprintf(paste0("overlap (%.3f) is too high for window_size %d: the ",
                        "hop would be < 1 sample. Reduce overlap or increase ",
                        "window_size."), overlap, window_size), call. = FALSE)
  }
  n_windows <- floor((n - window_size) / step_size) + 1

  # Generate window function
  window_func <- .getWindow(window_type, window_size)
  window_power <- sum(window_func^2)  # For Parseval's theorem compliance

  # Initialize output
  n_freqs <- floor(window_size / 2) + 1

  # Compute STFT. The compiled radix-2 kernel handles power-of-two windows
  # (the common case, e.g. the default 256); otherwise a batched mvfft path is
  # used. Both reproduce the pure-R per-window loop (retained as
  # .spectrogramLoop() for parity tests) exactly.
  power_matrix <- .spectrogramPower(signal, window_func, step_size, sr,
                                    window_power, window_size, n_freqs,
                                    n_windows)

  starts <- (seq_len(n_windows) - 1L) * step_size + 1L
  times <- (starts + (starts + window_size - 1L)) / 2 / sr

  # Frequency vector
  frequencies <- seq(0, sr / 2, length.out = n_freqs)

  list(
    power = power_matrix,
    frequencies = frequencies,
    times = times,
    sampling_rate = sr,
    window_size = window_size,
    overlap = overlap
  )
}

# STFT power dispatcher: compiled radix-2 kernel for power-of-two windows,
# batched mvfft otherwise. Both match .spectrogramLoop() to floating-point
# precision.
.spectrogramPower <- function(signal, window_func, step_size, sr,
                              window_power, window_size, n_freqs, n_windows) {
  if (n_windows <= 0) return(matrix(NA_real_, nrow = n_freqs, ncol = 0))
  if (bitwAnd(window_size, window_size - 1L) == 0L) {
    return(cpp_stft_power(as.numeric(signal), as.numeric(window_func),
                          as.integer(step_size), sr, window_power))
  }
  # Non-power-of-two fallback: batch all windows through a single mvfft.
  idx <- outer(seq_len(window_size) - 1L,
               (seq_len(n_windows) - 1L) * step_size, "+") + 1L
  segs <- matrix(signal[idx], nrow = window_size) * window_func
  psd <- (Mod(stats::mvfft(segs))[seq_len(n_freqs), , drop = FALSE])^2 /
    (sr * window_power)
  if (n_freqs > 2) psd[2:(n_freqs - 1), ] <- 2 * psd[2:(n_freqs - 1), ]
  psd
}

# Pure-R per-window reference implementation (pre-0.3.4), retained for
# numerical-parity tests against the compiled/batched paths.
.spectrogramLoop <- function(signal, window_func, step_size, sr,
                             window_power, window_size, n_freqs, n_windows) {
  power_matrix <- matrix(NA_real_, nrow = n_freqs, ncol = n_windows)
  for (i in seq_len(n_windows)) {
    start <- (i - 1) * step_size + 1
    end <- start + window_size - 1
    segment <- signal[start:end] * window_func
    psd <- Mod(stats::fft(segment)[1:n_freqs])^2 / (sr * window_power)
    if (n_freqs > 2) psd[2:(n_freqs - 1)] <- 2 * psd[2:(n_freqs - 1)]
    power_matrix[, i] <- psd
  }
  power_matrix
}

#' Wavelet transform
#'
#' Computes the continuous wavelet transform using Morlet wavelets.
#'
#' @param x A PhysioExperiment object.
#' @param frequencies Numeric vector of frequencies to analyze.
#' @param n_cycles Number of wavelet cycles (can be scalar or vector).
#' @param channel Channel index to analyze.
#' @param sample Sample index (for 3D data).
#' @param normalization Wavelet normalization method: \code{"L2"} (default,
#'   divides by square root of sum of squared absolute values) or \code{"L1"}
#'   (divides by sum of absolute values). L2 normalization preserves energy
#'   across frequencies and is preferred for power comparisons.
#' @return A list with the following components:
#'   \item{power}{Power matrix (frequency x time)}
#'   \item{phase}{Phase matrix in radians (frequency x time)}
#'   \item{frequencies}{Numeric vector of analyzed frequencies in Hz}
#'   \item{times}{Numeric vector of time points in seconds}
#'   \item{sampling_rate}{The sampling rate used}
#'   \item{n_cycles}{Number of wavelet cycles per frequency}
#' @references Torrence, C. & Compo, G.P. (1998). "A practical guide to wavelet
#'   analysis." Bulletin of the American Meteorological Society, 79(1), 61-78.
#' @seealso [spectrogram()] for STFT-based time-frequency analysis,
#'   [bandPower()] for band power extraction, [hilbertTransform()] for
#'   analytic signal computation.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
#'   samplingRate = 100
#' )
#'
#' # Compute wavelet transform (1-30 Hz)
#' wt <- waveletTransform(pe, frequencies = seq(1, 30), channel = 1)
#'
#' # Access power and phase
#' dim(wt$power)  # frequency x time
waveletTransform <- function(x, frequencies = seq(1, 40, by = 1),
                              n_cycles = 7, channel = 1L, sample = 1L,
                              normalization = c("L2", "L1")) {
  stopifnot(inherits(x, "PhysioExperiment"))
  normalization <- match.arg(normalization)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  # Extract signal
  if (length(dims) == 2) {
    signal <- data[, channel]
  } else if (length(dims) == 3) {
    signal <- data[, channel, sample]
  } else {
    stop("Data must be 2D or 3D", call. = FALSE)
  }

  n <- length(signal)
  n_freqs <- length(frequencies)

  # Handle n_cycles
  if (length(n_cycles) == 1) {
    n_cycles <- rep(n_cycles, n_freqs)
  }

  # Initialize output
  power_matrix <- matrix(NA_real_, nrow = n_freqs, ncol = n)
  phase_matrix <- matrix(NA_real_, nrow = n_freqs, ncol = n)

  # Compute wavelet transform for each frequency
  for (f_idx in seq_len(n_freqs)) {
    freq <- frequencies[f_idx]
    cycles <- n_cycles[f_idx]

    # Generate Morlet wavelet
    wavelet <- .morletWavelet(freq, sr, cycles, normalization = normalization)

    # Convolve signal with wavelet
    conv_result <- .convolveComplex(signal, wavelet)

    # Extract power and phase
    power_matrix[f_idx, ] <- Mod(conv_result)^2
    phase_matrix[f_idx, ] <- Arg(conv_result)
  }

  # Time vector
  times <- (seq_len(n) - 1) / sr

  list(
    power = power_matrix,
    phase = phase_matrix,
    frequencies = frequencies,
    times = times,
    sampling_rate = sr,
    n_cycles = n_cycles
  )
}

#' Generate Morlet wavelet
#' @noRd
.morletWavelet <- function(freq, sr, n_cycles,
                           normalization = c("L2", "L1")) {
  normalization <- match.arg(normalization)

  # Standard deviation of Gaussian
  sigma_t <- n_cycles / (2 * pi * freq)

  # Wavelet duration (6 sigma on each side)
  wavelet_duration <- 6 * sigma_t
  n_samples <- as.integer(2 * wavelet_duration * sr) + 1

  # Time vector centered at 0
  t <- seq(-wavelet_duration, wavelet_duration, length.out = n_samples)

  # Complex Morlet wavelet
  gaussian <- exp(-t^2 / (2 * sigma_t^2))
  sinusoid <- exp(2i * pi * freq * t)

  wavelet <- gaussian * sinusoid

  # Normalize
  if (normalization == "L1") {
    wavelet / sum(Mod(wavelet))
  } else {
    # L2 normalization for energy preservation
    wavelet / sqrt(sum(Mod(wavelet)^2))
  }
}

#' Complex convolution
#' @noRd
.convolveComplex <- function(signal, wavelet) {
  n_signal <- length(signal)
  n_wavelet <- length(wavelet)

  # Pad signal for convolution
  n_fft <- n_signal + n_wavelet - 1
  n_fft <- 2^ceiling(log2(n_fft))  # Pad to power of 2

  signal_padded <- c(signal, rep(0, n_fft - n_signal))
  wavelet_padded <- c(wavelet, rep(0, n_fft - n_wavelet))

  # FFT convolution
  fft_signal <- stats::fft(signal_padded)
  fft_wavelet <- stats::fft(wavelet_padded)
  conv_result <- stats::fft(fft_signal * fft_wavelet, inverse = TRUE) / n_fft

  # Extract valid part (centered)
  half_wavelet <- floor(n_wavelet / 2)
  conv_result[(half_wavelet + 1):(half_wavelet + n_signal)]
}

#' Get window function
#' @noRd
.getWindow <- function(type, n) {
  t <- seq(0, 1, length.out = n)

  switch(type,
    hanning = 0.5 * (1 - cos(2 * pi * t)),
    hamming = 0.54 - 0.46 * cos(2 * pi * t),
    blackman = 0.42 - 0.5 * cos(2 * pi * t) + 0.08 * cos(4 * pi * t),
    rectangular = rep(1, n)
  )
}

#' Compute band power
#'
#' Extracts power in specified frequency bands.
#'
#' @param x A PhysioExperiment object.
#' @param bands Named list of frequency bands. Each element should be c(low, high).
#'   Default includes standard EEG bands.
#' @param method Method: "welch" (PSD) or "wavelet".
#' @param relative If TRUE, returns relative power (proportion of total).
#' @return A `data.frame` with one row per channel and columns for the channel
#'   name and power in each frequency band. If `relative = TRUE`, values
#'   represent the proportion of total power.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [spectrogram()] for full time-frequency decomposition,
#'   [waveletTransform()] for wavelet-based power, [fftSignals()] for raw FFT,
#'   [plotPSD()] for power spectral density visualization.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(2560), nrow = 256, ncol = 10)),
#'   colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
#'   samplingRate = 256
#' )
#'
#' # Compute band power for standard EEG bands
#' bp <- bandPower(pe)
#' head(bp)
#'
#' # Compute relative band power
#' bp_rel <- bandPower(pe, relative = TRUE)
bandPower <- function(x, bands = NULL, method = c("welch", "wavelet"),
                      relative = FALSE) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  # Default EEG bands
  if (is.null(bands)) {
    bands <- list(
      delta = c(0.5, 4),
      theta = c(4, 8),
      alpha = c(8, 13),
      beta = c(13, 30),
      gamma = c(30, 100)
    )
  }

  sr <- samplingRate(x)
  assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  # Flatten to 2D
  if (length(dims) == 3) {
    data <- apply(data, c(1, 2), mean)
    dims <- dim(data)
  }

  n_channels <- dims[2]
  n_bands <- length(bands)
  band_names <- names(bands)

  # Initialize result
  result <- matrix(NA_real_, nrow = n_channels, ncol = n_bands)
  colnames(result) <- band_names

  for (ch in seq_len(n_channels)) {
    signal <- data[, ch]

    if (method == "welch") {
      # Welch's method for PSD estimation
      psd <- .welchPSD(signal, sr)

      for (b in seq_len(n_bands)) {
        freq_range <- bands[[b]]
        freq_idx <- which(psd$frequencies >= freq_range[1] &
                            psd$frequencies <= freq_range[2])
        result[ch, b] <- sum(psd$power[freq_idx])
      }

    } else if (method == "wavelet") {
      # Wavelet-based power
      freqs <- seq(0.5, min(100, sr / 2 - 1), by = 0.5)
      wt <- waveletTransform(x, frequencies = freqs, channel = ch)

      for (b in seq_len(n_bands)) {
        freq_range <- bands[[b]]
        freq_idx <- which(wt$frequencies >= freq_range[1] &
                            wt$frequencies <= freq_range[2])
        result[ch, b] <- mean(wt$power[freq_idx, ])
      }
    }
  }

  # Convert to relative power if requested
  if (relative) {
    row_sums <- rowSums(result)
    result <- result / row_sums
  }

  # Create result data.frame
  ch_names <- channelNames(x)
  if (length(ch_names) != n_channels) {
    ch_names <- paste0("Ch", seq_len(n_channels))
  }

  result_df <- as.data.frame(result)
  result_df$channel <- ch_names
  result_df <- result_df[, c("channel", band_names)]

  result_df
}

#' Welch's PSD estimation
#' @noRd
.welchPSD <- function(signal, sr, nperseg = 256, noverlap = NULL) {
  n <- length(signal)

  if (is.null(noverlap)) {
    noverlap <- floor(nperseg / 2)
  }

  step <- nperseg - noverlap
  n_segments <- floor((n - noverlap) / step)

  if (n_segments < 1) {
    nperseg <- n
    n_segments <- 1
    step <- nperseg
  }

  n_freqs <- floor(nperseg / 2) + 1
  psd <- numeric(n_freqs)

  window <- .getWindow("hanning", nperseg)
  window_sum <- sum(window^2)

  for (i in seq_len(n_segments)) {
    start <- (i - 1) * step + 1
    end <- start + nperseg - 1
    if (end > n) break

    segment <- signal[start:end] * window
    fft_result <- stats::fft(segment)
    psd <- psd + Mod(fft_result[1:n_freqs])^2
  }

  # Normalize
  psd <- psd / (n_segments * sr * window_sum)

  # Double power for one-sided spectrum (except DC and Nyquist)
  psd[2:(n_freqs - 1)] <- 2 * psd[2:(n_freqs - 1)]

  frequencies <- seq(0, sr / 2, length.out = n_freqs)

  list(power = psd, frequencies = frequencies)
}

#' Hilbert transform for instantaneous amplitude/phase
#'
#' Computes the analytic signal using the Hilbert transform.
#' The analytic signal can be used to extract instantaneous amplitude
#' and phase.
#'
#' @param x A PhysioExperiment object.
#' @param output_assay Name for the output assay.
#' @return A `PhysioExperiment` object with an additional assay (default
#'   `"analytic"`) containing the complex-valued analytic signal.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [instantaneousAmplitude()] to extract the signal envelope,
#'   [instantaneousPhase()] to extract the instantaneous phase,
#'   [plv()] for phase-based connectivity.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
#'   samplingRate = 100
#' )
#'
#' # Compute Hilbert transform
#' pe <- hilbertTransform(pe)
#'
#' # Extract amplitude and phase
#' pe <- instantaneousAmplitude(pe)
#' pe <- instantaneousPhase(pe)
hilbertTransform <- function(x, output_assay = "analytic") {
  stopifnot(inherits(x, "PhysioExperiment"))

  assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  hilbert_1d <- function(signal) {
    n <- length(signal)
    fft_signal <- stats::fft(signal)

    # Create Hilbert transform filter
    h <- numeric(n)
    h[1] <- 1
    if (n %% 2 == 0) {
      h[2:(n / 2)] <- 2
      h[n / 2 + 1] <- 1
    } else {
      h[2:((n + 1) / 2)] <- 2
    }

    # Apply filter and inverse FFT
    analytic <- stats::fft(fft_signal * h, inverse = TRUE) / n
    analytic
  }

  if (length(dims) == 2) {
    result <- apply(data, 2, hilbert_1d)
  } else if (length(dims) == 3) {
    result <- array(NA_complex_, dim = dims)
    for (s in seq_len(dims[3])) {
      # Use drop=FALSE to preserve 2D matrix when dims[3]=1
      slice <- data[, , s, drop = FALSE]
      dim(slice) <- dims[1:2]
      result[, , s] <- apply(slice, 2, hilbert_1d)
    }
  } else {
    stop("Data must be 2D or 3D", call. = FALSE)
  }

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- result
  SummarizedExperiment::assays(x) <- assays

  .recordProv(x, input_assay = assay_name, output_assay = output_assay,
              .package = "PhysioAnalysis")
}

#' Extract instantaneous amplitude (envelope)
#'
#' Extracts the instantaneous amplitude (envelope) from the analytic signal.
#'
#' @param x A PhysioExperiment object with analytic signal.
#' @param assay_name Name of the analytic signal assay.
#' @param output_assay Name for the output assay.
#' @return A `PhysioExperiment` object with an additional assay (default
#'   `"amplitude"`) containing the real-valued instantaneous amplitude
#'   (envelope) of the signal.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [hilbertTransform()] which must be called first,
#'   [instantaneousPhase()] for the companion phase extraction,
#'   [waveletTransform()] for alternative time-frequency decomposition.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
#'   samplingRate = 100
#' )
#'
#' # First compute Hilbert transform, then extract amplitude
#' pe <- hilbertTransform(pe)
#' pe <- instantaneousAmplitude(pe)
instantaneousAmplitude <- function(x, assay_name = "analytic",
                                    output_assay = "amplitude") {
  stopifnot(inherits(x, "PhysioExperiment"))

  data <- SummarizedExperiment::assay(x, assay_name)
  amplitude <- Mod(data)

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- amplitude
  SummarizedExperiment::assays(x) <- assays

  .recordProv(x, input_assay = assay_name, output_assay = output_assay,
              .package = "PhysioAnalysis")
}

#' Extract instantaneous phase
#'
#' Extracts the instantaneous phase from the analytic signal.
#'
#' @param x A PhysioExperiment object with analytic signal.
#' @param assay_name Name of the analytic signal assay.
#' @param output_assay Name for the output assay.
#' @return A `PhysioExperiment` object with an additional assay (default
#'   `"phase"`) containing instantaneous phase values in radians
#'   (range \eqn{[-\pi, \pi]}).
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [hilbertTransform()] which must be called first,
#'   [instantaneousAmplitude()] for the companion amplitude extraction,
#'   [plv()] for phase-based connectivity analysis.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
#'   samplingRate = 100
#' )
#'
#' # First compute Hilbert transform, then extract phase
#' pe <- hilbertTransform(pe)
#' pe <- instantaneousPhase(pe)
instantaneousPhase <- function(x, assay_name = "analytic",
                                output_assay = "phase") {
  stopifnot(inherits(x, "PhysioExperiment"))

  data <- SummarizedExperiment::assay(x, assay_name)
  phase <- Arg(data)

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- phase
  SummarizedExperiment::assays(x) <- assays

  .recordProv(x, input_assay = assay_name, output_assay = output_assay,
              .package = "PhysioAnalysis")
}

#' Plot spectrogram
#'
#' Creates a visualization of a spectrogram result.
#'
#' @param spec Spectrogram result from spectrogram().
#' @param freq_range Optional frequency range to display.
#' @param log_power If TRUE, displays log power (dB scale).
#' @return A `ggplot` object showing the spectrogram as a filled raster plot
#'   with time on the x-axis and frequency on the y-axis.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [spectrogram()] to compute the spectrogram data, [plotPSD()] for
#'   power spectral density plots, [plotSignal()] for time-domain visualization.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
#'   samplingRate = 250
#' )
#'
#' # Compute spectrogram
#' spec <- spectrogram(pe, channel = 1)
#'
#' # Plot with frequency range filter
#' plotSpectrogram(spec, freq_range = c(1, 50))
plotSpectrogram <- function(spec, freq_range = NULL, log_power = TRUE) {
  power <- spec$power
  freqs <- spec$frequencies
  times <- spec$times

  if (!is.null(freq_range)) {
    freq_idx <- which(freqs >= freq_range[1] & freqs <= freq_range[2])
    power <- power[freq_idx, ]
    freqs <- freqs[freq_idx]
  }

  if (log_power) {
    power <- 10 * log10(power + 1e-10)
  }

  # Create data frame for ggplot
  plot_df <- expand.grid(time = times, frequency = freqs)
  plot_df$power <- as.vector(t(power))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = frequency, fill = power)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::labs(x = "Time (s)", y = "Frequency (Hz)",
                  fill = if (log_power) "Power (dB)" else "Power",
                  title = "Spectrogram")

  p
}


#' Autocorrelation function (ACF)
#'
#' Computes the sample autocorrelation function of a one-dimensional signal at
#' lags \code{0, 1, \dots, lag_max} using the biased estimator (dividing each
#' autocovariance by \eqn{N}), with the mean removed by default -- the same
#' definition as \code{stats::acf()} and \code{statsmodels.tsa.acf(adjusted =
#' FALSE)}. The autocorrelation measures how similar a signal is to a
#' time-shifted copy of itself, and is the basic tool for detecting periodicity,
#' rhythmicity, and the decorrelation time of a physiological signal.
#'
#' @param x A numeric vector (the time series).
#' @param lag_max Maximum lag (default 30); the returned vector has
#'   \code{lag_max + 1} values (lags 0 through \code{lag_max}).
#' @param demean If \code{TRUE} (default), subtract the mean before computing the
#'   autocovariance (the standard definition).
#' @return A numeric vector of autocorrelation coefficients at lags
#'   \code{0:lag_max}; the lag-0 value is always 1.
#'
#' @references Box, G.E.P., Jenkins, G.M. (1976). Time Series Analysis:
#'   Forecasting and Control. Holden-Day.
#' @seealso \code{\link[stats]{acf}}
#' @export
#' @examples
#' set.seed(1)
#' autocorrelation(sin(seq(0, 20 * pi, length.out = 500)), lag_max = 20)
autocorrelation <- function(x, lag_max = 30L, demean = TRUE) {
  x <- as.numeric(x)
  n <- length(x)
  lag_max <- as.integer(lag_max)
  if (is.na(lag_max) || lag_max < 1L) stop("`lag_max` must be an integer >= 1.", call. = FALSE)
  if (lag_max >= n) stop("`lag_max` must be smaller than the signal length.", call. = FALSE)
  xc <- if (isTRUE(demean)) x - mean(x) else x
  c0 <- sum(xc^2)
  if (c0 <= 0) stop("`x` has zero variance; the ACF is undefined.", call. = FALSE)
  vapply(0:lag_max,
         function(k) sum(xc[1:(n - k)] * xc[(k + 1L):n]) / c0,
         numeric(1))                              # biased, demeaned ACF (= stats::acf)
}


#' Partial autocorrelation function (PACF)
#'
#' Computes the sample partial autocorrelation function of a one-dimensional
#' signal at lags \code{0, 1, \dots, lag_max}. The partial autocorrelation at lag
#' \eqn{k} is the correlation between \eqn{x_t} and \eqn{x_{t+k}} with the linear
#' effect of the intervening lags removed, obtained here by the Durbin-Levinson
#' recursion on the biased, demeaned autocorrelation sequence (via
#' \code{\link{autocorrelation}}) -- the same definition as \code{stats::pacf}
#' and \code{statsmodels.tsa.pacf(method = "ldb")}. Where the ACF \emph{decays},
#' the PACF \emph{cuts off} after the autoregressive order, so it is the standard
#' tool for choosing AR model order.
#'
#' @param x A numeric vector (the time series).
#' @param lag_max Maximum lag (default 20); the returned vector has
#'   \code{lag_max + 1} values (lags 0 through \code{lag_max}).
#' @return A numeric vector of partial autocorrelation coefficients at lags
#'   \code{0:lag_max}; the lag-0 value is 1 by convention and the lag-1 value
#'   equals the lag-1 autocorrelation.
#'
#' @references Box, G.E.P., Jenkins, G.M. (1976). Time Series Analysis:
#'   Forecasting and Control. Holden-Day.
#' @seealso \code{\link{autocorrelation}}, \code{\link[stats]{pacf}}
#' @export
#' @examples
#' set.seed(1)
#' partialAutocorrelation(as.numeric(arima.sim(list(ar = 0.6), 500)), lag_max = 10)
partialAutocorrelation <- function(x, lag_max = 20L) {
  x <- as.numeric(x)
  lag_max <- as.integer(lag_max)
  if (is.na(lag_max) || lag_max < 1L) stop("`lag_max` must be an integer >= 1.", call. = FALSE)
  r <- autocorrelation(x, lag_max = lag_max, demean = TRUE)   # biased demeaned ACF, r[1] = lag 0 = 1
  pacf <- numeric(lag_max)
  a <- numeric(0); v <- r[1]
  for (k in seq_len(lag_max)) {                               # Durbin-Levinson recursion
    if (k == 1L) {
      kk <- r[2] / r[1]
    } else {
      kk <- (r[k + 1L] - sum(a * r[k:2])) / v
      a <- c(a - kk * rev(a), kk)
    }
    if (k == 1L) a <- kk
    pacf[k] <- kk
    v <- v * (1 - kk^2)
  }
  c(1, pacf)                                                  # lags 0..lag_max (lag-0 = 1, matches statsmodels)
}


#' Autoregressive model fit by Yule-Walker (AR(p))
#'
#' Fits an autoregressive model of order \code{order} to a one-dimensional signal
#' by the Yule-Walker method: the AR coefficients solve the Yule-Walker equations
#' formed from the signal's biased demeaned autocovariance (obtained via the
#' Durbin-Levinson recursion on \code{\link{autocorrelation}}), and the innovation
#' variance is the residual (one-step prediction) variance. This is the same
#' definition as \code{statsmodels.regression.linear_model.yule_walker(method =
#' "mle")}. AR modelling underlies forecasting, parametric (AR) spectra, and --
#' with the order chosen from the \code{\link{partialAutocorrelation}} cut-off --
#' compact descriptions of oscillatory dynamics.
#'
#' @param x A numeric vector (the time series).
#' @param order Autoregressive order \eqn{p} (default 4).
#' @return A list with \code{ar} (the length-\code{order} AR coefficient vector
#'   \eqn{\phi_1, \dots, \phi_p}), \code{var_pred} (the innovation / one-step
#'   prediction variance), and \code{order}.
#'
#' @references Box, G.E.P., Jenkins, G.M. (1976). Time Series Analysis:
#'   Forecasting and Control. Holden-Day.
#' @seealso \code{\link{autocorrelation}}, \code{\link{partialAutocorrelation}},
#'   \code{\link[stats]{ar.yw}}
#' @export
#' @examples
#' set.seed(1)
#' arYuleWalker(as.numeric(arima.sim(list(ar = c(0.5, -0.3)), 500)), order = 2)
arYuleWalker <- function(x, order = 4L) {
  x <- as.numeric(x)
  order <- as.integer(order)
  if (is.na(order) || order < 1L) stop("`order` must be an integer >= 1.", call. = FALSE)
  r <- autocorrelation(x, lag_max = order, demean = TRUE)     # biased demeaned ACF, r[1] = 1
  a <- numeric(0); v <- r[1]
  for (k in seq_len(order)) {                                 # Durbin-Levinson: AR coeffs at order p
    if (k == 1L) {
      kk <- r[2] / r[1]
    } else {
      kk <- (r[k + 1L] - sum(a * r[k:2])) / v
      a <- c(a - kk * rev(a), kk)
    }
    if (k == 1L) a <- kk
    v <- v * (1 - kk^2)
  }
  gamma0 <- sum((x - mean(x))^2) / length(x)                  # biased variance (MLE convention)
  list(ar = a, var_pred = gamma0 * v, order = order)
}


#' Ljung-Box test for autocorrelation (white-noise / portmanteau test)
#'
#' Applies the Ljung-Box portmanteau test (Ljung & Box, 1978) to a
#' one-dimensional signal: the statistic
#' \eqn{Q = N(N+2)\sum_{k=1}^{h} \hat\rho_k^2 / (N-k)} aggregates the first
#' \code{lag} squared autocorrelations (from \code{\link{autocorrelation}}), and
#' under the null hypothesis that the series is white noise \eqn{Q} follows a
#' chi-squared distribution with \code{lag} degrees of freedom. It is the
#' standard formal test for the presence of autocorrelation (and, applied to
#' model residuals, for model adequacy). Same definition as
#' \code{statsmodels.stats.diagnostic.acorr_ljungbox} and \code{stats::Box.test(type
#' = "Ljung-Box")}.
#'
#' @param x A numeric vector (the time series).
#' @param lag Number of lags \eqn{h} to include (default 10), and the degrees of
#'   freedom of the reference chi-squared.
#' @return A list with \code{statistic} (the Ljung-Box \eqn{Q}), \code{p_value}
#'   (the chi-squared upper-tail probability), and \code{df} (\code{= lag}).
#'
#' @references Ljung, G.M., Box, G.E.P. (1978). On a measure of lack of fit in
#'   time series models. \emph{Biometrika}, 65(2), 297-303.
#' @seealso \code{\link{autocorrelation}}, \code{\link[stats]{Box.test}}
#' @export
#' @examples
#' set.seed(1)
#' ljungBoxTest(rnorm(500), lag = 10)$p_value   # white noise -> large p
ljungBoxTest <- function(x, lag = 10L) {
  x <- as.numeric(x)
  lag <- as.integer(lag)
  if (is.na(lag) || lag < 1L) stop("`lag` must be an integer >= 1.", call. = FALSE)
  n <- length(x)
  if (lag >= n) stop("`lag` must be smaller than the signal length.", call. = FALSE)
  rho <- autocorrelation(x, lag_max = lag, demean = TRUE)[-1]   # ACF at lags 1..lag
  q <- n * (n + 2) * sum(rho^2 / (n - seq_len(lag)))            # Ljung-Box statistic
  list(statistic = q,
       p_value   = stats::pchisq(q, df = lag, lower.tail = FALSE),
       df        = lag)
}


#' Augmented Dickey-Fuller (ADF) unit-root test statistic
#'
#' Computes the Augmented Dickey-Fuller test statistic for a unit root
#' (non-stationarity) in a one-dimensional signal. The test regresses the first
#' difference on a constant, the lagged level, and \code{lag} lagged first
#' differences, \eqn{\Delta y_t = \alpha + \beta y_{t-1} + \sum_{i=1}^{L}
#' \gamma_i \Delta y_{t-i} + \varepsilon_t}, and returns the t-statistic on
#' \eqn{\beta}. A large negative statistic (below the Dickey-Fuller critical
#' value) rejects the unit-root null in favour of stationarity. The statistic is
#' the same deterministic OLS quantity as
#' \code{statsmodels.tsa.stattools.adfuller(..., autolag = None)}.
#'
#' Note: only the \emph{statistic} is computed here (it is an exact regression
#' quantity). The critical values and p-value require the Dickey-Fuller /
#' MacKinnon reference tables (an approximation), so they are not returned;
#' compare the statistic against the tabulated critical value for the decision.
#'
#' @param x A numeric vector (the time series).
#' @param lag Number of augmenting lagged differences \eqn{L} (default 1).
#' @return A list with \code{statistic} (the ADF t-statistic on the lagged
#'   level), \code{n} (the regression sample size), and \code{lag}.
#'
#' @references Dickey, D.A., Fuller, W.A. (1979). Distribution of the estimators
#'   for autoregressive time series with a unit root. \emph{JASA}, 74, 427-431.
#' @seealso \code{\link{autocorrelation}}, \code{\link{ljungBoxTest}}
#' @export
#' @examples
#' set.seed(1)
#' adfTest(cumsum(rnorm(500)), lag = 1)$statistic   # random walk -> not far below 0
adfTest <- function(x, lag = 1L) {
  x <- as.numeric(x)
  lag <- as.integer(lag)
  if (is.na(lag) || lag < 0L) stop("`lag` must be an integer >= 0.", call. = FALSE)
  dy <- diff(x); n <- length(dy)
  if (n <= lag + 2L) stop("Signal too short for the requested `lag`.", call. = FALSE)
  y_lag <- x[(lag + 1L):(length(x) - 1L)]          # y_{t-1}
  y <- dy[(lag + 1L):n]                            # delta y_t
  xmat <- cbind(1, y_lag)
  for (i in seq_len(lag)) xmat <- cbind(xmat, dy[(lag + 1L - i):(n - i)])   # lagged differences
  fit <- stats::lm.fit(xmat, y)
  b <- fit$coefficients[2]
  res <- fit$residuals
  s2 <- sum(res^2) / (length(y) - ncol(xmat))      # residual variance
  se_b <- sqrt(s2 * solve(crossprod(xmat))[2, 2])  # standard error of the lagged-level coefficient
  list(statistic = unname(b / se_b), n = length(y), lag = lag)
}


#' KPSS test for stationarity
#'
#' Computes the KPSS statistic (Kwiatkowski, Phillips, Schmidt & Shin, 1992) for
#' level stationarity of a one-dimensional signal. The residuals from regressing
#' the signal on a constant are cumulatively summed, and the statistic is
#' \eqn{\eta = N^{-2}\sum_t S_t^2 / \hat\lambda^2}, where \eqn{S_t} are the
#' partial sums and \eqn{\hat\lambda^2} is a Bartlett (Newey-West) long-run
#' variance estimate with truncation \code{lag}. Unlike the Augmented
#' Dickey-Fuller test, the KPSS \emph{null} is stationarity: a \emph{small}
#' statistic (below the critical value) means the series is consistent with
#' stationarity. The statistic is the same deterministic quantity as
#' \code{statsmodels.tsa.stattools.kpss(regression = "c", nlags = lag)}.
#'
#' Note: only the \emph{statistic} is computed here (an exact quantity). The KPSS
#' critical values and p-value come from the reference tables (an approximation)
#' and are not returned.
#'
#' @param x A numeric vector (the time series).
#' @param lag Bartlett truncation lag \eqn{L} for the long-run variance (default
#'   10).
#' @return A list with \code{statistic} (the KPSS \eqn{\eta}), \code{n} (the
#'   signal length), and \code{lag}.
#'
#' @references Kwiatkowski, D., Phillips, P.C.B., Schmidt, P., Shin, Y. (1992).
#'   Testing the null hypothesis of stationarity against the alternative of a
#'   unit root. \emph{Journal of Econometrics}, 54(1-3), 159-178.
#' @seealso \code{\link{adfTest}}, \code{\link{autocorrelation}}
#' @export
#' @examples
#' set.seed(1)
#' kpssTest(rnorm(500), lag = 10)$statistic   # stationary -> small statistic
kpssTest <- function(x, lag = 10L) {
  x <- as.numeric(x)
  lag <- as.integer(lag)
  if (is.na(lag) || lag < 0L) stop("`lag` must be an integer >= 0.", call. = FALSE)
  n <- length(x)
  if (n <= lag + 1L) stop("Signal too short for the requested `lag`.", call. = FALSE)
  e <- x - mean(x)                                   # residuals from regression on a constant
  s <- cumsum(e)                                     # partial sums
  lrv <- sum(e^2) / n                                # gamma_0
  for (k in seq_len(lag)) {                          # Bartlett (Newey-West) long-run variance
    gs <- sum(e[1:(n - k)] * e[(k + 1L):n]) / n
    lrv <- lrv + 2 * (1 - k / (lag + 1)) * gs
  }
  list(statistic = sum(s^2) / (n^2 * lrv), n = n, lag = lag)
}


#' Autoregressive (parametric) spectral density
#'
#' Computes the parametric AR spectral density of a one-dimensional signal from
#' an autoregressive model of order \code{order} fitted by Yule-Walker (via
#' \code{\link{arYuleWalker}}): \eqn{S(f) = \sigma^2 / |1 - \sum_k \phi_k
#' e^{-i 2\pi f k}|^2}, evaluated on \code{n_freq} normalized frequencies in
#' \eqn{[0, 0.5]}. Unlike the periodogram / Welch estimate, the AR spectrum is
#' smooth and low-variance, with sharp peaks at resonant frequencies, using only
#' \code{order} parameters. It matches \code{stats::spec.ar} (including its
#' small-sample innovation-variance convention \eqn{N/(N-p-1)}).
#'
#' @param x A numeric vector (the time series).
#' @param order Autoregressive order \eqn{p} (default 8).
#' @param n_freq Number of frequencies in \eqn{[0, 0.5]} (default 500).
#' @return A list with \code{freq} (normalized frequencies, cycles/sample),
#'   \code{spec} (the AR spectral density), and \code{order}.
#'
#' @references Percival, D.B., Walden, A.T. (1993). Spectral Analysis for
#'   Physical Applications. Cambridge University Press.
#' @seealso \code{\link{arYuleWalker}}, \code{\link[stats]{spec.ar}}
#' @export
#' @examples
#' set.seed(1)
#' s <- arSpectrum(as.numeric(arima.sim(list(ar = c(0.6, -0.3)), 500)), order = 2)
#' s$freq[which.max(s$spec)]   # peak (normalized) frequency
arSpectrum <- function(x, order = 8L, n_freq = 500L) {
  x <- as.numeric(x)
  order <- as.integer(order); n_freq <- as.integer(n_freq)
  if (is.na(order) || order < 1L) stop("`order` must be an integer >= 1.", call. = FALSE)
  if (is.na(n_freq) || n_freq < 2L) stop("`n_freq` must be an integer >= 2.", call. = FALSE)
  fit <- arYuleWalker(x, order = order)
  n <- length(x)
  var_df <- fit$var_pred * n / (n - order - 1L)          # spec.ar / ar.yw small-sample variance
  freq <- seq(0, 0.5, length.out = n_freq)
  z <- exp(-1i * 2 * pi * outer(seq_len(order), freq))
  denom <- 1 - colSums(fit$ar * z)
  list(freq = freq, spec = var_df / (Mod(denom)^2), order = order)
}


#' Signal distribution moments (mean, SD, skewness, kurtosis)
#'
#' Computes the first four moments of a one-dimensional signal's amplitude
#' distribution: the mean, the (population) standard deviation, the skewness, and
#' the excess kurtosis. Skewness measures asymmetry (negative = a longer left
#' tail); excess kurtosis measures tailedness relative to a Gaussian (0 =
#' Gaussian, positive = heavier tails / more peaked). These are the standard
#' shape descriptors used, for example, to flag EEG artifacts (high kurtosis) or
#' characterize amplitude asymmetry. Definitions match
#' \code{scipy.stats.skew} (biased) and \code{scipy.stats.kurtosis}
#' (Fisher / excess, biased), with the population standard deviation
#' (\code{ddof = 0}).
#'
#' @param x A numeric vector (the time series).
#' @return A list with \code{mean}, \code{sd} (population), \code{skewness}, and
#'   \code{kurtosis} (excess).
#'
#' @references Joanes, D.N., Gill, C.A. (1998). Comparing measures of sample
#'   skewness and kurtosis. \emph{The Statistician}, 47(1), 183-189.
#' @seealso \code{\link[stats]{sd}}
#' @export
#' @examples
#' set.seed(1)
#' signalMoments(rnorm(1000))   # ~ 0 skewness, ~ 0 excess kurtosis for Gaussian
signalMoments <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 2L) stop("`x` must have at least 2 values.", call. = FALSE)
  mu <- mean(x)
  cen <- x - mu
  s2 <- sum(cen^2) / n                               # population variance (ddof = 0)
  sdp <- sqrt(s2)
  if (sdp == 0) stop("`x` has zero variance; skewness/kurtosis are undefined.", call. = FALSE)
  list(mean = mu,
       sd = sdp,
       skewness = (sum(cen^3) / n) / sdp^3,          # scipy.stats.skew (biased)
       kurtosis = (sum(cen^4) / n) / sdp^4 - 3)      # scipy.stats.kurtosis (Fisher/excess, biased)
}


#' Jarque-Bera test for normality
#'
#' Applies the Jarque-Bera goodness-of-fit test for normality (Jarque & Bera,
#' 1980) to a one-dimensional signal: the statistic
#' \eqn{JB = \frac{N}{6}\left(S^2 + \frac{K^2}{4}\right)}, where \eqn{S} is the
#' skewness and \eqn{K} the excess kurtosis (from \code{\link{signalMoments}}),
#' is compared to a chi-squared distribution with 2 degrees of freedom. It is the
#' formal, moment-based counterpart of the descriptive near-Gaussian check: a
#' large statistic rejects normality. Same definition as
#' \code{scipy.stats.jarque_bera}.
#'
#' @param x A numeric vector (the time series).
#' @return A list with \code{statistic} (the JB statistic), \code{p_value} (the
#'   chi-squared upper-tail probability), and \code{df} (2).
#'
#' @references Jarque, C.M., Bera, A.K. (1980). Efficient tests for normality,
#'   homoscedasticity and serial independence of regression residuals.
#'   \emph{Economics Letters}, 6(3), 255-259.
#' @seealso \code{\link{signalMoments}}
#' @export
#' @examples
#' set.seed(1)
#' jarqueBeraTest(rnorm(1000))$p_value   # Gaussian -> large p (fail to reject)
jarqueBeraTest <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 2L) stop("`x` must have at least 2 values.", call. = FALSE)
  m <- signalMoments(x)                              # skewness + excess kurtosis
  jb <- n / 6 * (m$skewness^2 + (m$kurtosis^2) / 4)  # Jarque-Bera statistic
  list(statistic = jb,
       p_value   = stats::pchisq(jb, df = 2, lower.tail = FALSE),
       df        = 2L)
}

#' Median Absolute Deviation (robust dispersion)
#'
#' Computes the median absolute deviation (MAD) of a one-dimensional signal --
#' \eqn{\mathrm{MAD} = c \cdot \mathrm{median}(|x - \mathrm{median}(x)|)} -- a
#' robust measure of dispersion that, unlike the standard deviation, is not
#' inflated by a small number of extreme values (spikes, movement artifacts).
#' With the default \code{constant = 1} it returns the raw MAD; set
#' \code{constant = 1 / stats::qnorm(0.75)} (approximately 1.4826) to obtain the
#' normal-consistent robust estimate of the standard deviation (the scale MAD and
#' SD share for Gaussian data).
#'
#' The robust dispersion companion of \code{\link{signalMoments}} (whose SD is the
#' classical, non-robust dispersion). Comparing the normal-consistent MAD with the
#' SD is a quick artifact check: they agree for clean, near-Gaussian data and
#' diverge (MAD-scale below SD) when heavy-tailed artifacts inflate the SD.
#'
#' @param x A numeric vector (the time series).
#' @param constant Scale factor (default \code{1}, the raw MAD). Use
#'   \code{1 / stats::qnorm(0.75)} for the normal-consistent robust SD estimate.
#' @return A single numeric value, the (scaled) median absolute deviation.
#'
#' @seealso \code{\link{signalMoments}} for the classical moments (mean, SD,
#'   skewness, kurtosis).
#'
#' @export
#' @examples
#' medianAbsDev(c(1, 2, 3, 4, 100))                         # robust to the outlier
#' set.seed(1); medianAbsDev(rnorm(1000), 1 / stats::qnorm(0.75))  # ~1 (robust SD)
medianAbsDev <- function(x, constant = 1) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 1L) stop("`x` must have at least one value.", call. = FALSE)
  stopifnot(is.numeric(constant), length(constant) == 1L, is.finite(constant))
  constant * stats::median(abs(x - stats::median(x)))
}

#' Interquartile Range (robust spread)
#'
#' Computes the interquartile range (IQR) of a one-dimensional signal -- the
#' spread of the middle 50%, \eqn{Q_3 - Q_1} -- a robust, quantile-based measure
#' of dispersion that ignores the tails entirely. The quartiles use R's default
#' type-7 quantile (the same linear interpolation as
#' \code{scipy.stats.iqr} and \code{numpy.percentile}), so the result matches both
#' bit-for-bit.
#'
#' A quantile-based robust spread that complements \code{\link{medianAbsDev}} (a
#' deviation-based robust scale) and \code{\link{signalMoments}} (the classical
#' SD). For approximately Gaussian data, \eqn{\mathrm{IQR} / (2\,\Phi^{-1}(0.75))}
#' (\eqn{\approx \mathrm{IQR}/1.349}) is a robust estimate of the SD; a value well
#' below the SD flags a skewed or heavy-tailed distribution.
#'
#' @param x A numeric vector (the time series).
#' @param type Quantile algorithm passed to \code{\link[stats]{quantile}}
#'   (default \code{7}, matching base R, scipy and numpy).
#' @return A single numeric value, the interquartile range.
#'
#' @seealso \code{\link{medianAbsDev}} for the deviation-based robust scale,
#'   \code{\link{signalMoments}} for the classical moments.
#'
#' @export
#' @examples
#' interquartileRange(c(1, 2, 3, 4, 100))   # robust to the outlier
#' interquartileRange(rnorm(1000)) / 1.349  # ~1 (robust SD estimate)
interquartileRange <- function(x, type = 7L) {
  x <- as.numeric(x)
  if (length(x) < 1L) stop("`x` must have at least one value.", call. = FALSE)
  q <- stats::quantile(x, probs = c(0.25, 0.75), type = type, names = FALSE)
  unname(q[2L] - q[1L])
}

#' Trimmed Mean (robust location)
#'
#' Computes the trimmed mean of a one-dimensional signal -- the arithmetic mean
#' after discarding the most extreme \code{trim} fraction of the values from each
#' end. A robust measure of central tendency that, unlike the ordinary mean, is
#' not pulled by a heavy tail or a few outliers, while (unlike the median) it
#' still uses all of the retained data. It matches base R
#' \code{mean(x, trim = ...)} and \code{scipy.stats.trim_mean} bit-for-bit, both
#' of which discard \eqn{\lfloor n\,\mathrm{trim}\rfloor} values from each end.
#'
#' A robust location estimator that completes the robust toolkit alongside
#' \code{\link{medianAbsDev}} and \code{\link{interquartileRange}} (robust
#' dispersion). On a skewed distribution the trimmed mean lies between the
#' ordinary mean and the median, and moves toward the median as \code{trim}
#' increases -- a measure of how much the tail inflates the mean.
#'
#' @param x A numeric vector (the time series).
#' @param trim Fraction (0 to 0.5) of values trimmed from EACH end (default 0.1).
#' @return A single numeric value, the trimmed mean.
#'
#' @seealso \code{\link{signalMoments}} for the ordinary mean and SD,
#'   \code{\link{medianAbsDev}} and \code{\link{interquartileRange}} for robust
#'   dispersion.
#'
#' @export
#' @examples
#' trimmedMean(c(1, 2, 3, 4, 100), trim = 0.2)  # 3 (robust); the mean is 22
#' set.seed(1); trimmedMean(rexp(1000), trim = 0.2)  # below the tail-inflated mean
trimmedMean <- function(x, trim = 0.1) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 1L) stop("`x` must have at least one value.", call. = FALSE)
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 0.5)
    stop("`trim` must be a single number in [0, 0.5).", call. = FALSE)
  k <- as.integer(floor(n * trim))
  xs <- sort(x)
  mean(xs[(k + 1L):(n - k)])
}

#' Kolmogorov-Smirnov Test for Normality
#'
#' Tests whether a one-dimensional signal is normally distributed by comparing
#' its empirical cumulative distribution function (ECDF) to that of a normal
#' distribution fitted to the sample (mean and SD). The statistic is the
#' Kolmogorov-Smirnov distance \eqn{D = \sup_x |F_n(x) - \Phi((x-\mu)/\sigma)|},
#' the largest vertical gap between the two CDFs; the two-sided asymptotic
#' (Kolmogorov) p-value follows. The \eqn{D} statistic reproduces
#' \code{stats::ks.test} and \code{scipy.stats.kstest} bit-for-bit.
#'
#' A distribution-shape (empirical-CDF) normality test, complementary to the
#' moment-based \code{\link{jarqueBeraTest}}: KS responds to the overall shape of
#' the distribution, Jarque-Bera to its skewness and kurtosis (the tails). The two
#' can disagree -- on a long, mildly non-Gaussian signal Jarque-Bera's power often
#' rejects while KS does not.
#'
#' @param x A numeric vector (the time series).
#' @return A list with \code{statistic} (the KS distance \eqn{D}),
#'   \code{p_value} (two-sided asymptotic), and \code{n} (sample size).
#'
#' @section Caveat:
#' Because the normal parameters are estimated from the same sample, the naive
#' p-value is anti-conservative (a Lilliefors correction gives an exact test); the
#' \eqn{D} statistic itself is exact and is what this function certifies.
#'
#' @seealso \code{\link{jarqueBeraTest}} for the moment-based normality test,
#'   \code{\link{signalMoments}} for the skewness and kurtosis.
#'
#' @export
#' @examples
#' set.seed(1); ksNormalityTest(rnorm(500))$p_value    # Gaussian -> large p
#' set.seed(1); ksNormalityTest(rexp(500))$statistic   # skewed -> large D
ksNormalityTest <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 2L) stop("`x` must have at least 2 values.", call. = FALSE)
  mu <- mean(x); sg <- stats::sd(x)
  if (!is.finite(sg) || sg <= 0) stop("`x` has (near) zero variance.", call. = FALSE)
  xs <- sort(x)
  F0 <- stats::pnorm(xs, mu, sg)
  D <- max(max((1:n) / n - F0), max(F0 - (0:(n - 1L)) / n))
  tt <- sqrt(n) * D
  k <- 1:100
  p <- 2 * sum((-1)^(k - 1L) * exp(-2 * k^2 * tt^2))   # asymptotic two-sided Kolmogorov
  list(statistic = D, p_value = min(max(p, 0), 1), n = n)
}

#' Paired-Samples t-Test
#'
#' Computes the two-sided paired-samples t-test comparing two matched vectors --
#' the one-sample t-test of the within-pair differences \eqn{d = x - y}:
#' \eqn{t = \bar{d} / (s_d / \sqrt{n})} on \eqn{n - 1} degrees of freedom. The
#' standard test for before/after or condition-A/condition-B measurements on the
#' same subjects. Reproduces \code{stats::t.test(x, y, paired = TRUE)} and
#' \code{scipy.stats.ttest_rel} bit-for-bit.
#'
#' @param x A numeric vector, OR (when \code{y} is \code{NULL}) a two-column
#'   matrix whose columns are the paired samples.
#' @param y A numeric vector matched to \code{x}, or \code{NULL} if \code{x} is a
#'   two-column matrix.
#' @return A list with \code{statistic} (the t-statistic of \eqn{x - y}),
#'   \code{p_value} (two-sided) and \code{df}.
#'
#' @seealso \code{\link{signalMoments}}, \code{\link{jarqueBeraTest}}.
#'
#' @export
#' @examples
#' set.seed(1); pairedTTest(rnorm(20, 1), rnorm(20))$statistic  # x shifted up -> positive t
pairedTTest <- function(x, y = NULL) {
  if (is.null(y)) {
    x <- as.matrix(x)
    if (ncol(x) != 2L) stop("`x` must have 2 columns when `y` is NULL.", call. = FALSE)
    y <- as.numeric(x[, 2L]); x <- as.numeric(x[, 1L])
  } else {
    x <- as.numeric(x); y <- as.numeric(y)
  }
  if (length(x) != length(y)) stop("`x` and `y` must have the same length.", call. = FALSE)
  d <- x - y; n <- length(d)
  if (n < 2L) stop("need at least 2 paired observations.", call. = FALSE)
  sdd <- stats::sd(d)
  if (!is.finite(sdd) || sdd <= 0) stop("zero within-pair variance; t is undefined.", call. = FALSE)
  t <- mean(d) / (sdd / sqrt(n))
  df <- n - 1L
  list(statistic = t, p_value = 2 * stats::pt(-abs(t), df), df = df)
}

#' Correlation Test
#'
#' Tests the association between two matched numeric vectors. Returns the
#' correlation coefficient, its two-sided p-value (from the exact t-test on the
#' Pearson/Spearman correlation) and the degrees of freedom. For
#' \code{method = "pearson"} it reproduces \code{stats::cor.test} and
#' \code{scipy.stats.pearsonr} bit-for-bit.
#'
#' @param x A numeric vector, OR (when \code{y} is \code{NULL}) a matrix whose
#'   first two columns are the paired variables.
#' @param y A numeric vector matched to \code{x}, or \code{NULL} if \code{x} is a
#'   matrix.
#' @param method Correlation method: \code{"pearson"} (default) or
#'   \code{"spearman"} (rank correlation; the t-based p is the large-sample
#'   approximation).
#' @return A list with \code{statistic} (the correlation coefficient),
#'   \code{p_value} (two-sided) and \code{df} (\eqn{n - 2}).
#'
#' @seealso \code{\link{pairedTTest}}, \code{\link{signalMoments}}.
#'
#' @export
#' @examples
#' set.seed(1); x <- rnorm(50); correlationTest(x, x + rnorm(50))$statistic  # ~0.7
correlationTest <- function(x, y = NULL, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (is.null(y)) {
    x <- as.matrix(x)
    if (ncol(x) < 2L) stop("`x` needs >= 2 columns when `y` is NULL.", call. = FALSE)
    y <- as.numeric(x[, 2L]); x <- as.numeric(x[, 1L])
  } else {
    x <- as.numeric(x); y <- as.numeric(y)
  }
  if (length(x) != length(y)) stop("`x` and `y` must have the same length.", call. = FALSE)
  n <- length(x)
  if (n < 3L) stop("need at least 3 paired observations.", call. = FALSE)
  r <- stats::cor(x, y, method = method)
  df <- n - 2L
  if (abs(r) >= 1) {
    p <- 0
  } else {
    tstat <- r * sqrt(df / (1 - r^2))
    p <- 2 * stats::pt(-abs(tstat), df)
  }
  list(statistic = r, p_value = p, df = df)
}

#' Two-Sample t-Test
#'
#' Computes the two-sided independent two-sample t-test comparing two groups. By
#' default it is Welch's t-test (unequal variances, Welch-Satterthwaite df) --
#' matching \code{stats::t.test} and \code{scipy.stats.ttest_ind(equal_var =
#' FALSE)}; set \code{var_equal = TRUE} for Student's pooled-variance test. The
#' standard test for a difference between two independent groups (e.g. older vs
#' younger, patients vs controls). Reproduces both references bit-for-bit.
#'
#' @param x A numeric vector (group 1), OR (when \code{y} is \code{NULL}) a list
#'   of two numeric vectors.
#' @param y A numeric vector (group 2), or \code{NULL} if \code{x} is a two-element
#'   list.
#' @param var_equal Logical; \code{FALSE} (default) for Welch's t-test, \code{TRUE}
#'   for Student's pooled-variance t-test.
#' @return A list with \code{statistic} (the t-statistic of \eqn{x - y}),
#'   \code{p_value} (two-sided) and \code{df}.
#'
#' @seealso \code{\link{pairedTTest}} for matched samples.
#'
#' @export
#' @examples
#' set.seed(1); twoSampleTTest(rnorm(20, 1), rnorm(30))$statistic  # group 1 higher -> positive t
twoSampleTTest <- function(x, y = NULL, var_equal = FALSE) {
  if (is.null(y)) {
    if (!is.list(x) || length(x) != 2L) stop("`x` must be a list of two vectors when `y` is NULL.", call. = FALSE)
    y <- as.numeric(x[[2L]]); x <- as.numeric(x[[1L]])
  } else {
    x <- as.numeric(x); y <- as.numeric(y)
  }
  nx <- length(x); ny <- length(y)
  if (nx < 2L || ny < 2L) stop("each group needs at least 2 observations.", call. = FALSE)
  mx <- mean(x); my <- mean(y); vx <- stats::var(x); vy <- stats::var(y)
  if (var_equal) {
    sp2 <- ((nx - 1L) * vx + (ny - 1L) * vy) / (nx + ny - 2L)
    se <- sqrt(sp2 * (1 / nx + 1 / ny)); df <- nx + ny - 2L
  } else {
    se <- sqrt(vx / nx + vy / ny)
    df <- (vx / nx + vy / ny)^2 / ((vx / nx)^2 / (nx - 1L) + (vy / ny)^2 / (ny - 1L))
  }
  t <- (mx - my) / se
  list(statistic = t, p_value = 2 * stats::pt(-abs(t), df), df = df)
}

#' One-way ANOVA (omnibus F-test) across groups
#'
#' Classic one-way analysis of variance (equal-variance omnibus F-test)
#' comparing the means of `k >= 2` independent groups. The F-statistic and its
#' p-value reproduce base R `stats::oneway.test(var.equal = TRUE)` and
#' `scipy.stats.f_oneway` bit-for-bit. This is the standard many-group
#' generalisation of the independent two-sample t-test
#' (\code{\link{twoSampleTTest}}).
#'
#' @param x Either a list of `k >= 2` numeric group vectors, or the first
#'   group's numeric vector (with further groups passed through `...`).
#' @param ... Additional numeric group vectors when `x` is a single vector.
#' @return A list with `statistic` (the F-ratio), `p_value`, `df1`
#'   (between-groups degrees of freedom, `k - 1`), `df2` (within-groups,
#'   `N - k`), `k` (number of groups) and `n` (total sample size).
#' @seealso \code{\link{twoSampleTTest}}, \code{\link{pairedTTest}}
#' @export
#' @examples
#' oneWayAnova(list(rnorm(20), rnorm(20, 1), rnorm(20, 2)))
oneWayAnova <- function(x, ...) {
  groups <- if (is.list(x) && !is.data.frame(x)) x else c(list(x), list(...))
  groups <- lapply(groups, function(g) { g <- as.numeric(g); g[is.finite(g)] })
  k <- length(groups)
  if (k < 2L) stop("one-way ANOVA needs at least 2 groups.", call. = FALSE)
  ni <- vapply(groups, length, integer(1))
  if (any(ni < 2L)) stop("each group needs at least 2 observations.", call. = FALSE)
  N  <- sum(ni)
  mi <- vapply(groups, mean, numeric(1))
  gm <- sum(vapply(groups, sum, numeric(1))) / N
  ssb <- sum(ni * (mi - gm)^2)                                   # between-groups SS
  ssw <- sum(vapply(groups, function(g) sum((g - mean(g))^2), numeric(1)))  # within-groups SS
  df1 <- k - 1L; df2 <- N - k
  Fv  <- (ssb / df1) / (ssw / df2)
  list(statistic = Fv, p_value = stats::pf(Fv, df1, df2, lower.tail = FALSE),
       df1 = df1, df2 = df2, k = k, n = N)
}

#' Friedman rank-sum test (non-parametric repeated-measures omnibus)
#'
#' The non-parametric repeated-measures omnibus test: compares `k >= 2` matched
#' treatments (columns) measured on the same `N` blocks/subjects (rows), by
#' ranking within each block and testing whether the rank totals differ. It is
#' the within-subject, distribution-free counterpart of the one-way ANOVA
#' (\code{\link{oneWayAnova}}). The tie-corrected chi-square statistic and its
#' p-value reproduce base R `stats::friedman.test` bit-for-bit and
#' `scipy.stats.friedmanchisquare` to machine precision.
#'
#' @param x A blocks-by-treatments numeric matrix (rows = subjects/blocks,
#'   columns = conditions/treatments), or a list of `k` equal-length numeric
#'   vectors (one per treatment, aligned by block). Rows with any missing value
#'   are dropped (complete blocks only).
#' @param ... Additional treatment vectors when `x` is a single vector.
#' @return A list with `statistic` (the tie-corrected Friedman chi-square),
#'   `p_value`, `df` (`k - 1`), `k` (treatments) and `n` (complete blocks).
#' @seealso \code{\link{oneWayAnova}}, \code{\link{pairedTTest}}
#' @export
#' @examples
#' m <- matrix(rnorm(40) + rep(c(0, 0.5, 1, 1.5), each = 10), nrow = 10)
#' friedmanTest(m)
friedmanTest <- function(x, ...) {
  M <- if (is.matrix(x)) x
       else if (is.list(x) && !is.data.frame(x)) do.call(cbind, x)
       else do.call(cbind, c(list(x), list(...)))
  M <- matrix(as.numeric(M), nrow = nrow(M))
  M <- M[stats::complete.cases(M), , drop = FALSE]
  N <- nrow(M); k <- ncol(M)
  if (k < 2L) stop("Friedman test needs at least 2 treatments (columns).", call. = FALSE)
  if (N < 2L) stop("Friedman test needs at least 2 complete blocks (rows).", call. = FALSE)
  R  <- t(apply(M, 1L, rank))                                 # rank within each block
  Rj <- colSums(R)
  ties <- sum(apply(M, 1L, function(v) { tt <- table(v); sum(tt^3 - tt) }))  # tie correction
  Q  <- (12 * sum((Rj - N * (k + 1) / 2)^2)) / (N * k * (k + 1) - ties / (k - 1))
  list(statistic = Q, p_value = stats::pchisq(Q, k - 1L, lower.tail = FALSE),
       df = k - 1L, k = k, n = N)
}

#' Kruskal-Wallis rank-sum test (non-parametric between-groups omnibus)
#'
#' The non-parametric between-groups omnibus test: pools all observations,
#' ranks them, and tests whether the rank totals differ across `k >= 2`
#' independent groups. It is the distribution-free counterpart of the one-way
#' ANOVA (\code{\link{oneWayAnova}}) and, unlike it, makes no normality or
#' equal-variance assumption. The tie-corrected statistic and its p-value
#' reproduce base R `stats::kruskal.test` bit-for-bit (the H is accumulated in
#' the same order) and `scipy.stats.kruskal` to machine precision.
#'
#' @param x Either a list of `k >= 2` numeric group vectors, or the first
#'   group's numeric vector (with further groups passed through `...`).
#' @param ... Additional numeric group vectors when `x` is a single vector.
#' @return A list with `statistic` (the tie-corrected H), `p_value`, `df`
#'   (`k - 1`), `k` (number of groups) and `n` (total sample size).
#' @seealso \code{\link{oneWayAnova}}, \code{\link{friedmanTest}}
#' @export
#' @examples
#' kruskalTest(list(rnorm(20), rnorm(20, 1), rnorm(20, 2)))
kruskalTest <- function(x, ...) {
  groups <- if (is.list(x) && !is.data.frame(x)) x else c(list(x), list(...))
  groups <- lapply(groups, function(g) { g <- as.numeric(g); g[is.finite(g)] })
  k <- length(groups)
  if (k < 2L) stop("Kruskal-Wallis test needs at least 2 groups.", call. = FALSE)
  ni <- vapply(groups, length, integer(1))
  if (any(ni < 1L)) stop("each group needs at least 1 observation.", call. = FALSE)
  x_all <- unlist(groups, use.names = FALSE)
  n <- length(x_all)
  r <- rank(x_all)
  idx <- rep(seq_len(k), ni)
  # accumulate exactly as base R stats::kruskal.test (identical FP order -> bit-for-bit)
  S    <- sum(tapply(r, idx, sum)^2 / tapply(r, idx, length))
  TIES <- table(x_all)
  H    <- ((12 * S / (n * (n + 1)) - 3 * (n + 1)) / (1 - sum(TIES^3 - TIES) / (n^3 - n)))
  list(statistic = H, p_value = stats::pchisq(H, k - 1L, lower.tail = FALSE),
       df = k - 1L, k = k, n = n)
}

#' Simple linear regression (ordinary least squares)
#'
#' Univariate OLS regression of `y` on `x`: the slope, intercept, R-squared,
#' and the two-sided t-test of the slope. The slope, R-squared and p-value
#' reproduce base R `stats::lm` / `summary.lm` and `scipy.stats.linregress`
#' bit-for-bit (to machine precision against `lm`'s QR path).
#'
#' @param x Either the numeric predictor vector (with `y` the response), or a
#'   two-column matrix/data frame (column 1 = `x`, column 2 = `y`), or a list of
#'   two numeric vectors. Incomplete `(x, y)` pairs are dropped.
#' @param y The numeric response vector when `x` is a single predictor vector.
#' @return A list with `slope`, `intercept`, `r_squared`, `statistic` (the
#'   slope t-statistic), `p_value` (two-sided), `df` (`n - 2`) and `n`.
#' @seealso \code{\link{correlationTest}}
#' @export
#' @examples
#' linearRegression(1:20, 2 * (1:20) + rnorm(20))
linearRegression <- function(x, y = NULL) {
  if (is.null(y)) {
    if (is.matrix(x) || is.data.frame(x)) {
      if (ncol(x) < 2L) stop("a single-argument `x` must have >= 2 columns.", call. = FALSE)
      xx <- as.numeric(x[, 1L]); yy <- as.numeric(x[, 2L])
    } else if (is.list(x) && length(x) >= 2L) {
      xx <- as.numeric(x[[1L]]); yy <- as.numeric(x[[2L]])
    } else stop("provide `x` and `y`, or a two-column `x`.", call. = FALSE)
  } else { xx <- as.numeric(x); yy <- as.numeric(y) }
  ok <- is.finite(xx) & is.finite(yy); xx <- xx[ok]; yy <- yy[ok]
  n <- length(xx)
  if (n < 3L) stop("need at least 3 complete (x, y) pairs.", call. = FALSE)
  mx <- mean(xx); my <- mean(yy)
  sxx <- sum((xx - mx)^2); sxy <- sum((xx - mx) * (yy - my)); syy <- sum((yy - my)^2)
  if (sxx == 0) stop("`x` has zero variance.", call. = FALSE)
  slope <- sxy / sxx; intercept <- my - slope * mx
  r2 <- (sxy * sxy) / (sxx * syy)
  df <- n - 2L
  sse <- sum((yy - (intercept + slope * xx))^2)
  se_slope <- sqrt((sse / df) / sxx)
  t <- slope / se_slope
  list(slope = slope, intercept = intercept, r_squared = r2, statistic = t,
       p_value = 2 * stats::pt(-abs(t), df), df = df, n = n)
}

#' Kendall's tau-b rank correlation
#'
#' Kendall's tau-b, the concordance-based rank correlation: the normalised
#' difference between concordant and discordant pairs, tie-corrected. It is the
#' third rank-correlation method alongside Pearson's r
#' (\code{\link{correlationTest}}, method "pearson") and Spearman's rho (method
#' "spearman"), and is the robust choice for ordinal or heavily-tied data and
#' small samples. The tau-b statistic and the normal-approximation p-value
#' reproduce base R `stats::cor.test(method = "kendall")` and
#' `scipy.stats.kendalltau` bit-for-bit.
#'
#' @param x Either the numeric predictor vector (with `y` the response), or a
#'   two-column matrix/data frame (column 1 = `x`, column 2 = `y`), or a list of
#'   two numeric vectors. Incomplete `(x, y)` pairs are dropped.
#' @param y The second numeric vector when `x` is a single vector.
#' @return A list with `statistic` (tau-b), `p_value` (two-sided,
#'   normal-approximation), `z` (the approximation's z-statistic) and `n`.
#' @seealso \code{\link{correlationTest}}
#' @export
#' @examples
#' kendallTau(c(1, 2, 3, 4, 5), c(2, 1, 4, 3, 5))
kendallTau <- function(x, y = NULL) {
  if (is.null(y)) {
    if (is.matrix(x) || is.data.frame(x)) {
      if (ncol(x) < 2L) stop("a single-argument `x` must have >= 2 columns.", call. = FALSE)
      xx <- as.numeric(x[, 1L]); yy <- as.numeric(x[, 2L])
    } else if (is.list(x) && length(x) >= 2L) {
      xx <- as.numeric(x[[1L]]); yy <- as.numeric(x[[2L]])
    } else stop("provide `x` and `y`, or a two-column `x`.", call. = FALSE)
  } else { xx <- as.numeric(x); yy <- as.numeric(y) }
  ok <- is.finite(xx) & is.finite(yy); xx <- xx[ok]; yy <- yy[ok]
  n <- length(xx)
  if (n < 3L) stop("need at least 3 complete (x, y) pairs.", call. = FALSE)
  tau <- stats::cor(xx, yy, method = "kendall")
  # normal-approximation z with tie corrections, exactly as base R cor.test(method="kendall")
  xt <- table(xx[duplicated(xx)]) + 1
  yt <- table(yy[duplicated(yy)]) + 1
  T0 <- n * (n - 1) / 2
  T1 <- sum(xt * (xt - 1)) / 2; T2 <- sum(yt * (yt - 1)) / 2
  S  <- tau * sqrt((T0 - T1) * (T0 - T2))
  v0 <- n * (n - 1) * (2 * n + 5)
  vt <- sum(xt * (xt - 1) * (2 * xt + 5)); vu <- sum(yt * (yt - 1) * (2 * yt + 5))
  v1 <- sum(xt * (xt - 1)) * sum(yt * (yt - 1))
  v2 <- sum(xt * (xt - 1) * (xt - 2)) * sum(yt * (yt - 1) * (yt - 2))
  var_S <- (v0 - vt - vu) / 18 + v1 / (2 * n * (n - 1)) + v2 / (9 * n * (n - 1) * (n - 2))
  z <- S / sqrt(var_S)
  list(statistic = tau, p_value = 2 * stats::pnorm(-abs(z)), z = z, n = n)
}

#' Partial correlation controlling for covariates
#'
#' The correlation between `x` and `y` after linearly removing one or more
#' covariates `z` from both -- the standard tool for asking whether an
#' association survives adjustment for a confound. With `method = "spearman"`
#' the variables are rank-transformed first (a rank / robust partial
#' correlation). The coefficient reproduces base R's residual method and the
#' `ppcor` package to machine precision.
#'
#' @param x Either the first numeric variable (with `y` and `z` given), or a
#'   matrix/data frame whose column 1 is `x`, column 2 is `y` and the remaining
#'   columns are the covariates. Incomplete rows are dropped.
#' @param y The second numeric variable when `x` is a single vector.
#' @param z The covariate(s) to control for: a numeric vector or matrix.
#' @param method "pearson" (default) or "spearman" (rank-based).
#' @return A list with `statistic` (the partial correlation), `p_value`
#'   (two-sided t-test), `df` (`n - 2 - k`), `n` and `k` (number of covariates).
#' @seealso \code{\link{correlationTest}}, \code{\link{kendallTau}}
#' @export
#' @examples
#' x <- rnorm(50); z <- rnorm(50); y <- 0.5 * z + rnorm(50)
#' partialCorrelation(x, y, z)
partialCorrelation <- function(x, y = NULL, z = NULL, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (is.null(y) && is.null(z)) {
    m <- as.matrix(x)
    if (ncol(m) < 3L) stop("a single-argument `x` must have >= 3 columns (x, y, covariate[s]).", call. = FALSE)
    xx <- as.numeric(m[, 1L]); yy <- as.numeric(m[, 2L]); zz <- m[, 3:ncol(m), drop = FALSE]
  } else {
    xx <- as.numeric(x); yy <- as.numeric(y); zz <- as.matrix(z)
  }
  storage.mode(zz) <- "double"
  ok <- is.finite(xx) & is.finite(yy) & apply(is.finite(zz), 1L, all)
  xx <- xx[ok]; yy <- yy[ok]; zz <- zz[ok, , drop = FALSE]
  n <- length(xx); k <- ncol(zz)
  if (n <= k + 2L) stop("not enough complete observations for the covariates.", call. = FALSE)
  if (method == "spearman") {
    xx <- rank(xx); yy <- rank(yy); zz <- apply(zz, 2L, rank); dim(zz) <- c(n, k)
  }
  qrZ <- qr(cbind(1, zz))                       # residualise x and y on [1, covariates]
  rx <- qr.resid(qrZ, xx); ry <- qr.resid(qrZ, yy)
  r  <- stats::cor(rx, ry)
  df <- n - 2L - k
  t  <- r * sqrt(df / (1 - r^2))
  list(statistic = r, p_value = 2 * stats::pt(-abs(t), df), df = df, n = n, k = k)
}

#' Tukey HSD post-hoc test (all pairwise comparisons after one-way ANOVA)
#'
#' Tukey's Honest Significant Difference test: all pairwise mean comparisons
#' following a one-way ANOVA (\code{\link{oneWayAnova}}), with p-values and
#' confidence intervals from the studentized-range distribution and single-step
#' control of the family-wise error rate. Reproduces base R
#' \code{stats::TukeyHSD} and \code{scipy.stats.tukey_hsd} (machine precision).
#'
#' @param x Either a named list of `k >= 2` numeric group vectors, or the first
#'   group's numeric vector (with further groups passed through `...` -- use a
#'   named list to label the comparisons).
#' @param ... Additional numeric group vectors when `x` is a single vector.
#' @param conf_level Confidence level for the intervals (default 0.95).
#' @return A data frame with one row per pair: `comparison`, `diff` (mean
#'   difference), `lwr`/`upr` (CI) and `p_adj` (family-wise-adjusted p-value),
#'   carrying `k`, `df` and `n` as attributes.
#' @seealso \code{\link{oneWayAnova}}
#' @export
#' @examples
#' tukeyHSD(list(a = rnorm(20), b = rnorm(20, 1), c = rnorm(20, 2)))
tukeyHSD <- function(x, ..., conf_level = 0.95) {
  groups <- if (is.list(x) && !is.data.frame(x)) x else c(list(x), list(...))
  if (is.null(names(groups)) || any(!nzchar(names(groups))))
    names(groups) <- paste0("g", seq_along(groups))
  groups <- lapply(groups, function(g) as.numeric(g[is.finite(g)]))
  k <- length(groups)
  if (k < 2L) stop("Tukey HSD needs at least 2 groups.", call. = FALSE)
  ni <- vapply(groups, length, integer(1)); N <- sum(ni)
  mi <- vapply(groups, mean, numeric(1)); nm <- names(groups)
  gm <- sum(vapply(groups, sum, numeric(1))) / N
  ssw <- sum(vapply(groups, function(g) sum((g - mean(g))^2), numeric(1)))
  df <- N - k; mse <- ssw / df
  q <- stats::qtukey(conf_level, k, df)
  rows <- list()
  for (j in 2:k) for (i in 1:(j - 1L)) {                 # base R order: (i, j), i < j
    d  <- unname(mi[j] - mi[i])
    se <- sqrt((mse / 2) * (1 / ni[j] + 1 / ni[i]))
    w  <- q * se
    p  <- stats::ptukey(abs(d) / se, k, df, lower.tail = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      comparison = paste(nm[j], nm[i], sep = "-"),
      diff = d, lwr = d - w, upr = d + w, p_adj = p, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  attr(out, "k") <- k; attr(out, "df") <- df; attr(out, "n") <- N
  out
}

#' Dunn's test (post-hoc pairwise comparisons after Kruskal-Wallis)
#'
#' Dunn's test: all pairwise rank-based comparisons following a Kruskal-Wallis
#' test (\code{\link{kruskalTest}}), using the shared tie-corrected rank variance
#' and a chosen multiple-comparison adjustment. It is the non-parametric
#' counterpart of Tukey's HSD (\code{\link{tukeyHSD}}). The z-statistics and
#' adjusted p-values reproduce base R's Dunn formula, the `PMCMRplus` package and
#' `scikit_posthocs.posthoc_dunn` to machine precision.
#'
#' @param x Either a named list of `k >= 2` numeric group vectors, or the first
#'   group's numeric vector (with further groups passed through `...`).
#' @param ... Additional numeric group vectors when `x` is a single vector.
#' @param p_adjust Multiple-comparison adjustment: "bonferroni" (default),
#'   "none", "holm" or "BH".
#' @return A data frame with one row per pair: `comparison`, `z`,
#'   `p_value` (two-sided, unadjusted) and `p_adj`, carrying `k`, `n` as attributes.
#' @seealso \code{\link{kruskalTest}}, \code{\link{tukeyHSD}}
#' @export
#' @examples
#' dunnTest(list(a = rnorm(20), b = rnorm(20, 1), c = rnorm(20, 2)))
dunnTest <- function(x, ..., p_adjust = c("bonferroni", "none", "holm", "BH")) {
  p_adjust <- match.arg(p_adjust)
  groups <- if (is.list(x) && !is.data.frame(x)) x else c(list(x), list(...))
  if (is.null(names(groups)) || any(!nzchar(names(groups))))
    names(groups) <- paste0("g", seq_along(groups))
  groups <- lapply(groups, function(g) as.numeric(g[is.finite(g)]))
  k <- length(groups); nm <- names(groups)
  if (k < 2L) stop("Dunn's test needs at least 2 groups.", call. = FALSE)
  y <- unlist(groups, use.names = FALSE)
  grp <- rep(seq_len(k), vapply(groups, length, integer(1)))
  N <- length(y); r <- rank(y)
  tt <- table(y); sigma <- sqrt((N * (N + 1) / 12) - sum(tt^3 - tt) / (12 * (N - 1)))
  Rbar <- vapply(seq_len(k), function(i) mean(r[grp == i]), numeric(1))
  ni <- vapply(groups, length, integer(1))
  rows <- list()
  for (j in 2:k) for (i in 1:(j - 1L)) {                 # (i, j), i < j; report j-i
    z <- (Rbar[j] - Rbar[i]) / (sigma * sqrt(1 / ni[j] + 1 / ni[i]))
    rows[[length(rows) + 1L]] <- data.frame(
      comparison = paste(nm[j], nm[i], sep = "-"),
      z = z, p_value = 2 * stats::pnorm(-abs(z)), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  attr(out, "k") <- k; attr(out, "n") <- N
  out
}

#' Nemenyi post-hoc test (all pairwise comparisons after a Friedman test)
#'
#' The Nemenyi test: all pairwise comparisons of the treatments following a
#' Friedman repeated-measures test (\code{\link{friedmanTest}}), from the
#' within-block mean ranks and the studentized-range distribution with
#' single-step control of the family-wise error rate. It is the repeated-measures
#' counterpart of Tukey's HSD (\code{\link{tukeyHSD}}) and Dunn's test
#' (\code{\link{dunnTest}}). The statistics and p-values reproduce base R's
#' Nemenyi formula, the `PMCMRplus` package and
#' `scikit_posthocs.posthoc_nemenyi_friedman` to machine precision.
#'
#' @param x A blocks-by-treatments numeric matrix (rows = subjects/blocks,
#'   columns = conditions/treatments); rows with any missing value are dropped.
#' @return A data frame with one row per pair: `comparison`, `statistic` and
#'   `p_adj`, carrying `k` (treatments) and `n` (complete blocks) as attributes.
#' @seealso \code{\link{friedmanTest}}, \code{\link{tukeyHSD}}, \code{\link{dunnTest}}
#' @export
#' @examples
#' m <- matrix(rnorm(40) + rep(c(0, 0.5, 1, 1.5), each = 10), nrow = 10)
#' colnames(m) <- c("a", "b", "c", "d"); nemenyiTest(m)
nemenyiTest <- function(x) {
  M <- as.matrix(x)
  M <- M[stats::complete.cases(M), , drop = FALSE]
  storage.mode(M) <- "double"
  N <- nrow(M); k <- ncol(M)
  if (k < 2L) stop("Nemenyi test needs at least 2 treatments (columns).", call. = FALSE)
  if (N < 2L) stop("Nemenyi test needs at least 2 complete blocks (rows).", call. = FALSE)
  nm <- colnames(M); if (is.null(nm)) nm <- paste0("t", seq_len(k))
  Rbar <- colMeans(t(apply(M, 1L, rank)))
  se <- sqrt(k * (k + 1) / (6 * N))
  rows <- list()
  for (j in 2:k) for (i in 1:(j - 1L)) {                 # (i, j), i < j; report j vs i
    q <- abs(Rbar[j] - Rbar[i]) / se
    rows[[length(rows) + 1L]] <- data.frame(
      comparison = paste(nm[j], nm[i], sep = " vs "),
      statistic = unname(q),
      p_adj = stats::ptukey(q * sqrt(2), k, Inf, lower.tail = FALSE),
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  attr(out, "k") <- k; attr(out, "n") <- N
  out
}

#' Multiple linear regression (ordinary least squares)
#'
#' Multivariable OLS regression of a response on two or more predictors: the
#' coefficient table (estimate, standard error, t-statistic, p-value) plus the
#' overall model fit (R-squared, adjusted R-squared and the F-test). The
#' many-predictor generalisation of \code{\link{linearRegression}}. Computed via
#' the same QR decomposition as base R, so the coefficients, R-squared and F
#' reproduce `stats::lm` / `summary.lm` bit-for-bit and `statsmodels` OLS to
#' machine precision.
#'
#' @param x A matrix or data frame whose first column is the response and whose
#'   remaining columns are the predictors. Incomplete rows are dropped.
#' @return A data frame with one row per term (`term`, `estimate`, `std_error`,
#'   `statistic`, `p_value`), carrying `r_squared`, `adj_r_squared`,
#'   `f_statistic`, `df1`, `df2`, `f_p_value` and `n` as attributes.
#' @seealso \code{\link{linearRegression}}, \code{\link{partialCorrelation}}
#' @export
#' @examples
#' x <- matrix(rnorm(120), ncol = 3); colnames(x) <- c("y", "a", "b")
#' multipleRegression(x)
multipleRegression <- function(x) {
  m <- as.matrix(x)
  if (ncol(m) < 3L) stop("need a response and >= 2 predictors (>= 3 columns).", call. = FALSE)
  m <- m[stats::complete.cases(m), , drop = FALSE]
  storage.mode(m) <- "double"
  y  <- m[, 1L]; P <- m[, -1L, drop = FALSE]
  nm <- colnames(m); if (is.null(nm)) nm <- c("y", paste0("x", seq_len(ncol(P))))
  n  <- length(y); X <- cbind(`(Intercept)` = 1, P)
  p  <- ncol(X); df2 <- n - p
  if (df2 < 1L) stop("not enough observations for the predictors.", call. = FALSE)
  qrX <- qr(X)
  b   <- qr.coef(qrX, y)
  resid <- y - X %*% b; rss <- sum(resid^2)
  sigma2 <- rss / df2
  xtxinv <- chol2inv(qr.R(qrX))                 # same covariance basis as stats::lm
  se  <- sqrt(diag(sigma2 * xtxinv))
  tval <- as.numeric(b) / se
  tss <- sum((y - mean(y))^2)
  r2  <- 1 - rss / tss; adjr2 <- 1 - (1 - r2) * (n - 1) / df2
  df1 <- p - 1L
  fst <- ((tss - rss) / df1) / sigma2
  out <- data.frame(term = c("(Intercept)", colnames(P)),
                    estimate = as.numeric(b), std_error = as.numeric(se),
                    statistic = as.numeric(tval),
                    p_value = 2 * stats::pt(-abs(as.numeric(tval)), df2),
                    stringsAsFactors = FALSE)
  attr(out, "r_squared") <- r2; attr(out, "adj_r_squared") <- adjr2
  attr(out, "f_statistic") <- fst; attr(out, "df1") <- df1; attr(out, "df2") <- df2
  attr(out, "f_p_value") <- stats::pf(fst, df1, df2, lower.tail = FALSE); attr(out, "n") <- n
  out
}

#' Cronbach's alpha (internal-consistency reliability)
#'
#' Cronbach's alpha, the standard index of a scale's internal-consistency
#' reliability: how consistently a set of `k` items measures a single construct.
#' Computed from the item variances and the total-score variance. The raw alpha
#' and the mean inter-item correlation reproduce base R, the `psych` package
#' (`psych::alpha` raw_alpha) and `pingouin.cronbach_alpha` bit-for-bit.
#'
#' @param x A respondents-by-items numeric matrix or data frame (rows =
#'   respondents, columns = scale items). Rows with any missing value are dropped.
#' @return A list with `alpha` (raw Cronbach's alpha), `average_r` (mean
#'   inter-item correlation), `k` (items) and `n` (respondents).
#' @export
#' @examples
#' cronbachAlpha(matrix(round(runif(150, 1, 5)), ncol = 5))
cronbachAlpha <- function(x) {
  X <- as.matrix(x)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  storage.mode(X) <- "double"
  k <- ncol(X); n <- nrow(X)
  if (k < 2L) stop("need at least 2 items (columns).", call. = FALSE)
  if (n < 2L) stop("need at least 2 respondents (rows).", call. = FALSE)
  vi <- apply(X, 2L, stats::var)                 # per-item sample variance
  vt <- stats::var(rowSums(X))                   # total-score variance
  alpha <- (k / (k - 1)) * (1 - sum(vi) / vt)
  R <- stats::cor(X)
  list(alpha = alpha, average_r = mean(R[upper.tri(R)]), k = k, n = n)
}

#' Intraclass correlation coefficients (ICC)
#'
#' The six McGraw & Wong (1996) intraclass correlation coefficients from a
#' subjects-by-raters (or subjects-by-repeated-measurements) matrix: ICC1, ICC2,
#' ICC3 (single measurement) and ICC1k, ICC2k, ICC3k (average of `k`
#' measurements). The standard index of test-retest / inter-rater reliability,
#' computed from a two-way ANOVA. Reproduces the `psych` package
#' (`psych::ICC(lmer = FALSE)`) to machine precision and
#' `pingouin.intraclass_corr` bit-for-bit.
#'
#' @param x A subjects-by-raters numeric matrix (rows = targets/subjects, columns
#'   = raters or repeated measurements). Rows with any missing value are dropped.
#' @return A data frame with one row per ICC form (`form`, `ICC`), carrying `n`
#'   (subjects) and `k` (raters/measurements) as attributes.
#' @export
#' @examples
#' m <- matrix(rnorm(30), ncol = 3) + rnorm(10)
#' intraclassCorrelation(m)
intraclassCorrelation <- function(x) {
  M <- as.matrix(x)
  M <- M[stats::complete.cases(M), , drop = FALSE]
  storage.mode(M) <- "double"
  n <- nrow(M); k <- ncol(M)
  if (k < 2L) stop("need at least 2 raters/measurements (columns).", call. = FALSE)
  if (n < 2L) stop("need at least 2 subjects (rows).", call. = FALSE)
  gm <- mean(M)
  SST <- sum((M - gm)^2)
  SSR <- k * sum((rowMeans(M) - gm)^2)               # between-subjects
  SSC <- n * sum((colMeans(M) - gm)^2)               # between-raters
  SSE <- SST - SSR - SSC
  MSR <- SSR / (n - 1); MSC <- SSC / (k - 1)
  MSE <- SSE / ((n - 1) * (k - 1)); MSW <- (SST - SSR) / (n * (k - 1))
  icc <- c(
    ICC1  = (MSR - MSW) / (MSR + (k - 1) * MSW),
    ICC2  = (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n),
    ICC3  = (MSR - MSE) / (MSR + (k - 1) * MSE),
    ICC1k = (MSR - MSW) / MSR,
    ICC2k = (MSR - MSE) / (MSR + (MSC - MSE) / n),
    ICC3k = (MSR - MSE) / MSR)
  out <- data.frame(form = names(icc), ICC = as.numeric(icc), stringsAsFactors = FALSE)
  attr(out, "n") <- n; attr(out, "k") <- k
  out
}

#' Levene's test for homogeneity of variance
#'
#' Levene's test of whether `k >= 2` groups have equal variance, as a one-way
#' ANOVA on the absolute deviations from each group's centre. With
#' `center = "median"` (default) it is the robust Brown-Forsythe variant,
#' appropriate when the data may be non-normal; `center = "mean"` is the classic
#' Levene test. The F-statistic and p reproduce `scipy.stats.levene` and the
#' `car` package (`car::leveneTest`) bit-for-bit.
#'
#' @param x Either a list of `k >= 2` numeric group vectors, or the first
#'   group's numeric vector (with further groups passed through `...`).
#' @param ... Additional numeric group vectors when `x` is a single vector.
#' @param center "median" (default, Brown-Forsythe) or "mean" (classic Levene).
#' @return A list with `statistic` (F), `p_value`, `df1` (`k - 1`), `df2`
#'   (`N - k`), `k`, `n` and `center`.
#' @seealso \code{\link{oneWayAnova}}
#' @export
#' @examples
#' leveneTest(list(rnorm(20), rnorm(20, 0, 2), rnorm(20, 0, 0.5)))
leveneTest <- function(x, ..., center = c("median", "mean")) {
  center <- match.arg(center)
  groups <- if (is.list(x) && !is.data.frame(x)) x else c(list(x), list(...))
  groups <- lapply(groups, function(g) { g <- as.numeric(g); g[is.finite(g)] })
  k <- length(groups)
  if (k < 2L) stop("Levene's test needs at least 2 groups.", call. = FALSE)
  ctr <- if (center == "median") stats::median else mean
  z <- lapply(groups, function(g) abs(g - ctr(g)))           # |deviation from group centre|
  ni <- vapply(z, length, integer(1)); N <- sum(ni)
  mi <- vapply(z, mean, numeric(1)); gm <- sum(vapply(z, sum, numeric(1))) / N
  ssb <- sum(ni * (mi - gm)^2)
  ssw <- sum(vapply(z, function(zz) sum((zz - mean(zz))^2), numeric(1)))
  df1 <- k - 1L; df2 <- N - k
  Fv <- (ssb / df1) / (ssw / df2)
  list(statistic = Fv, p_value = stats::pf(Fv, df1, df2, lower.tail = FALSE),
       df1 = df1, df2 = df2, k = k, n = N, center = center)
}
