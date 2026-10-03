#' Connectivity Analysis for PhysioExperiment
#'
#' Functions for computing functional connectivity between channels including
#' coherence, phase synchrony measures, and correlation-based metrics.

#' Compute coherence between channels
#'
#' Calculates the magnitude-squared coherence between pairs of channels,
#' which measures the linear correlation between signals as a function of
#' frequency.
#'
#' @param x A PhysioExperiment object.
#' @param channels Integer vector of channel indices to analyze. If NULL, uses all.
#' @param freq_range Numeric vector of length 2 specifying frequency range (Hz).
#' @param nperseg Number of samples per segment for Welch's method. Default is 256.
#' @param noverlap Number of overlapping samples. Default is nperseg/2.
#' @param assay_name Input assay name. If NULL, uses default assay.
#' @return A list with components:
#'   \item{coherence}{3D array (freq x channel x channel) of coherence values}
#'   \item{frequencies}{Frequency vector}
#'   \item{channel_names}{Channel names}
#' @details
#' Coherence is computed using Welch's averaged periodogram method.
#' Values range from 0 (no linear relationship) to 1 (perfect linear relationship).
#' @references Nolte, G., et al. (2004). "Identifying true brain interaction
#'   from EEG data using the imaginary part of coherency." Clinical
#'   Neurophysiology, 115(10), 2292-2307. doi:10.1016/j.clinph.2004.04.029
#' @seealso [crossSpectrum()] for the underlying cross-spectral density,
#'   [plv()] for phase-based connectivity, [connectivityMatrix()] for a
#'   unified connectivity interface.
#' @export
#' @examples
#' # Create example with 4 channels
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
#'   colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
#'   samplingRate = 256
#' )
#'
#' # Compute coherence
#' coh <- coherence(pe, freq_range = c(1, 50))
#' dim(coh$coherence)  # frequencies x channels x channels
coherence <- function(x, channels = NULL, freq_range = NULL,
                      nperseg = 256L, noverlap = NULL,
                      assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  # Handle different dimensions
  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) >= 3) {
    signal_data <- data[, , 1]
  }

  # Get channels
  if (is.null(channels)) {
    channels <- seq_len(ncol(signal_data))
  }
  n_channels <- length(channels)
  signal_data <- signal_data[, channels, drop = FALSE]

  ch_names <- channelNames(x)
  if (length(ch_names) >= max(channels)) {
    ch_names <- ch_names[channels]
  } else {
    ch_names <- paste0("Ch", channels)
  }

  if (is.null(noverlap)) {
    noverlap <- as.integer(nperseg / 2)
  }

  # Compute cross-spectral density for all pairs
  n_freqs <- floor(nperseg / 2) + 1
  freqs <- seq(0, sr / 2, length.out = n_freqs)

  # Apply frequency range filter
  if (!is.null(freq_range)) {
    freq_idx <- which(freqs >= freq_range[1] & freqs <= freq_range[2])
    freqs <- freqs[freq_idx]
  } else {
    freq_idx <- seq_len(n_freqs)
  }

  coh_array <- array(NA_real_, dim = c(length(freq_idx), n_channels, n_channels))

  # Compute PSD for each channel and cross-spectrum for each pair
  psd_list <- list()
  for (i in seq_len(n_channels)) {
    psd_list[[i]] <- .welchPSD2(signal_data[, i], nperseg, noverlap, sr)
  }

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      if (i == j) {
        coh_array[, i, j] <- 1  # Perfect coherence with self
      } else {
        csd <- .crossSpectralDensity(signal_data[, i], signal_data[, j],
                                     nperseg, noverlap, sr)
        pxx <- psd_list[[i]]$psd
        pyy <- psd_list[[j]]$psd

        # Magnitude-squared coherence: |Pxy|^2 / (Pxx * Pyy)
        coh <- (Mod(csd$csd)^2) / (pxx * pyy)
        coh_array[, i, j] <- coh[freq_idx]
      }
    }
  }

  list(
    coherence = coh_array,
    frequencies = freqs,
    channel_names = ch_names
  )
}

#' Compute cross-spectral density
#'
#' Calculates the cross-spectral density between channel pairs.
#'
#' @param x A PhysioExperiment object.
#' @param channels Integer vector of channel indices. If NULL, uses all.
#' @param nperseg Number of samples per segment. Default is 256.
#' @param noverlap Number of overlapping samples.
#' @param assay_name Input assay name.
#' @return A list with the following components:
#'   \item{csd}{3D complex array of cross-spectral density values
#'     (frequency x channel x channel)}
#'   \item{frequencies}{Numeric vector of frequencies in Hz}
#'   \item{channel_names}{Character vector of channel names}
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [coherence()] for magnitude-squared coherence,
#'   [connectivityMatrix()] for a unified connectivity interface,
#'   [spectrogram()] for single-channel spectral analysis.
#' @export
#' @examples
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(2000), nrow = 500, ncol = 4)),
#'   samplingRate = 256
#' )
#'
#' # Compute cross-spectral density
#' csd <- crossSpectrum(pe)
crossSpectrum <- function(x, channels = NULL, nperseg = 256L,
                          noverlap = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) >= 3) {
    signal_data <- data[, , 1]
  }

  if (is.null(channels)) {
    channels <- seq_len(ncol(signal_data))
  }
  n_channels <- length(channels)
  signal_data <- signal_data[, channels, drop = FALSE]

  ch_names <- channelNames(x)
  if (length(ch_names) >= max(channels)) {
    ch_names <- ch_names[channels]
  } else {
    ch_names <- paste0("Ch", channels)
  }

  if (is.null(noverlap)) {
    noverlap <- as.integer(nperseg / 2)
  }

  n_freqs <- floor(nperseg / 2) + 1
  freqs <- seq(0, sr / 2, length.out = n_freqs)

  csd_array <- array(complex(real = 0, imaginary = 0),
                     dim = c(n_freqs, n_channels, n_channels))

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      csd <- .crossSpectralDensity(signal_data[, i], signal_data[, j],
                                   nperseg, noverlap, sr)
      csd_array[, i, j] <- csd$csd
    }
  }

  list(
    csd = csd_array,
    frequencies = freqs,
    channel_names = ch_names
  )
}

#' Compute Phase Locking Value (PLV)
#'
#' Calculates the Phase Locking Value between channel pairs, measuring the
#' consistency of phase difference across time.
#'
#' @param x A PhysioExperiment object.
#' @param freq_band Numeric vector of length 2 specifying frequency band (Hz).
#' @param channels Integer vector of channel indices.
#' @param assay_name Input assay name.
#' @return A numeric matrix of PLV values (channel x channel) with values
#'   ranging from 0 (random phase relationship) to 1 (perfect phase locking).
#'   Row and column names are set to channel names when available.
#' @details
#' PLV measures the consistency of the phase difference between two signals.
#' A value of 1 indicates perfect phase locking, while 0 indicates random
#' phase relationship.
#'
#' The signals are first bandpass filtered to the specified frequency band,
#' then the analytic signal is computed using the Hilbert transform.
#' @references Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in
#'   brain signals." Human Brain Mapping, 8(4), 194-208.
#'   doi:10.1002/(SICI)1097-0193(1999)8:4<194::AID-HBM4>3.0.CO;2-C
#' @seealso [pli()] for phase lag index, [wPLI()] for weighted phase lag index,
#'   [coherence()] for frequency-domain coherence, [connectivityMatrix()] for
#'   a unified connectivity interface.
#' @export
#' @examples
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
#'   colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
#'   samplingRate = 256
#' )
#'
#' # Compute PLV in alpha band (8-12 Hz)
#' plv_matrix <- plv(pe, freq_band = c(8, 12))
plv <- function(x, freq_band, channels = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (length(freq_band) != 2) {
    stop("freq_band must be a vector of length 2", call. = FALSE)
  }

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) >= 3) {
    signal_data <- data[, , 1]
  }

  if (is.null(channels)) {
    channels <- seq_len(ncol(signal_data))
  }
  n_channels <- length(channels)
  signal_data <- signal_data[, channels, drop = FALSE]

  n_samples <- nrow(signal_data)

  # Bandpass filter and extract phase for each channel
  phases <- matrix(NA_real_, nrow = n_samples, ncol = n_channels)

  for (i in seq_len(n_channels)) {
    # Bandpass filter
    filtered <- .bandpassFilter(signal_data[, i], freq_band[1], freq_band[2], sr)
    # Hilbert transform for phase
    analytic <- .hilbertTransformSignal(filtered)
    phases[, i] <- Arg(analytic)
  }

  # Compute PLV for each pair
  plv_matrix <- matrix(1, nrow = n_channels, ncol = n_channels)

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      if (i != j) {
        phase_diff <- phases[, i] - phases[, j]
        # PLV = |mean(exp(i * phase_diff))|
        plv_matrix[i, j] <- abs(mean(exp(complex(imaginary = phase_diff))))
      }
    }
  }

  ch_names <- channelNames(x)
  if (length(ch_names) >= max(channels)) {
    rownames(plv_matrix) <- ch_names[channels]
    colnames(plv_matrix) <- ch_names[channels]
  }

  plv_matrix
}

#' Compute Phase Lag Index (PLI)
#'
#' Calculates the Phase Lag Index between channel pairs, a measure of
#' asymmetry in the phase difference distribution.
#'
#' @param x A PhysioExperiment object.
#' @param freq_band Numeric vector of length 2 specifying frequency band (Hz).
#' @param channels Integer vector of channel indices.
#' @param assay_name Input assay name.
#' @return A numeric matrix of PLI values (channel x channel) with values
#'   ranging from 0 (no asymmetry in phase distribution) to 1 (maximal
#'   asymmetry). Row and column names are set to channel names when available.
#' @details
#' PLI measures the asymmetry of the distribution of phase differences.
#' Unlike PLV, PLI is insensitive to volume conduction effects that lead
#' to zero-lag synchronization.
#'
#' PLI = |mean(sign(Im(S_xy)))|, where S_xy is the cross-spectrum.
#' @references Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in
#'   brain signals." Human Brain Mapping, 8(4), 194-208.
#'   doi:10.1002/(SICI)1097-0193(1999)8:4<194::AID-HBM4>3.0.CO;2-C
#' @seealso [plv()] for phase locking value, [wPLI()] for weighted phase lag
#'   index, [coherence()] for frequency-domain coherence,
#'   [connectivityMatrix()] for a unified connectivity interface.
#' @export
#' @examples
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
#'   colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
#'   samplingRate = 256
#' )
#'
#' # Compute PLI in alpha band (8-12 Hz)
#' pli_matrix <- pli(pe, freq_band = c(8, 12))
pli <- function(x, freq_band, channels = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (length(freq_band) != 2) {
    stop("freq_band must be a vector of length 2", call. = FALSE)
  }

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) >= 3) {
    signal_data <- data[, , 1]
  }

  if (is.null(channels)) {
    channels <- seq_len(ncol(signal_data))
  }
  n_channels <- length(channels)
  signal_data <- signal_data[, channels, drop = FALSE]

  n_samples <- nrow(signal_data)

  # Bandpass filter and extract phase for each channel
  phases <- matrix(NA_real_, nrow = n_samples, ncol = n_channels)

  for (i in seq_len(n_channels)) {
    filtered <- .bandpassFilter(signal_data[, i], freq_band[1], freq_band[2], sr)
    analytic <- .hilbertTransformSignal(filtered)
    phases[, i] <- Arg(analytic)
  }

  # Compute PLI for each pair
  pli_matrix <- matrix(0, nrow = n_channels, ncol = n_channels)

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      if (i != j) {
        phase_diff <- phases[, i] - phases[, j]
        # PLI = |mean(sign(sin(phase_diff)))|
        pli_matrix[i, j] <- abs(mean(sign(sin(phase_diff))))
      }
    }
  }

  ch_names <- channelNames(x)
  if (length(ch_names) >= max(channels)) {
    rownames(pli_matrix) <- ch_names[channels]
    colnames(pli_matrix) <- ch_names[channels]
  }

  pli_matrix
}

#' Compute weighted Phase Lag Index (wPLI)
#'
#' Calculates the weighted Phase Lag Index, which is less sensitive to
#' noise than standard PLI.
#'
#' @param x A PhysioExperiment object.
#' @param freq_band Numeric vector of length 2 specifying frequency band (Hz).
#' @param channels Integer vector of channel indices.
#' @param assay_name Input assay name.
#' @return A numeric matrix of wPLI values (channel x channel) with values
#'   ranging from 0 to 1. Row and column names are set to channel names
#'   when available.
#' @details
#' wPLI weights the contribution of each phase difference by the magnitude
#' of the imaginary component, reducing the influence of noise sources.
#' @references Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in
#'   brain signals." Human Brain Mapping, 8(4), 194-208.
#'   doi:10.1002/(SICI)1097-0193(1999)8:4<194::AID-HBM4>3.0.CO;2-C
#' @seealso [plv()] for phase locking value, [pli()] for phase lag index,
#'   [coherence()] for frequency-domain coherence, [connectivityMatrix()] for
#'   a unified connectivity interface.
#' @export
#' @examples
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
#'   samplingRate = 256
#' )
#'
#' # Compute wPLI in theta band (4-8 Hz)
#' wpli_matrix <- wPLI(pe, freq_band = c(4, 8))
wPLI <- function(x, freq_band, channels = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (length(freq_band) != 2) {
    stop("freq_band must be a vector of length 2", call. = FALSE)
  }

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate required", call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) >= 3) {
    signal_data <- data[, , 1]
  }

  if (is.null(channels)) {
    channels <- seq_len(ncol(signal_data))
  }
  n_channels <- length(channels)
  signal_data <- signal_data[, channels, drop = FALSE]

  n_samples <- nrow(signal_data)

  # Bandpass filter and compute cross-spectrum phases
  filtered <- matrix(NA_real_, nrow = n_samples, ncol = n_channels)

  for (i in seq_len(n_channels)) {
    filtered[, i] <- .bandpassFilter(signal_data[, i], freq_band[1], freq_band[2], sr)
  }

  # Compute wPLI for each pair using cross-spectrum
  wpli_matrix <- matrix(0, nrow = n_channels, ncol = n_channels)

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      if (i != j) {
        analytic_i <- .hilbertTransformSignal(filtered[, i])
        analytic_j <- .hilbertTransformSignal(filtered[, j])

        # Cross-spectrum at each time point
        cross <- analytic_i * Conj(analytic_j)
        imag_cross <- Im(cross)

        # wPLI = |sum(|Im(cross)| * sign(Im(cross)))| / sum(|Im(cross)|)
        wpli_matrix[i, j] <- abs(sum(abs(imag_cross) * sign(imag_cross))) /
          sum(abs(imag_cross))
      }
    }
  }

  ch_names <- channelNames(x)
  if (length(ch_names) >= max(channels)) {
    rownames(wpli_matrix) <- ch_names[channels]
    colnames(wpli_matrix) <- ch_names[channels]
  }

  wpli_matrix
}

#' Compute correlation matrix between channels
#'
#' Calculates the Pearson correlation coefficient between all channel pairs.
#'
#' @param x A PhysioExperiment object.
#' @param channels Integer vector of channel indices.
#' @param method Correlation method: "pearson", "spearman", or "kendall".
#' @param assay_name Input assay name.
#' @return A numeric correlation matrix (channel x channel) with values
#'   ranging from -1 to 1. Row and column names are set to channel names
#'   when available.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [coherence()] for frequency-domain coherence, [plv()] for
#'   phase-based connectivity, [connectivityMatrix()] for a unified
#'   connectivity interface.
#' @export
#' @examples
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
#'   colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
#'   samplingRate = 256
#' )
#'
#' # Compute Pearson correlation
#' cor_matrix <- correlationMatrix(pe)
correlationMatrix <- function(x, channels = NULL,
                               method = c("pearson", "spearman", "kendall"),
                               assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) == 3) {
    signal_data <- data[, , 1]
  } else if (length(dims) == 4) {
    signal_data <- data[, , 1, 1]
  } else {
    stop("Unsupported array dimension", call. = FALSE)
  }

  if (is.null(channels)) {
    channels <- seq_len(ncol(signal_data))
  }
  signal_data <- signal_data[, channels, drop = FALSE]

  cor_mat <- stats::cor(signal_data, method = method, use = "pairwise.complete.obs")

  ch_names <- channelNames(x)
  if (length(ch_names) >= max(channels)) {
    rownames(cor_mat) <- ch_names[channels]
    colnames(cor_mat) <- ch_names[channels]
  }

  cor_mat
}

#' Compute connectivity matrix for a frequency band
#'
#' High-level function to compute connectivity using various metrics.
#'
#' @param x A PhysioExperiment object.
#' @param method Connectivity method: "coherence", "plv", "pli", "wpli", "correlation".
#' @param freq_band Frequency band for phase-based methods.
#' @param channels Integer vector of channel indices.
#' @param assay_name Input assay name.
#' @return A numeric matrix (channel x channel) of connectivity values.
#'   For `"coherence"`, returns the mean coherence across frequencies.
#'   For phase-based methods (`"plv"`, `"pli"`, `"wpli"`), returns
#'   values between 0 and 1. For `"correlation"`, returns values between
#'   -1 and 1.
#' @references Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in
#'   brain signals." Human Brain Mapping, 8(4), 194-208.
#'   doi:10.1002/(SICI)1097-0193(1999)8:4<194::AID-HBM4>3.0.CO;2-C
#' @seealso [coherence()], [plv()], [pli()], [wPLI()],
#'   [correlationMatrix()] for individual connectivity methods.
#' @export
#' @examples
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
#'   samplingRate = 256
#' )
#'
#' # Compute PLV connectivity in alpha band
#' conn <- connectivityMatrix(pe, method = "plv", freq_band = c(8, 12))
connectivityMatrix <- function(x, method = c("coherence", "plv", "pli", "wpli", "correlation"),
                               freq_band = NULL, channels = NULL, assay_name = NULL) {
  method <- match.arg(method)

  if (method %in% c("plv", "pli", "wpli") && is.null(freq_band)) {
    stop("freq_band required for phase-based connectivity methods", call. = FALSE)
  }

  result <- switch(method,
    "coherence" = {
      coh <- coherence(x, channels = channels, freq_range = freq_band,
                       assay_name = assay_name)
      # Return mean coherence across frequencies
      apply(coh$coherence, c(2, 3), mean, na.rm = TRUE)
    },
    "plv" = plv(x, freq_band = freq_band, channels = channels, assay_name = assay_name),
    "pli" = pli(x, freq_band = freq_band, channels = channels, assay_name = assay_name),
    "wpli" = wPLI(x, freq_band = freq_band, channels = channels, assay_name = assay_name),
    "correlation" = correlationMatrix(x, channels = channels, assay_name = assay_name)
  )

  result
}

#' Internal: Welch PSD estimation
#' @noRd
.welchPSD2 <- function(x, nperseg, noverlap, sr) {
  n <- length(x)
  step <- nperseg - noverlap
  n_segments <- floor((n - noverlap) / step)

  if (n_segments < 1) {
    n_segments <- 1
    nperseg <- n
    step <- n
  }

  n_freqs <- floor(nperseg / 2) + 1
  psd <- numeric(n_freqs)

  # Hanning window
  window <- 0.5 * (1 - cos(2 * pi * seq(0, nperseg - 1) / (nperseg - 1)))

  for (seg in seq_len(n_segments)) {
    start <- (seg - 1) * step + 1
    end <- start + nperseg - 1

    if (end > n) break

    segment <- x[start:end] * window
    fft_result <- stats::fft(segment)
    psd <- psd + (Mod(fft_result[1:n_freqs])^2)
  }

  psd <- psd / n_segments
  psd <- psd / (sr * sum(window^2))

  list(
    psd = psd,
    frequencies = seq(0, sr / 2, length.out = n_freqs)
  )
}

#' Internal: Cross-spectral density estimation
#' @noRd
.crossSpectralDensity <- function(x, y, nperseg, noverlap, sr) {
  n <- length(x)
  step <- nperseg - noverlap
  n_segments <- floor((n - noverlap) / step)

  if (n_segments < 1) {
    n_segments <- 1
    nperseg <- n
    step <- n
  }

  n_freqs <- floor(nperseg / 2) + 1
  csd <- complex(real = rep(0, n_freqs), imaginary = rep(0, n_freqs))

  # Hanning window
  window <- 0.5 * (1 - cos(2 * pi * seq(0, nperseg - 1) / (nperseg - 1)))

  for (seg in seq_len(n_segments)) {
    start <- (seg - 1) * step + 1
    end <- start + nperseg - 1

    if (end > n) break

    seg_x <- x[start:end] * window
    seg_y <- y[start:end] * window

    fft_x <- stats::fft(seg_x)
    fft_y <- stats::fft(seg_y)

    # Cross-spectrum: X * conj(Y)
    csd <- csd + fft_x[1:n_freqs] * Conj(fft_y[1:n_freqs])
  }

  csd <- csd / n_segments
  csd <- csd / (sr * sum(window^2))

  list(
    csd = csd,
    frequencies = seq(0, sr / 2, length.out = n_freqs)
  )
}

#' Internal: Simple bandpass filter
#' @noRd
.bandpassFilter <- function(x, low, high, sr) {
  n <- length(x)
  nyq <- sr / 2

  # Design filter coefficients (simple Butterworth)
  order <- 4
  Wn <- c(low / nyq, high / nyq)

  # Clamp to valid range
  Wn <- pmax(pmin(Wn, 0.99), 0.01)

  tryCatch({
    bf <- signal::butter(order, Wn, type = "pass")
    signal::filtfilt(bf, x)
  }, error = function(e) {
    # Fallback: simple FFT-based filter
    fft_x <- stats::fft(x)
    freqs <- seq(0, sr, length.out = n)

    # Create filter mask
    mask <- (freqs >= low & freqs <= high) | (freqs >= (sr - high) & freqs <= (sr - low))
    fft_x[!mask] <- 0

    Re(stats::fft(fft_x, inverse = TRUE)) / n
  })
}

#' Internal: Hilbert transform for analytic signal
#' @noRd
.hilbertTransformSignal <- function(x) {
  n <- length(x)
  fft_x <- stats::fft(x)

  # Create Hilbert transform multiplier
  h <- rep(0, n)
  if (n > 0) {
    h[1] <- 1  # DC component
    if (n %% 2 == 0) {
      h[2:(n / 2)] <- 2
      h[n / 2 + 1] <- 1  # Nyquist
    } else {
      h[2:((n + 1) / 2)] <- 2
    }
  }

  stats::fft(fft_x * h, inverse = TRUE) / n
}
