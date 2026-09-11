#' Fast Fourier transform helper
#'
#' Computes the discrete Fourier transform along the time axis of the default
#' assay and stores the magnitude spectrum in a new assay named `"fft"`.
#'
#' @param x A `PhysioExperiment` object.
#' @return A `PhysioExperiment` object with an additional `"fft"` assay
#'   containing the magnitude spectrum (same dimensions as the input assay).
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [spectrogram()] for time-frequency analysis, [bandPower()] for
#'   frequency band power extraction, [plotPSD()] for power spectral density
#'   visualization.
#' @export
fftSignals <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))
  assay_name <- defaultAssay(x)
  if (is.na(assay_name)) {
    stop("No assays available for FFT", call. = FALSE)
  }
  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)
  fft_data <- data
  # Batch the per-channel FFT through a single mvfft call (column-wise FFT of
  # the whole matrix in one C call), avoiding the apply()/array-permute
  # overhead of a per-column loop. Results are identical to
  # apply(data, 2, function(v) Mod(stats::fft(v))).
  if (length(dims) == 2L) {
    fft_data[] <- Mod(stats::mvfft(data))
  } else if (length(dims) == 1L) {
    fft_data[] <- Mod(stats::fft(as.numeric(data)))
  } else {
    # 3D (time x channels x samples): mvfft each sample slice. Keep the slice
    # a matrix even for a single channel (mvfft requires a matrix).
    for (s in seq_len(dims[3])) {
      slice <- matrix(data[, , s], nrow = dims[1], ncol = dims[2])
      fft_data[, , s] <- Mod(stats::mvfft(slice))
    }
  }
  assays <- SummarizedExperiment::assays(x)
  assays[["fft"]] <- fft_data
  SummarizedExperiment::assays(x) <- assays
  .recordProv(x, input_assay = defaultAssay(x), output_assay = "fft",
              .package = "PhysioAnalysis")
}
