# Spectral Analysis with PhysioAnalysis

## Introduction

PhysioAnalysis provides a comprehensive set of spectral analysis tools
for physiological signal data stored in `PhysioExperiment` objects. This
vignette covers three core spectral analysis workflows:

1.  **Fast Fourier Transform (FFT)** – computing magnitude spectra
2.  **Power Spectral Density (PSD)** – estimating and visualizing signal
    power across frequencies
3.  **Time-frequency analysis** – spectrograms and wavelet transforms
    for non-stationary signals

## Setup

``` r

library(PhysioCore)
library(PhysioAnalysis)
```

## Creating Example Data

We begin by constructing a synthetic `PhysioExperiment` object with a
known frequency composition so that we can verify our spectral analyses.

``` r

set.seed(42)
sr <- 256  # sampling rate in Hz
duration <- 4  # seconds
n_time <- sr * duration
t <- seq(0, duration - 1/sr, by = 1/sr)

# Simulate a 4-channel signal with known frequency components:
#   Channel 1: 10 Hz sine (alpha)
#   Channel 2: 10 Hz + 25 Hz sine (alpha + beta)
#   Channel 3: White noise
#   Channel 4: 10 Hz sine + noise
signal <- matrix(0, nrow = n_time, ncol = 4)
signal[, 1] <- sin(2 * pi * 10 * t)
signal[, 2] <- sin(2 * pi * 10 * t) + 0.5 * sin(2 * pi * 25 * t)
signal[, 3] <- rnorm(n_time)
signal[, 4] <- sin(2 * pi * 10 * t) + 0.3 * rnorm(n_time)

pe <- PhysioExperiment(
  assays = list(raw = signal),
  colData = S4Vectors::DataFrame(
    label = c("Ch1_alpha", "Ch2_alpha_beta", "Ch3_noise", "Ch4_alpha_noise")
  ),
  samplingRate = sr
)
```

## FFT Analysis

The
[`fftSignals()`](https://x-biosignal.github.io/PhysioAnalysis/reference/fftSignals.md)
function computes the discrete Fourier transform along the time axis of
the default assay and stores the magnitude spectrum in a new assay named
`"fft"`.

``` r

pe_fft <- fftSignals(pe)

# The FFT result is stored as a new assay
assayNames(pe_fft)

# Access the magnitude spectrum
fft_data <- SummarizedExperiment::assay(pe_fft, "fft")
dim(fft_data)  # time x channels (frequency bins x channels)
```

### Interpreting the FFT Output

The magnitude spectrum returned by
[`fftSignals()`](https://x-biosignal.github.io/PhysioAnalysis/reference/fftSignals.md)
has the same dimensions as the input data. The first half of the rows
corresponds to frequencies from 0 Hz to the Nyquist frequency (sampling
rate / 2).

``` r

n <- nrow(fft_data)
freqs <- seq(0, sr / 2, length.out = floor(n / 2) + 1)

# Plot magnitude spectrum for channel 1 (should peak at 10 Hz)
plot(freqs, fft_data[1:(floor(n/2) + 1), 1],
     type = "l", xlab = "Frequency (Hz)", ylab = "Magnitude",
     main = "FFT Magnitude - Channel 1 (10 Hz)")
```

## Power Spectral Density

The
[`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md)
function provides a convenient way to compute and visualize the power
spectral density for one or more channels.

``` r

# Plot PSD for all channels with log scale
plotPSD(pe, channels = 1:4, log_scale = TRUE)

# Plot PSD for a specific frequency range
plotPSD(pe, channels = c(1, 2), freq_range = c(1, 50), log_scale = TRUE)

# Linear scale PSD
plotPSD(pe, channels = 1:2, log_scale = FALSE, freq_range = c(1, 50))
```

## Band Power Extraction

The
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md)
function extracts power within predefined or custom frequency bands.
This is useful for quantifying activity in standard EEG bands (delta,
theta, alpha, beta, gamma).

``` r

# Extract power in standard EEG bands
bp <- bandPower(pe, bands = list(
  delta = c(1, 4),
  theta = c(4, 8),
  alpha = c(8, 13),
  beta  = c(13, 30),
  gamma = c(30, 100)
))

# The result contains power values for each band and channel
bp
```

## Spectrogram (Short-Time Fourier Transform)

The
[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md)
function computes a time-frequency representation of the signal using
the Short-Time Fourier Transform (STFT). This is useful for analyzing
how the frequency content of a signal changes over time.

``` r

# Compute spectrogram for channel 1
spec <- spectrogram(pe, window_size = 128, overlap = 0.75, channel = 1)

# The result contains:
#   - power: frequency x time power matrix
#   - frequencies: frequency vector (Hz)
#   - times: time vector (seconds)
str(spec)

# Visualize the spectrogram
plotSpectrogram(spec, freq_range = c(1, 50))
```

### Window Parameters

The choice of window size and overlap affects the trade-off between time
and frequency resolution:

- **Larger windows** provide better frequency resolution but poorer time
  resolution.
- **Higher overlap** provides smoother temporal estimates at the cost of
  increased computation.
- Common window functions include Hanning (default), Hamming, and
  Blackman.

``` r

# Higher frequency resolution (larger window)
spec_high_freq <- spectrogram(pe, window_size = 512, overlap = 0.75, channel = 1)
plotSpectrogram(spec_high_freq, freq_range = c(1, 50))

# Higher time resolution (smaller window)
spec_high_time <- spectrogram(pe, window_size = 64, overlap = 0.75, channel = 1)
plotSpectrogram(spec_high_time, freq_range = c(1, 50))
```

## Wavelet Transform

The
[`waveletTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/waveletTransform.md)
function uses Morlet wavelets to compute a time-frequency
representation. Compared to STFT, wavelet analysis provides better
frequency resolution at low frequencies and better time resolution at
high frequencies.

``` r

# Compute wavelet transform for channel 1
wt <- waveletTransform(pe, frequencies = seq(1, 50, by = 0.5), channel = 1)

# The result contains the wavelet power at each frequency and time point
str(wt)
```

## Hilbert Transform and Instantaneous Measures

The Hilbert transform computes the analytic signal, which can be used to
extract instantaneous amplitude (envelope) and instantaneous phase.

``` r

# Compute the analytic signal
pe_hilbert <- hilbertTransform(pe)

# Extract instantaneous amplitude (envelope)
pe_amp <- instantaneousAmplitude(pe_hilbert)

# Extract instantaneous phase
pe_phase <- instantaneousPhase(pe_hilbert)
```

## Summary

| Function | Purpose | Output |
|----|----|----|
| [`fftSignals()`](https://x-biosignal.github.io/PhysioAnalysis/reference/fftSignals.md) | Magnitude spectrum via FFT | New assay |
| [`plotPSD()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotPSD.md) | Power spectral density visualization | ggplot |
| [`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md) | Power in frequency bands | List |
| [`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md) | STFT time-frequency analysis | List |
| [`waveletTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/waveletTransform.md) | Wavelet time-frequency analysis | List |
| [`hilbertTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/hilbertTransform.md) | Analytic signal via Hilbert transform | New assay |
| [`plotSpectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotSpectrogram.md) | Spectrogram visualization | ggplot |

## References

- Oppenheim, A.V. & Willsky, A.S. (1997). *Signals and Systems*. 2nd ed.
  Prentice Hall.
- Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
  Springer-Verlag New York.

## Session Info

``` r

sessionInfo()
```
