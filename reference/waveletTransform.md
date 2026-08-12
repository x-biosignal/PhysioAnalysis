# Wavelet transform

Computes the continuous wavelet transform using Morlet wavelets.

## Usage

``` r
waveletTransform(
  x,
  frequencies = seq(1, 40, by = 1),
  n_cycles = 7,
  channel = 1L,
  sample = 1L,
  normalization = c("L2", "L1")
)
```

## Arguments

- x:

  A PhysioExperiment object.

- frequencies:

  Numeric vector of frequencies to analyze.

- n_cycles:

  Number of wavelet cycles (can be scalar or vector).

- channel:

  Channel index to analyze.

- sample:

  Sample index (for 3D data).

- normalization:

  Wavelet normalization method: `"L2"` (default, divides by square root
  of sum of squared absolute values) or `"L1"` (divides by sum of
  absolute values). L2 normalization preserves energy across frequencies
  and is preferred for power comparisons.

## Value

A list with the following components:

- power:

  Power matrix (frequency x time)

- phase:

  Phase matrix in radians (frequency x time)

- frequencies:

  Numeric vector of analyzed frequencies in Hz

- times:

  Numeric vector of time points in seconds

- sampling_rate:

  The sampling rate used

- n_cycles:

  Number of wavelet cycles per frequency

## References

Torrence, C. & Compo, G.P. (1998). "A practical guide to wavelet
analysis." Bulletin of the American Meteorological Society, 79(1),
61-78.

## See also

[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md)
for STFT-based time-frequency analysis,
[`bandPower()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bandPower.md)
for band power extraction,
[`hilbertTransform()`](https://x-biosignal.github.io/PhysioAnalysis/reference/hilbertTransform.md)
for analytic signal computation.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 4), nrow = 500)),
  samplingRate = 100
)

# Compute wavelet transform (1-30 Hz)
wt <- waveletTransform(pe, frequencies = seq(1, 30), channel = 1)

# Access power and phase
dim(wt$power)  # frequency x time
#> [1]  30 500
```
