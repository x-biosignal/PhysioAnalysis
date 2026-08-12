# Compute cross-spectral density

Calculates the cross-spectral density between channel pairs.

## Usage

``` r
crossSpectrum(
  x,
  channels = NULL,
  nperseg = 256L,
  noverlap = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- channels:

  Integer vector of channel indices. If NULL, uses all.

- nperseg:

  Number of samples per segment. Default is 256.

- noverlap:

  Number of overlapping samples.

- assay_name:

  Input assay name.

## Value

A list with the following components:

- csd:

  3D complex array of cross-spectral density values (frequency x channel
  x channel)

- frequencies:

  Numeric vector of frequencies in Hz

- channel_names:

  Character vector of channel names

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioAnalysis/reference/coherence.md)
for magnitude-squared coherence,
[`connectivityMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/connectivityMatrix.md)
for a unified connectivity interface,
[`spectrogram()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spectrogram.md)
for single-channel spectral analysis.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(2000), nrow = 500, ncol = 4)),
  samplingRate = 256
)

# Compute cross-spectral density
csd <- crossSpectrum(pe)
```
