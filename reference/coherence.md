# Connectivity Analysis for PhysioExperiment

Functions for computing functional connectivity between channels
including coherence, phase synchrony measures, and correlation-based
metrics. Compute coherence between channels

## Usage

``` r
coherence(
  x,
  channels = NULL,
  freq_range = NULL,
  nperseg = 256L,
  noverlap = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- channels:

  Integer vector of channel indices to analyze. If NULL, uses all.

- freq_range:

  Numeric vector of length 2 specifying frequency range (Hz).

- nperseg:

  Number of samples per segment for Welch's method. Default is 256.

- noverlap:

  Number of overlapping samples. Default is nperseg/2.

- assay_name:

  Input assay name. If NULL, uses default assay.

## Value

A list with components:

- coherence:

  3D array (freq x channel x channel) of coherence values

- frequencies:

  Frequency vector

- channel_names:

  Channel names

## Details

Calculates the magnitude-squared coherence between pairs of channels,
which measures the linear correlation between signals as a function of
frequency.

Coherence is computed using Welch's averaged periodogram method. Values
range from 0 (no linear relationship) to 1 (perfect linear
relationship).

## References

Nolte, G., et al. (2004). "Identifying true brain interaction from EEG
data using the imaginary part of coherency." Clinical Neurophysiology,
115(10), 2292-2307. doi:10.1016/j.clinph.2004.04.029

## See also

[`crossSpectrum()`](https://x-biosignal.github.io/PhysioAnalysis/reference/crossSpectrum.md)
for the underlying cross-spectral density,
[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md)
for phase-based connectivity,
[`connectivityMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/connectivityMatrix.md)
for a unified connectivity interface.

## Examples

``` r
# Create example with 4 channels
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
  colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
  samplingRate = 256
)

# Compute coherence
coh <- coherence(pe, freq_range = c(1, 50))
dim(coh$coherence)  # frequencies x channels x channels
#> [1] 50  4  4
```
