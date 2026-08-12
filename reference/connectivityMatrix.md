# Compute connectivity matrix for a frequency band

High-level function to compute connectivity using various metrics.

## Usage

``` r
connectivityMatrix(
  x,
  method = c("coherence", "plv", "pli", "wpli", "correlation"),
  freq_band = NULL,
  channels = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- method:

  Connectivity method: "coherence", "plv", "pli", "wpli", "correlation".

- freq_band:

  Frequency band for phase-based methods.

- channels:

  Integer vector of channel indices.

- assay_name:

  Input assay name.

## Value

A numeric matrix (channel x channel) of connectivity values. For
`"coherence"`, returns the mean coherence across frequencies. For
phase-based methods (`"plv"`, `"pli"`, `"wpli"`), returns values between
0 and 1. For `"correlation"`, returns values between -1 and 1.

## References

Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in brain
signals." Human Brain Mapping, 8(4), 194-208.
doi:10.1002/(SICI)1097-0193(1999)8:4\<194::AID-HBM4\>3.0.CO;2-C

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioAnalysis/reference/coherence.md),
[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md),
[`pli()`](https://x-biosignal.github.io/PhysioAnalysis/reference/pli.md),
[`wPLI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/wPLI.md),
[`correlationMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correlationMatrix.md)
for individual connectivity methods.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
  samplingRate = 256
)

# Compute PLV connectivity in alpha band
conn <- connectivityMatrix(pe, method = "plv", freq_band = c(8, 12))
```
