# Compute weighted Phase Lag Index (wPLI)

Calculates the weighted Phase Lag Index, which is less sensitive to
noise than standard PLI.

## Usage

``` r
wPLI(x, freq_band, channels = NULL, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object.

- freq_band:

  Numeric vector of length 2 specifying frequency band (Hz).

- channels:

  Integer vector of channel indices.

- assay_name:

  Input assay name.

## Value

A numeric matrix of wPLI values (channel x channel) with values ranging
from 0 to 1. Row and column names are set to channel names when
available.

## Details

wPLI weights the contribution of each phase difference by the magnitude
of the imaginary component, reducing the influence of noise sources.

## References

Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in brain
signals." Human Brain Mapping, 8(4), 194-208.
doi:10.1002/(SICI)1097-0193(1999)8:4\<194::AID-HBM4\>3.0.CO;2-C

## See also

[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md)
for phase locking value,
[`pli()`](https://x-biosignal.github.io/PhysioAnalysis/reference/pli.md)
for phase lag index,
[`coherence()`](https://x-biosignal.github.io/PhysioAnalysis/reference/coherence.md)
for frequency-domain coherence,
[`connectivityMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/connectivityMatrix.md)
for a unified connectivity interface.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
  samplingRate = 256
)

# Compute wPLI in theta band (4-8 Hz)
wpli_matrix <- wPLI(pe, freq_band = c(4, 8))
```
