# Compute Phase Lag Index (PLI)

Calculates the Phase Lag Index between channel pairs, a measure of
asymmetry in the phase difference distribution.

## Usage

``` r
pli(x, freq_band, channels = NULL, assay_name = NULL)
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

A numeric matrix of PLI values (channel x channel) with values ranging
from 0 (no asymmetry in phase distribution) to 1 (maximal asymmetry).
Row and column names are set to channel names when available.

## Details

PLI measures the asymmetry of the distribution of phase differences.
Unlike PLV, PLI is insensitive to volume conduction effects that lead to
zero-lag synchronization.

PLI = \|mean(sign(Im(S_xy)))\|, where S_xy is the cross-spectrum.

## References

Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in brain
signals." Human Brain Mapping, 8(4), 194-208.
doi:10.1002/(SICI)1097-0193(1999)8:4\<194::AID-HBM4\>3.0.CO;2-C

## See also

[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md)
for phase locking value,
[`wPLI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/wPLI.md)
for weighted phase lag index,
[`coherence()`](https://x-biosignal.github.io/PhysioAnalysis/reference/coherence.md)
for frequency-domain coherence,
[`connectivityMatrix()`](https://x-biosignal.github.io/PhysioAnalysis/reference/connectivityMatrix.md)
for a unified connectivity interface.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(4000), nrow = 1000, ncol = 4)),
  colData = S4Vectors::DataFrame(label = c("Fz", "Cz", "Pz", "Oz")),
  samplingRate = 256
)

# Compute PLI in alpha band (8-12 Hz)
pli_matrix <- pli(pe, freq_band = c(8, 12))
```
