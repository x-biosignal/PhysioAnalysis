# Compute Phase Locking Value (PLV)

Calculates the Phase Locking Value between channel pairs, measuring the
consistency of phase difference across time.

## Usage

``` r
plv(x, freq_band, channels = NULL, assay_name = NULL)
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

A numeric matrix of PLV values (channel x channel) with values ranging
from 0 (random phase relationship) to 1 (perfect phase locking). Row and
column names are set to channel names when available.

## Details

PLV measures the consistency of the phase difference between two
signals. A value of 1 indicates perfect phase locking, while 0 indicates
random phase relationship.

The signals are first bandpass filtered to the specified frequency band,
then the analytic signal is computed using the Hilbert transform.

## References

Lachaux, J.-P., et al. (1999). "Measuring phase synchrony in brain
signals." Human Brain Mapping, 8(4), 194-208.
doi:10.1002/(SICI)1097-0193(1999)8:4\<194::AID-HBM4\>3.0.CO;2-C

## See also

[`pli()`](https://x-biosignal.github.io/PhysioAnalysis/reference/pli.md)
for phase lag index,
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

# Compute PLV in alpha band (8-12 Hz)
plv_matrix <- plv(pe, freq_band = c(8, 12))
```
