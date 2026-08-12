# Compute correlation matrix between channels

Calculates the Pearson correlation coefficient between all channel
pairs.

## Usage

``` r
correlationMatrix(
  x,
  channels = NULL,
  method = c("pearson", "spearman", "kendall"),
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- channels:

  Integer vector of channel indices.

- method:

  Correlation method: "pearson", "spearman", or "kendall".

- assay_name:

  Input assay name.

## Value

A numeric correlation matrix (channel x channel) with values ranging
from -1 to 1. Row and column names are set to channel names when
available.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioAnalysis/reference/coherence.md)
for frequency-domain coherence,
[`plv()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plv.md)
for phase-based connectivity,
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

# Compute Pearson correlation
cor_matrix <- correlationMatrix(pe)
```
