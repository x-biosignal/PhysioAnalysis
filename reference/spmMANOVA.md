# SPM one-way MANOVA (Hotelling's T^2 / chi-square vector field)

Computes a multivariate SPM field for a one-way MANOVA at every node,
for vector-valued waveforms (e.g. 3D joint angles). With two groups the
field is Hotelling's \\T^2\\ (reducing to the squared two-sample t-field
for a single component); with more than two groups it is the chi-square
field of the Bartlett-transformed Wilks' \\\Lambda\\. Field-level
significance uses random-field theory (an F-field threshold for \\T^2\\,
a chi-square-field threshold for the \\X^2\\ field).

## Usage

``` r
spmMANOVA(x, groups, vector_components = NULL, alpha = 0.05)
```

## Arguments

- x:

  A 3D array `time x obs x component`, or a list of component matrices
  (each `time x obs`), or a single `time x obs` matrix (one component).

- groups:

  A grouping factor, one value per observation.

- vector_components:

  Optional integer; only used to disambiguate a 2D matrix input
  (defaults to a single component).

- alpha:

  Significance level (default 0.05).

## Value

A list of class `"spm_result"` (`test_type = "manova"`) with the
statistic field (`T2` or `X2`), RFT `threshold`, significant `clusters`,
`p_values`, degrees of freedom, `fwhm`, and `resel_count`.

## References

Pataky 2016 (vector-field 1D SPM); Worsley 1994 (chi-square RFT).
spm1d.stats.manova1 / hotellings2.

## See also

[`spmTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmTTest.md),
[`spmRegression()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmRegression.md),
[`spmSnPM()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmSnPM.md)

## Examples

``` r
set.seed(1)
# two groups, 3 components (e.g. hip flexion/abduction/rotation)
arr <- array(rnorm(50 * 16 * 3), c(50, 16, 3))
arr[20:30, 1:8, ] <- arr[20:30, 1:8, ] + 1.2
spmMANOVA(arr, groups = rep(c("A", "B"), each = 8))
#> SPM Analysis Result
#> ==================
#> Test type: manova
#> Time points: 50
#> Alpha: 0.050
#> Threshold: 42.471
#> FWHM (smoothness): 2.21
#> Resel count: 22.15
#> 
#> Significant clusters: 2
#>   Cluster 1: [22-22] extent=1, p=0.0056
#>   Cluster 2: [30-30] extent=1, p=0.0056
```
