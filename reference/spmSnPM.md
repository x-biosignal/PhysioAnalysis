# Non-parametric (permutation) SPM

Builds the permutation distribution of the field-maximum statistic to
obtain a distribution-free critical threshold with strong FWER control
(Nichols & Holmes 2002), plus permutation cluster- and set-level
p-values. This is the validity cross-check for the parametric RFT
thresholds: as the field smoothness grows the permutation threshold
converges to the RFT threshold.

## Usage

``` r
spmSnPM(
  x,
  groups = NULL,
  group1 = NULL,
  group2 = NULL,
  statistic = c("t", "F", "T2"),
  n_permutations = 1000,
  alpha = 0.05,
  two_tailed = TRUE,
  vector_components = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment`/matrix (time x obs) for `t`/`F`, or a 3D array /
  list of component matrices for `T2`.

- groups:

  Grouping factor (one per observation) for `F` or `T2`, or a two-sample
  `t`.

- group1, group2:

  Observation indices for a two-sample `t` (`group2 = NULL` gives a
  one-sample sign-flip test).

- statistic:

  `"t"`, `"F"`, or `"T2"`.

- n_permutations:

  Number of random permutations (default 1000).

- alpha:

  Significance level (default 0.05).

- two_tailed:

  For `t`: use the two-tailed field maximum `max|t|`.

- vector_components:

  Passed to the `T2` array coercion.

- seed:

  Optional RNG seed for reproducibility.

## Value

A list of class `"spm_result"` (`test_type = "snpm"`) with the observed
statistic field, the permutation `threshold`, `clusters` with
permutation p-values, pointwise permutation `p_values`, and the
`perm_max` distribution.

## References

Nichols & Holmes 2002; Pataky 2016. spm1d.stats.nonparam.

## See also

[`spmTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmTTest.md),
[`spmAnova()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmAnova.md),
[`spmMANOVA()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmMANOVA.md)

## Examples

``` r
set.seed(1)
g1 <- matrix(rnorm(50 * 10), 50); g2 <- matrix(rnorm(50 * 10), 50)
g2[20:30, ] <- g2[20:30, ] + 1.5
spmSnPM(cbind(g1, g2), group1 = 1:10, group2 = 11:20,
        n_permutations = 200, seed = 1)
#> SPM Analysis Result
#> ==================
#> Test type: snpm
#> Time points: 50
#> Alpha: 0.050
#> Threshold: 4.125
#> Permutations: 200
#> 
#> Significant clusters: 2
#>   Cluster 1: [23-23] extent=1, p=0.0050
#>   Cluster 2: [28-28] extent=1, p=0.0050
```
