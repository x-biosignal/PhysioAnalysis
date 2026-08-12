# SPM linear-regression t-field

Computes an SPM{t} field testing, at every node, the slope of a linear
regression of the waveform on a continuous predictor. The node-wise
statistic is the ordinary least-squares slope t-value (identical to
`summary(lm(y ~ predictor))` at that node), and field-level significance
uses the same random-field-theory threshold as
[`spmTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmTTest.md)
with `df = n - 2`.

## Usage

``` r
spmRegression(x, predictor, alpha = 0.05, two_tailed = TRUE)
```

## Arguments

- x:

  A `PhysioExperiment` or a numeric matrix (time x observations).

- predictor:

  Numeric covariate, one value per observation (column).

- alpha:

  Significance level (default 0.05).

- two_tailed:

  Logical; two-tailed slope test (default `TRUE`).

## Value

A list of class `"spm_result"` (`test_type = "regression"`) with the `t`
field, RFT `threshold`, significant `clusters`, pointwise `p_values`,
`df`, `fwhm`, `resel_count`, and the fitted `slope`/`intercept` fields.

## References

Pataky 2016; Friston et al. 2007 (RFT). spm1d.stats.regress.

## See also

[`spmTTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmTTest.md),
[`spmMANOVA()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmMANOVA.md),
[`spmSnPM()`](https://x-biosignal.github.io/PhysioAnalysis/reference/spmSnPM.md)

## Examples

``` r
set.seed(1)
pred <- rnorm(20)
data <- matrix(rnorm(100 * 20), nrow = 100)
data[40:60, ] <- data[40:60, ] + outer(rep(1.5, 21), pred)  # slope 30% window
spmRegression(data, pred)
#> SPM Analysis Result
#> ==================
#> Test type: regression
#> Time points: 100
#> Alpha: 0.050
#> Threshold: 4.374
#> FWHM (smoothness): 2.24
#> Resel count: 44.12
#> 
#> Significant clusters: 4
#>   Cluster 1: [40-40] extent=1, p=0.0114
#>   Cluster 2: [42-44] extent=3, p=0.0000
#>   Cluster 3: [46-57] extent=12, p=0.0000
#>   Cluster 4: [59-60] extent=2, p=0.0001
```
