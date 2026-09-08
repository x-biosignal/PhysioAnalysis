# Correlation Test

Tests the association between two matched numeric vectors. Returns the
correlation coefficient, its two-sided p-value (from the exact t-test on
the Pearson/Spearman correlation) and the degrees of freedom. For
`method = "pearson"` it reproduces
[`stats::cor.test`](https://rdrr.io/r/stats/cor.test.html) and
`scipy.stats.pearsonr` bit-for-bit.

## Usage

``` r
correlationTest(x, y = NULL, method = c("pearson", "spearman"))
```

## Arguments

- x:

  A numeric vector, OR (when `y` is `NULL`) a matrix whose first two
  columns are the paired variables.

- y:

  A numeric vector matched to `x`, or `NULL` if `x` is a matrix.

- method:

  Correlation method: `"pearson"` (default) or `"spearman"` (rank
  correlation; the t-based p is the large-sample approximation).

## Value

A list with `statistic` (the correlation coefficient), `p_value`
(two-sided) and `df` (\\n - 2\\).

## See also

[`pairedTTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/pairedTTest.md),
[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md).

## Examples

``` r
set.seed(1); x <- rnorm(50); correlationTest(x, x + rnorm(50))$statistic  # ~0.7
#> [1] 0.6339331
```
