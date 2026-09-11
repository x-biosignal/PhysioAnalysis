# Simple linear regression (ordinary least squares)

Univariate OLS regression of `y` on `x`: the slope, intercept,
R-squared, and the two-sided t-test of the slope. The slope, R-squared
and p-value reproduce base R
[`stats::lm`](https://rdrr.io/r/stats/lm.html) / `summary.lm` and
`scipy.stats.linregress` bit-for-bit (to machine precision against
`lm`'s QR path).

## Usage

``` r
linearRegression(x, y = NULL)
```

## Arguments

- x:

  Either the numeric predictor vector (with `y` the response), or a
  two-column matrix/data frame (column 1 = `x`, column 2 = `y`), or a
  list of two numeric vectors. Incomplete `(x, y)` pairs are dropped.

- y:

  The numeric response vector when `x` is a single predictor vector.

## Value

A list with `slope`, `intercept`, `r_squared`, `statistic` (the slope
t-statistic), `p_value` (two-sided), `df` (`n - 2`) and `n`.

## See also

[`correlationTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/correlationTest.md)

## Examples

``` r
linearRegression(1:20, 2 * (1:20) + rnorm(20))
#> $slope
#> [1] 2.065174
#> 
#> $intercept
#> [1] -0.7632061
#> 
#> $r_squared
#> [1] 0.9941564
#> 
#> $statistic
#> [1] 55.33819
#> 
#> $p_value
#> [1] 1.477463e-21
#> 
#> $df
#> [1] 18
#> 
#> $n
#> [1] 20
#> 
```
