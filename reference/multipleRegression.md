# Multiple linear regression (ordinary least squares)

Multivariable OLS regression of a response on two or more predictors:
the coefficient table (estimate, standard error, t-statistic, p-value)
plus the overall model fit (R-squared, adjusted R-squared and the
F-test). The many-predictor generalisation of
[`linearRegression`](https://x-biosignal.github.io/PhysioAnalysis/reference/linearRegression.md).
Computed via the same QR decomposition as base R, so the coefficients,
R-squared and F reproduce [`stats::lm`](https://rdrr.io/r/stats/lm.html)
/ `summary.lm` bit-for-bit and `statsmodels` OLS to machine precision.

## Usage

``` r
multipleRegression(x)
```

## Arguments

- x:

  A matrix or data frame whose first column is the response and whose
  remaining columns are the predictors. Incomplete rows are dropped.

## Value

A data frame with one row per term (`term`, `estimate`, `std_error`,
`statistic`, `p_value`), carrying `r_squared`, `adj_r_squared`,
`f_statistic`, `df1`, `df2`, `f_p_value` and `n` as attributes.

## See also

[`linearRegression`](https://x-biosignal.github.io/PhysioAnalysis/reference/linearRegression.md),
[`partialCorrelation`](https://x-biosignal.github.io/PhysioAnalysis/reference/partialCorrelation.md)

## Examples

``` r
x <- matrix(rnorm(120), ncol = 3); colnames(x) <- c("y", "a", "b")
multipleRegression(x)
#>          term     estimate std_error   statistic   p_value
#> 1 (Intercept) -0.107750779 0.1605314 -0.67121305 0.5062539
#> 2           a  0.008988436 0.1431897  0.06277294 0.9502852
#> 3           b -0.004637306 0.1866392 -0.02484637 0.9803110
```
