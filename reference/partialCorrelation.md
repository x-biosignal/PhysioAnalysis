# Partial correlation controlling for covariates

The correlation between `x` and `y` after linearly removing one or more
covariates `z` from both – the standard tool for asking whether an
association survives adjustment for a confound. With
`method = "spearman"` the variables are rank-transformed first (a rank /
robust partial correlation). The coefficient reproduces base R's
residual method and the `ppcor` package to machine precision.

## Usage

``` r
partialCorrelation(x, y = NULL, z = NULL, method = c("pearson", "spearman"))
```

## Arguments

- x:

  Either the first numeric variable (with `y` and `z` given), or a
  matrix/data frame whose column 1 is `x`, column 2 is `y` and the
  remaining columns are the covariates. Incomplete rows are dropped.

- y:

  The second numeric variable when `x` is a single vector.

- z:

  The covariate(s) to control for: a numeric vector or matrix.

- method:

  "pearson" (default) or "spearman" (rank-based).

## Value

A list with `statistic` (the partial correlation), `p_value` (two-sided
t-test), `df` (`n - 2 - k`), `n` and `k` (number of covariates).

## See also

[`correlationTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/correlationTest.md),
[`kendallTau`](https://x-biosignal.github.io/PhysioAnalysis/reference/kendallTau.md)

## Examples

``` r
x <- rnorm(50); z <- rnorm(50); y <- 0.5 * z + rnorm(50)
partialCorrelation(x, y, z)
#> $statistic
#> [1] -0.05034463
#> 
#> $p_value
#> [1] 0.7311977
#> 
#> $df
#> [1] 47
#> 
#> $n
#> [1] 50
#> 
#> $k
#> [1] 1
#> 
```
