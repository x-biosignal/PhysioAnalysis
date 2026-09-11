# Kendall's tau-b rank correlation

Kendall's tau-b, the concordance-based rank correlation: the normalised
difference between concordant and discordant pairs, tie-corrected. It is
the third rank-correlation method alongside Pearson's r
([`correlationTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/correlationTest.md),
method "pearson") and Spearman's rho (method "spearman"), and is the
robust choice for ordinal or heavily-tied data and small samples. The
tau-b statistic and the normal-approximation p-value reproduce base R
`stats::cor.test(method = "kendall")` and `scipy.stats.kendalltau`
bit-for-bit.

## Usage

``` r
kendallTau(x, y = NULL)
```

## Arguments

- x:

  Either the numeric predictor vector (with `y` the response), or a
  two-column matrix/data frame (column 1 = `x`, column 2 = `y`), or a
  list of two numeric vectors. Incomplete `(x, y)` pairs are dropped.

- y:

  The second numeric vector when `x` is a single vector.

## Value

A list with `statistic` (tau-b), `p_value` (two-sided,
normal-approximation), `z` (the approximation's z-statistic) and `n`.

## See also

[`correlationTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/correlationTest.md)

## Examples

``` r
kendallTau(c(1, 2, 3, 4, 5), c(2, 1, 4, 3, 5))
#> $statistic
#> [1] 0.6
#> 
#> $p_value
#> [1] 0.1416447
#> 
#> $z
#> [1] 1.469694
#> 
#> $n
#> [1] 5
#> 
```
