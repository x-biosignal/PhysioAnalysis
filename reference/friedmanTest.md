# Friedman rank-sum test (non-parametric repeated-measures omnibus)

The non-parametric repeated-measures omnibus test: compares `k >= 2`
matched treatments (columns) measured on the same `N` blocks/subjects
(rows), by ranking within each block and testing whether the rank totals
differ. It is the within-subject, distribution-free counterpart of the
one-way ANOVA
([`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md)).
The tie-corrected chi-square statistic and its p-value reproduce base R
[`stats::friedman.test`](https://rdrr.io/r/stats/friedman.test.html)
bit-for-bit and `scipy.stats.friedmanchisquare` to machine precision.

## Usage

``` r
friedmanTest(x, ...)
```

## Arguments

- x:

  A blocks-by-treatments numeric matrix (rows = subjects/blocks, columns
  = conditions/treatments), or a list of `k` equal-length numeric
  vectors (one per treatment, aligned by block). Rows with any missing
  value are dropped (complete blocks only).

- ...:

  Additional treatment vectors when `x` is a single vector.

## Value

A list with `statistic` (the tie-corrected Friedman chi-square),
`p_value`, `df` (`k - 1`), `k` (treatments) and `n` (complete blocks).

## See also

[`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md),
[`pairedTTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/pairedTTest.md)

## Examples

``` r
m <- matrix(rnorm(40) + rep(c(0, 0.5, 1, 1.5), each = 10), nrow = 10)
friedmanTest(m)
#> $statistic
#> [1] 5.16
#> 
#> $p_value
#> [1] 0.1604491
#> 
#> $df
#> [1] 3
#> 
#> $k
#> [1] 4
#> 
#> $n
#> [1] 10
#> 
```
