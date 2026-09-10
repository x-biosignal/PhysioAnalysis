# Kruskal-Wallis rank-sum test (non-parametric between-groups omnibus)

The non-parametric between-groups omnibus test: pools all observations,
ranks them, and tests whether the rank totals differ across `k >= 2`
independent groups. It is the distribution-free counterpart of the
one-way ANOVA
([`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md))
and, unlike it, makes no normality or equal-variance assumption. The
tie-corrected statistic and its p-value reproduce base R
[`stats::kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html)
bit-for-bit (the H is accumulated in the same order) and
`scipy.stats.kruskal` to machine precision.

## Usage

``` r
kruskalTest(x, ...)
```

## Arguments

- x:

  Either a list of `k >= 2` numeric group vectors, or the first group's
  numeric vector (with further groups passed through `...`).

- ...:

  Additional numeric group vectors when `x` is a single vector.

## Value

A list with `statistic` (the tie-corrected H), `p_value`, `df`
(`k - 1`), `k` (number of groups) and `n` (total sample size).

## See also

[`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md),
[`friedmanTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/friedmanTest.md)

## Examples

``` r
kruskalTest(list(rnorm(20), rnorm(20, 1), rnorm(20, 2)))
#> $statistic
#> [1] 27.2118
#> 
#> $p_value
#> [1] 1.233196e-06
#> 
#> $df
#> [1] 2
#> 
#> $k
#> [1] 3
#> 
#> $n
#> [1] 60
#> 
```
