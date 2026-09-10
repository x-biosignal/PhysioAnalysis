# One-way ANOVA (omnibus F-test) across groups

Classic one-way analysis of variance (equal-variance omnibus F-test)
comparing the means of `k >= 2` independent groups. The F-statistic and
its p-value reproduce base R `stats::oneway.test(var.equal = TRUE)` and
`scipy.stats.f_oneway` bit-for-bit. This is the standard many-group
generalisation of the independent two-sample t-test
([`twoSampleTTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/twoSampleTTest.md)).

## Usage

``` r
oneWayAnova(x, ...)
```

## Arguments

- x:

  Either a list of `k >= 2` numeric group vectors, or the first group's
  numeric vector (with further groups passed through `...`).

- ...:

  Additional numeric group vectors when `x` is a single vector.

## Value

A list with `statistic` (the F-ratio), `p_value`, `df1` (between-groups
degrees of freedom, `k - 1`), `df2` (within-groups, `N - k`), `k`
(number of groups) and `n` (total sample size).

## See also

[`twoSampleTTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/twoSampleTTest.md),
[`pairedTTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/pairedTTest.md)

## Examples

``` r
oneWayAnova(list(rnorm(20), rnorm(20, 1), rnorm(20, 2)))
#> $statistic
#> [1] 37.63678
#> 
#> $p_value
#> [1] 3.805846e-11
#> 
#> $df1
#> [1] 2
#> 
#> $df2
#> [1] 57
#> 
#> $k
#> [1] 3
#> 
#> $n
#> [1] 60
#> 
```
