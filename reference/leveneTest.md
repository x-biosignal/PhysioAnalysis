# Levene's test for homogeneity of variance

Levene's test of whether `k >= 2` groups have equal variance, as a
one-way ANOVA on the absolute deviations from each group's centre. With
`center = "median"` (default) it is the robust Brown-Forsythe variant,
appropriate when the data may be non-normal; `center = "mean"` is the
classic Levene test. The F-statistic and p reproduce
`scipy.stats.levene` and the `car` package (`car::leveneTest`)
bit-for-bit.

## Usage

``` r
leveneTest(x, ..., center = c("median", "mean"))
```

## Arguments

- x:

  Either a list of `k >= 2` numeric group vectors, or the first group's
  numeric vector (with further groups passed through `...`).

- ...:

  Additional numeric group vectors when `x` is a single vector.

- center:

  "median" (default, Brown-Forsythe) or "mean" (classic Levene).

## Value

A list with `statistic` (F), `p_value`, `df1` (`k - 1`), `df2`
(`N - k`), `k`, `n` and `center`.

## See also

[`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md)

## Examples

``` r
leveneTest(list(rnorm(20), rnorm(20, 0, 2), rnorm(20, 0, 0.5)))
#> $statistic
#> [1] 16.09617
#> 
#> $p_value
#> [1] 2.871424e-06
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
#> $center
#> [1] "median"
#> 
```
