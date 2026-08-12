# Rank-biserial correlation for two independent samples

The rank-biserial correlation derived from the Mann-Whitney U statistic,
equal to `2*U1/(n1*n2) - 1` and to `1 - 2*U2/(n1*n2)` (Kerby 2014),
where `U1` counts the pairs in which the first sample exceeds the
second. For two independent samples this equals Cliff's delta; the
confidence interval is Cliff's (1993) consistent-variance interval.

## Usage

``` r
rankBiserial(x, y, conf_level = 0.95)
```

## Arguments

- x, y:

  Numeric vectors for the two samples.

- conf_level:

  Confidence level (default 0.95).

## Value

A list with `r` (rank-biserial correlation), `ci_lower`, `ci_upper`, the
Mann-Whitney `u` for `x`, and the sample sizes.

## References

Kerby, D.S. (2014). Compr Psychol 3:11.IT.3.1. Cliff, N. (1993). Psychol
Bull 114(3):494-509.

## See also

[`cliffsDelta()`](https://x-biosignal.github.io/PhysioAnalysis/reference/cliffsDelta.md),
[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md)

## Examples

``` r
rankBiserial(c(5, 6, 7, 8), c(1, 2, 3, 9))
#> $r
#> [1] 0.5
#> 
#> $ci_lower
#> [1] -0.4878146
#> 
#> $ci_upper
#> [1] 0.9263176
#> 
#> $u
#> [1] 12
#> 
#> $n1
#> [1] 4
#> 
#> $n2
#> [1] 4
#> 
#> $conf_level
#> [1] 0.95
#> 
```
