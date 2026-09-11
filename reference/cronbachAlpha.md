# Cronbach's alpha (internal-consistency reliability)

Cronbach's alpha, the standard index of a scale's internal-consistency
reliability: how consistently a set of `k` items measures a single
construct. Computed from the item variances and the total-score
variance. The raw alpha and the mean inter-item correlation reproduce
base R, the `psych` package (`psych::alpha` raw_alpha) and
`pingouin.cronbach_alpha` bit-for-bit.

## Usage

``` r
cronbachAlpha(x)
```

## Arguments

- x:

  A respondents-by-items numeric matrix or data frame (rows =
  respondents, columns = scale items). Rows with any missing value are
  dropped.

## Value

A list with `alpha` (raw Cronbach's alpha), `average_r` (mean inter-item
correlation), `k` (items) and `n` (respondents).

## Examples

``` r
cronbachAlpha(matrix(round(runif(150, 1, 5)), ncol = 5))
#> $alpha
#> [1] 0.1695771
#> 
#> $average_r
#> [1] 0.04097421
#> 
#> $k
#> [1] 5
#> 
#> $n
#> [1] 30
#> 
```
