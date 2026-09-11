# Tukey HSD post-hoc test (all pairwise comparisons after one-way ANOVA)

Tukey's Honest Significant Difference test: all pairwise mean
comparisons following a one-way ANOVA
([`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md)),
with p-values and confidence intervals from the studentized-range
distribution and single-step control of the family-wise error rate.
Reproduces base R
[`stats::TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html) and
`scipy.stats.tukey_hsd` (machine precision).

## Usage

``` r
tukeyHSD(x, ..., conf_level = 0.95)
```

## Arguments

- x:

  Either a named list of `k >= 2` numeric group vectors, or the first
  group's numeric vector (with further groups passed through `...` – use
  a named list to label the comparisons).

- ...:

  Additional numeric group vectors when `x` is a single vector.

- conf_level:

  Confidence level for the intervals (default 0.95).

## Value

A data frame with one row per pair: `comparison`, `diff` (mean
difference), `lwr`/`upr` (CI) and `p_adj` (family-wise-adjusted
p-value), carrying `k`, `df` and `n` as attributes.

## See also

[`oneWayAnova`](https://x-biosignal.github.io/PhysioAnalysis/reference/oneWayAnova.md)

## Examples

``` r
tukeyHSD(list(a = rnorm(20), b = rnorm(20, 1), c = rnorm(20, 2)))
#>    comparison      diff        lwr      upr        p_adj
#> b         b-a 0.5290554 -0.2624564 1.320567 2.503816e-01
#> c         c-a 1.6022073  0.8106955 2.393719 2.715375e-05
#> c1        c-b 1.0731519  0.2816402 1.864664 5.228321e-03
```
