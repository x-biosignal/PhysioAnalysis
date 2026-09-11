# Nemenyi post-hoc test (all pairwise comparisons after a Friedman test)

The Nemenyi test: all pairwise comparisons of the treatments following a
Friedman repeated-measures test
([`friedmanTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/friedmanTest.md)),
from the within-block mean ranks and the studentized-range distribution
with single-step control of the family-wise error rate. It is the
repeated-measures counterpart of Tukey's HSD
([`tukeyHSD`](https://x-biosignal.github.io/PhysioAnalysis/reference/tukeyHSD.md))
and Dunn's test
([`dunnTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/dunnTest.md)).
The statistics and p-values reproduce base R's Nemenyi formula, the
`PMCMRplus` package and `scikit_posthocs.posthoc_nemenyi_friedman` to
machine precision.

## Usage

``` r
nemenyiTest(x)
```

## Arguments

- x:

  A blocks-by-treatments numeric matrix (rows = subjects/blocks, columns
  = conditions/treatments); rows with any missing value are dropped.

## Value

A data frame with one row per pair: `comparison`, `statistic` and
`p_adj`, carrying `k` (treatments) and `n` (complete blocks) as
attributes.

## See also

[`friedmanTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/friedmanTest.md),
[`tukeyHSD`](https://x-biosignal.github.io/PhysioAnalysis/reference/tukeyHSD.md),
[`dunnTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/dunnTest.md)

## Examples

``` r
m <- matrix(rnorm(40) + rep(c(0, 0.5, 1, 1.5), each = 10), nrow = 10)
colnames(m) <- c("a", "b", "c", "d"); nemenyiTest(m)
#>    comparison statistic      p_adj
#> b      b vs a 0.1732051 0.99815453
#> c      c vs a 1.9052559 0.22574088
#> c1     c vs b 2.0784610 0.16012588
#> d      d vs a 2.4248711 0.07245072
#> d1     d vs b 2.5980762 0.04626766
#> d2     d vs c 0.5196152 0.95441134
```
