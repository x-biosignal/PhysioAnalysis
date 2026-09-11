# Dunn's test (post-hoc pairwise comparisons after Kruskal-Wallis)

Dunn's test: all pairwise rank-based comparisons following a
Kruskal-Wallis test
([`kruskalTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/kruskalTest.md)),
using the shared tie-corrected rank variance and a chosen
multiple-comparison adjustment. It is the non-parametric counterpart of
Tukey's HSD
([`tukeyHSD`](https://x-biosignal.github.io/PhysioAnalysis/reference/tukeyHSD.md)).
The z-statistics and adjusted p-values reproduce base R's Dunn formula,
the `PMCMRplus` package and `scikit_posthocs.posthoc_dunn` to machine
precision.

## Usage

``` r
dunnTest(x, ..., p_adjust = c("bonferroni", "none", "holm", "BH"))
```

## Arguments

- x:

  Either a named list of `k >= 2` numeric group vectors, or the first
  group's numeric vector (with further groups passed through `...`).

- ...:

  Additional numeric group vectors when `x` is a single vector.

- p_adjust:

  Multiple-comparison adjustment: "bonferroni" (default), "none", "holm"
  or "BH".

## Value

A data frame with one row per pair: `comparison`, `z`, `p_value`
(two-sided, unadjusted) and `p_adj`, carrying `k`, `n` as attributes.

## See also

[`kruskalTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/kruskalTest.md),
[`tukeyHSD`](https://x-biosignal.github.io/PhysioAnalysis/reference/tukeyHSD.md)

## Examples

``` r
dunnTest(list(a = rnorm(20), b = rnorm(20, 1), c = rnorm(20, 2)))
#>    comparison        z      p_value        p_adj
#> b         b-a 3.032947 2.421777e-03 7.265332e-03
#> c         c-a 5.359716 8.335280e-08 2.500584e-07
#> c1        c-b 2.326769 1.997758e-02 5.993274e-02
```
