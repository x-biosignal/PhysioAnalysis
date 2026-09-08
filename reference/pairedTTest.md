# Paired-Samples t-Test

Computes the two-sided paired-samples t-test comparing two matched
vectors – the one-sample t-test of the within-pair differences \\d = x -
y\\: \\t = \bar{d} / (s_d / \sqrt{n})\\ on \\n - 1\\ degrees of freedom.
The standard test for before/after or condition-A/condition-B
measurements on the same subjects. Reproduces
`stats::t.test(x, y, paired = TRUE)` and `scipy.stats.ttest_rel`
bit-for-bit.

## Usage

``` r
pairedTTest(x, y = NULL)
```

## Arguments

- x:

  A numeric vector, OR (when `y` is `NULL`) a two-column matrix whose
  columns are the paired samples.

- y:

  A numeric vector matched to `x`, or `NULL` if `x` is a two-column
  matrix.

## Value

A list with `statistic` (the t-statistic of \\x - y\\), `p_value`
(two-sided) and `df`.

## See also

[`signalMoments`](https://x-biosignal.github.io/PhysioAnalysis/reference/signalMoments.md),
[`jarqueBeraTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/jarqueBeraTest.md).

## Examples

``` r
set.seed(1); pairedTTest(rnorm(20, 1), rnorm(20))$statistic  # x shifted up -> positive t
#> [1] 2.578616
```
