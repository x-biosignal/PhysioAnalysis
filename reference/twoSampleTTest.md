# Two-Sample t-Test

Computes the two-sided independent two-sample t-test comparing two
groups. By default it is Welch's t-test (unequal variances,
Welch-Satterthwaite df) – matching
[`stats::t.test`](https://rdrr.io/r/stats/t.test.html) and
`scipy.stats.ttest_ind(equal_var = FALSE)`; set `var_equal = TRUE` for
Student's pooled-variance test. The standard test for a difference
between two independent groups (e.g. older vs younger, patients vs
controls). Reproduces both references bit-for-bit.

## Usage

``` r
twoSampleTTest(x, y = NULL, var_equal = FALSE)
```

## Arguments

- x:

  A numeric vector (group 1), OR (when `y` is `NULL`) a list of two
  numeric vectors.

- y:

  A numeric vector (group 2), or `NULL` if `x` is a two-element list.

- var_equal:

  Logical; `FALSE` (default) for Welch's t-test, `TRUE` for Student's
  pooled-variance t-test.

## Value

A list with `statistic` (the t-statistic of \\x - y\\), `p_value`
(two-sided) and `df`.

## See also

[`pairedTTest`](https://x-biosignal.github.io/PhysioAnalysis/reference/pairedTTest.md)
for matched samples.

## Examples

``` r
set.seed(1); twoSampleTTest(rnorm(20, 1), rnorm(30))$statistic  # group 1 higher -> positive t
#> [1] 4.566756
```
