# Intraclass correlation coefficients (ICC)

The six McGraw & Wong (1996) intraclass correlation coefficients from a
subjects-by-raters (or subjects-by-repeated-measurements) matrix: ICC1,
ICC2, ICC3 (single measurement) and ICC1k, ICC2k, ICC3k (average of `k`
measurements). The standard index of test-retest / inter-rater
reliability, computed from a two-way ANOVA. Reproduces the `psych`
package (`psych::ICC(lmer = FALSE)`) to machine precision and
`pingouin.intraclass_corr` bit-for-bit.

## Usage

``` r
intraclassCorrelation(x)
```

## Arguments

- x:

  A subjects-by-raters numeric matrix (rows = targets/subjects, columns
  = raters or repeated measurements). Rows with any missing value are
  dropped.

## Value

A data frame with one row per ICC form (`form`, `ICC`), carrying `n`
(subjects) and `k` (raters/measurements) as attributes.

## Examples

``` r
m <- matrix(rnorm(30), ncol = 3) + rnorm(10)
intraclassCorrelation(m)
#>    form       ICC
#> 1  ICC1 0.6073458
#> 2  ICC2 0.6028622
#> 3  ICC3 0.5828942
#> 4 ICC1k 0.8227048
#> 5 ICC2k 0.8199513
#> 6 ICC3k 0.8074116
```
