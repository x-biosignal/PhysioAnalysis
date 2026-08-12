# Compute effect size (Cohen's d / Hedges g)

Calculates the standardized mean difference at each time point and
channel, with an exact non-central-t confidence interval (Cumming 2014)
and an optional Hedges small-sample bias correction.

## Usage

``` r
effectSize(
  x,
  condition1 = NULL,
  condition2 = NULL,
  pooled = TRUE,
  correction = c("none", "hedges"),
  conf_level = 0.95
)
```

## Arguments

- x:

  An epoched PhysioExperiment object (4D data).

- condition1:

  Indices for first condition.

- condition2:

  Indices for second condition. If NULL, computes d against zero.

- pooled:

  If TRUE (default for two-sample), uses pooled standard deviation.

- correction:

  `"none"` for Cohen's d (default), or `"hedges"` to apply the Hedges
  bias correction `J = 1 - 3/(4*df - 1)` to the estimate and interval.

- conf_level:

  Confidence level for the interval (default 0.95).

## Value

A list containing:

- d:

  Matrix of standardized mean differences (time x channel)

- ci_lower:

  Lower confidence limit for d (non-central t)

- ci_upper:

  Upper confidence limit for d (non-central t)

- correction:

  The bias correction applied

## References

Cumming, G. (2014). The New Statistics. Psychol Sci 25(1):7-29. Hedges,
L.V. (1981). J Educ Stat 6(2):107-128.

## See also

[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for significance testing,
[`bootstrapCI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bootstrapCI.md)
for bootstrap confidence intervals,
[`rankBiserial()`](https://x-biosignal.github.io/PhysioAnalysis/reference/rankBiserial.md)
and
[`cliffsDelta()`](https://x-biosignal.github.io/PhysioAnalysis/reference/cliffsDelta.md)
for rank-based effect sizes.

## Examples

``` r
# Create example epoched data
set.seed(123)
epochs <- array(rnorm(100 * 4 * 20 * 1), dim = c(100, 4, 20, 1))
pe <- PhysioExperiment(
  assays = list(epoched = epochs),
  samplingRate = 100
)
# Effect size for one-sample
result <- effectSize(pe, condition1 = 1:10)
# Hedges g between conditions
result2 <- effectSize(pe, condition1 = 1:10, condition2 = 11:20,
                      correction = "hedges")
```
