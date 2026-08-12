# Statistical Testing for PhysioExperiment

Functions for statistical analysis of physiological signal data,
including t-tests, ANOVA, cluster-based permutation tests, and effect
sizes. Pointwise t-test across epochs

## Usage

``` r
tTestEpochs(
  x,
  condition1 = NULL,
  condition2 = NULL,
  mu = 0,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  var.equal = FALSE
)
```

## Arguments

- x:

  An epoched PhysioExperiment object (4D data).

- condition1:

  Indices or logical vector for first condition epochs.

- condition2:

  Indices or logical vector for second condition epochs. If NULL,
  performs one-sample t-test against mu.

- mu:

  Value to test against for one-sample t-test (default: 0).

- paired:

  Logical; if TRUE, performs paired t-test.

- alternative:

  Alternative hypothesis: "two.sided", "less", or "greater".

- var.equal:

  Logical; if TRUE, assumes equal variances.

## Value

A list containing:

- t_values:

  Matrix of t-statistics (time x channel)

- p_values:

  Matrix of p-values (time x channel)

- df:

  Degrees of freedom

- n1, n2:

  Sample sizes for each condition

## Details

Performs t-tests at each time point and channel, comparing epochs
against a baseline or between two conditions.

## References

Maris, E. & Oostenveld, R. (2007). "Nonparametric statistical testing of
EEG- and MEG-data." *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

## See also

[`anovaEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/anovaEpochs.md)
for multi-group comparisons,
[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md)
for multiple comparison correction,
[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md)
for Cohen's d effect size,
[`correctPValues()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correctPValues.md)
for p-value correction methods.

## Examples

``` r
# Create example epoched data
set.seed(123)
epochs <- array(rnorm(100 * 4 * 20 * 1), dim = c(100, 4, 20, 1))
pe <- PhysioExperiment(
  assays = list(epoched = epochs),
  samplingRate = 100
)
# One-sample t-test against zero
result <- tTestEpochs(pe)
# Two-sample t-test comparing conditions
result2 <- tTestEpochs(pe, condition1 = 1:10, condition2 = 11:20)
```
