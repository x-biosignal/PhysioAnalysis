# ANOVA across conditions

Performs one-way ANOVA at each time point and channel across multiple
conditions.

## Usage

``` r
anovaEpochs(x, groups)
```

## Arguments

- x:

  An epoched PhysioExperiment object (4D data).

- groups:

  Factor or character vector indicating group membership for each epoch.
  Can also be a column name from epoch_info metadata.

## Value

A list containing:

- f_values:

  Matrix of F-statistics (time x channel)

- p_values:

  Matrix of p-values (time x channel)

- df_between:

  Between-group degrees of freedom

- df_within:

  Within-group degrees of freedom

- group_means:

  Array of group means (time x channel x group)

## References

Maris, E. & Oostenveld, R. (2007). "Nonparametric statistical testing of
EEG- and MEG-data." *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

## See also

[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for pairwise comparisons,
[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md)
for multiple comparison correction,
[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md)
for Cohen's d effect size.

## Examples

``` r
# Create example epoched data with conditions
set.seed(123)
epochs <- array(rnorm(100 * 4 * 30 * 1), dim = c(100, 4, 30, 1))
pe <- PhysioExperiment(
  assays = list(epoched = epochs),
  samplingRate = 100,
  metadata = list(
    epoch_info = S4Vectors::DataFrame(
      condition = rep(c("A", "B", "C"), each = 10)
    )
  )
)
# ANOVA across conditions
result <- anovaEpochs(pe, groups = "condition")
```
