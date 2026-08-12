# Bootstrap confidence interval for ERP

Computes bootstrap confidence intervals for averaged epochs.

## Usage

``` r
bootstrapCI(
  x,
  n_bootstrap = 1000L,
  ci_level = 0.95,
  condition = NULL,
  seed = NULL
)
```

## Arguments

- x:

  An epoched PhysioExperiment object (4D data).

- n_bootstrap:

  Number of bootstrap iterations (default: 1000).

- ci_level:

  Confidence interval level (default: 0.95).

- condition:

  Epoch indices to include. If NULL, uses all epochs.

- seed:

  Random seed for reproducibility.

## Value

A list containing:

- mean:

  Mean across epochs (time x channel)

- ci_lower:

  Lower CI bound (time x channel)

- ci_upper:

  Upper CI bound (time x channel)

- se:

  Standard error (time x channel)

## References

Maris, E. & Oostenveld, R. (2007). "Nonparametric statistical testing of
EEG- and MEG-data." *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

## See also

[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md)
for Cohen's d effect size,
[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for parametric significance testing,
[`plotERP()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotERP.md)
for ERP visualization.

## Examples

``` r
# Create example epoched data
set.seed(123)
epochs <- array(rnorm(100 * 4 * 20 * 1), dim = c(100, 4, 20, 1))
pe <- PhysioExperiment(
  assays = list(epoched = epochs),
  samplingRate = 100
)
# Bootstrap CI
result <- bootstrapCI(pe, n_bootstrap = 500)
```
