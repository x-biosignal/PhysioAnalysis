# SPM paired t-test

Performs SPM analysis for paired/repeated measures data.

## Usage

``` r
spmPairedTTest(x, condition1, condition2, alpha = 0.05, two_tailed = TRUE)
```

## Arguments

- x:

  A PhysioExperiment object or matrix (time x observations).

- condition1:

  Indices for first condition.

- condition2:

  Indices for second condition.

- alpha:

  Significance level.

- two_tailed:

  Logical; if TRUE, performs two-tailed test.

## Value

A list of class "spm_result".

## Examples

``` r
# Pre-post intervention comparison
set.seed(123)
pre <- matrix(rnorm(100 * 15), nrow = 100)
post <- pre + 0.8  # Effect across all time points
post[30:50, ] <- post[30:50, ] + 0.5  # Additional effect

data <- cbind(pre, post)
pe <- PhysioExperiment(assays = list(values = data), samplingRate = 100)

result <- spmPairedTTest(pe, condition1 = 1:15, condition2 = 16:30)
```
