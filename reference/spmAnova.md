# SPM ANOVA (F-test)

Performs SPM F-test for comparing multiple groups.

## Usage

``` r
spmAnova(x, groups, alpha = 0.05)
```

## Arguments

- x:

  A PhysioExperiment object or matrix (time x observations).

- groups:

  Factor or list indicating group membership.

- alpha:

  Significance level.

## Value

A list of class "spm_result" containing F-statistics and clusters.

## Examples

``` r
# Three-group comparison
set.seed(123)
g1 <- matrix(rnorm(100 * 8), nrow = 100)
g2 <- matrix(rnorm(100 * 8), nrow = 100) + 0.5
g3 <- matrix(rnorm(100 * 8), nrow = 100) + 1.0

data <- cbind(g1, g2, g3)
groups <- factor(rep(c("A", "B", "C"), each = 8))

pe <- PhysioExperiment(assays = list(values = data), samplingRate = 100)
result <- spmAnova(pe, groups = groups)
```
