# Statistical Testing with PhysioAnalysis

## Introduction

PhysioAnalysis provides a suite of statistical testing functions
designed for epoched physiological signal data. This vignette covers:

1.  **Parametric tests** – pointwise t-tests and ANOVA across epochs
2.  **Cluster-based permutation tests** – non-parametric correction for
    multiple comparisons
3.  **Effect sizes and confidence intervals** – Cohen’s d and bootstrap
    CI
4.  **Multiple comparison correction** – FDR, Bonferroni, and Holm
    methods

All statistical functions operate on 4D epoched data (time x channel x
epoch x sample) produced by
[`epochData()`](https://x-biosignal.github.io/PhysioAnalysis/reference/epochData.md).

## Setup

``` r

library(PhysioCore)
library(PhysioAnalysis)
```

## Creating Epoched Test Data

We create synthetic epoched data with a known effect to demonstrate the
statistical testing functions.

``` r

set.seed(123)
n_time <- 100      # time points per epoch
n_channels <- 4    # number of channels
n_epochs <- 40     # total epochs (20 per condition)
sr <- 250          # sampling rate

# Create 4D array: time x channel x epoch x sample
epochs <- array(rnorm(n_time * n_channels * n_epochs * 1),
                dim = c(n_time, n_channels, n_epochs, 1))

# Add a signal to condition 2 (epochs 21-40) in channels 1-2,
# time points 30-60 (simulating an ERP effect)
epochs[30:60, 1:2, 21:40, 1] <- epochs[30:60, 1:2, 21:40, 1] + 1.0

pe <- PhysioExperiment(
  assays = list(epoched = epochs),
  colData = S4Vectors::DataFrame(
    label = paste0("Ch", 1:n_channels)
  ),
  samplingRate = sr,
  metadata = list(
    epoch_tmin = -0.2,
    epoch_tmax = 0.2,
    epoch_info = S4Vectors::DataFrame(
      condition = rep(c("control", "stimulus"), each = 20)
    )
  )
)
```

## Pointwise t-Tests

The
[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
function performs t-tests at each time point and channel. It supports
one-sample, two-sample (independent), and paired designs.

### One-Sample t-Test

Test whether epoch values differ significantly from zero.

``` r

# One-sample t-test against zero (all epochs)
result_one <- tTestEpochs(pe)

# The result contains t-values and p-values as time x channel matrices
str(result_one)
dim(result_one$t_values)  # 100 x 4
dim(result_one$p_values)  # 100 x 4
```

### Two-Sample t-Test

Compare two conditions (e.g., control vs. stimulus).

``` r

# Two-sample t-test: control (epochs 1-20) vs stimulus (epochs 21-40)
result_two <- tTestEpochs(pe,
                          condition1 = 1:20,
                          condition2 = 21:40)

# Find time points with significant effects (uncorrected)
sig_mask <- result_two$p_values < 0.05
cat("Number of significant time-channel pairs:", sum(sig_mask, na.rm = TRUE), "\n")
```

### Paired t-Test

When epochs are naturally paired (e.g., pre/post within subjects).

``` r

result_paired <- tTestEpochs(pe,
                             condition1 = 1:20,
                             condition2 = 21:40,
                             paired = TRUE)
```

### One-Sided Tests

``` r

# Test if stimulus > control
result_greater <- tTestEpochs(pe,
                              condition1 = 21:40,
                              condition2 = 1:20,
                              alternative = "greater")
```

## ANOVA Across Conditions

The
[`anovaEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/anovaEpochs.md)
function performs one-way ANOVA at each time point and channel across
multiple conditions.

``` r

# Create data with 3 conditions
set.seed(42)
epochs3 <- array(rnorm(n_time * n_channels * 30 * 1),
                 dim = c(n_time, n_channels, 30, 1))

# Add effects to condition B and C
epochs3[30:60, 1:2, 11:20, 1] <- epochs3[30:60, 1:2, 11:20, 1] + 0.8
epochs3[30:60, 1:2, 21:30, 1] <- epochs3[30:60, 1:2, 21:30, 1] + 1.5

pe3 <- PhysioExperiment(
  assays = list(epoched = epochs3),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:n_channels)),
  samplingRate = sr,
  metadata = list(
    epoch_info = S4Vectors::DataFrame(
      condition = rep(c("A", "B", "C"), each = 10)
    )
  )
)

# Run ANOVA using the condition column from epoch_info
anova_result <- anovaEpochs(pe3, groups = "condition")

# Results include F-values, p-values, and group means
str(anova_result)
```

## Cluster-Based Permutation Tests

The
[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md)
function implements the non-parametric cluster-based permutation test
described by Maris and Oostenveld (2007). This approach corrects for the
multiple comparisons problem while maintaining sensitivity to effects
that are distributed across time and channels.

### How It Works

1.  Compute pointwise t-statistics between conditions.
2.  Threshold the t-values to identify clusters of contiguous
    significant points.
3.  For each cluster, compute a cluster-level statistic (sum of
    t-values).
4.  Build a null distribution by permuting condition labels and
    repeating steps 1-3.
5.  Compare observed cluster statistics to the permutation distribution.

``` r

# Run cluster-based permutation test
cluster_result <- clusterPermutationTest(
  pe,
  condition1 = 1:20,
  condition2 = 21:40,
  n_permutations = 1000,
  cluster_threshold = 0.05,
  tail = 0,        # two-tailed
  seed = 42
)

# Examine results
cat("Number of clusters found:", length(cluster_result$clusters), "\n")
cat("Cluster p-values:", cluster_result$cluster_p, "\n")

# The cluster_mask indicates which time-channel points belong to
# significant clusters (p < 0.05)
sig_points <- sum(cluster_result$cluster_mask)
cat("Significant time-channel points:", sig_points, "\n")
```

### One-Sample Cluster Test

When only one condition is provided, the test uses sign-flipping
permutations to test against zero.

``` r

cluster_one <- clusterPermutationTest(
  pe,
  condition1 = 21:40,  # stimulus epochs only
  n_permutations = 500,
  seed = 42
)
```

### Examining the Permutation Distribution

``` r

# The permutation distribution is available for inspection
hist(cluster_result$permutation_distribution,
     breaks = 30,
     main = "Permutation Distribution of Max Cluster Statistic",
     xlab = "Max cluster statistic")

# Add observed cluster statistics
abline(v = abs(cluster_result$cluster_stats), col = "red", lwd = 2)
```

## Effect Sizes

The
[`effectSize()`](https://x-biosignal.github.io/PhysioAnalysis/reference/effectSize.md)
function computes Cohen’s d at each time point and channel, with 95%
confidence intervals.

``` r

# One-sample effect size (d = mean / sd)
es_one <- effectSize(pe, condition1 = 21:40)

# Two-sample effect size (pooled SD)
es_two <- effectSize(pe, condition1 = 1:20, condition2 = 21:40)

# Glass's delta (using control group SD)
es_glass <- effectSize(pe, condition1 = 21:40, condition2 = 1:20,
                       pooled = FALSE)

# Examine effect sizes with confidence intervals
str(es_two)
# d:        Cohen's d matrix (time x channel)
# ci_lower: Lower 95% CI
# ci_upper: Upper 95% CI
```

## Bootstrap Confidence Intervals

The
[`bootstrapCI()`](https://x-biosignal.github.io/PhysioAnalysis/reference/bootstrapCI.md)
function computes bootstrap confidence intervals for the mean across
epochs. This is particularly useful for ERP analysis.

``` r

# Bootstrap CI with 1000 iterations
boot_result <- bootstrapCI(pe,
                           n_bootstrap = 1000,
                           ci_level = 0.95,
                           condition = 21:40,  # stimulus epochs
                           seed = 42)

# The result contains the mean and CI bounds
str(boot_result)

# Plot the ERP with bootstrap CI for channel 1
plot(boot_result$times, boot_result$mean[, 1],
     type = "l", lwd = 2,
     xlab = "Time (s)", ylab = "Amplitude",
     main = "ERP with 95% Bootstrap CI")
polygon(c(boot_result$times, rev(boot_result$times)),
        c(boot_result$ci_lower[, 1], rev(boot_result$ci_upper[, 1])),
        col = rgb(0, 0, 1, 0.2), border = NA)
```

## Multiple Comparison Correction

The
[`correctPValues()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correctPValues.md)
function applies standard correction methods to p-value matrices or
vectors.

``` r

# Get uncorrected p-values from t-test
p_uncorrected <- result_two$p_values

# False Discovery Rate (Benjamini-Hochberg)
p_fdr <- correctPValues(p_uncorrected, method = "fdr")

# Bonferroni correction
p_bonf <- correctPValues(p_uncorrected, method = "bonferroni")

# Holm's step-down procedure
p_holm <- correctPValues(p_uncorrected, method = "holm")

# Compare number of significant results
cat("Uncorrected:", sum(p_uncorrected < 0.05, na.rm = TRUE), "\n")
cat("FDR:", sum(p_fdr < 0.05, na.rm = TRUE), "\n")
cat("Bonferroni:", sum(p_bonf < 0.05, na.rm = TRUE), "\n")
cat("Holm:", sum(p_holm < 0.05, na.rm = TRUE), "\n")
```

### Choosing a Correction Method

| Method | Approach | When to use |
|----|----|----|
| `"fdr"` / `"bh"` | Benjamini-Hochberg | Exploratory analyses; controls FDR |
| `"bonferroni"` | Bonferroni | Conservative; small number of tests |
| `"holm"` | Holm step-down | Moderate; uniformly more powerful |
|  |  | than Bonferroni |
| Cluster-based | Permutation distribution | Recommended for spatiotemporal data |

## Finding Significant Time Windows

The
[`findSignificantWindows()`](https://x-biosignal.github.io/PhysioAnalysis/reference/findSignificantWindows.md)
function identifies contiguous time periods with significant effects,
useful for summarizing results.

``` r

# Get the time vector
times <- seq(-0.2, 0.2, length.out = n_time)

# Find significant windows for channel 1 (FDR-corrected)
windows <- findSignificantWindows(
  p_values = p_fdr[, 1],
  times = times,
  alpha = 0.05,
  min_duration = 0.02  # at least 20 ms
)

# Display the significant windows
print(windows)
# Columns: start, end, duration, min_p
```

## Complete Workflow Example

A typical statistical analysis workflow for ERP data:

``` r

# 1. Start with epoched data
#    (assuming pe is already epoched with condition labels)

# 2. Compute pointwise t-test
t_result <- tTestEpochs(pe, condition1 = 1:20, condition2 = 21:40)

# 3. Apply cluster-based permutation test for correction
cluster_result <- clusterPermutationTest(
  pe,
  condition1 = 1:20,
  condition2 = 21:40,
  n_permutations = 1000,
  seed = 42
)

# 4. Compute effect sizes
es <- effectSize(pe, condition1 = 1:20, condition2 = 21:40)

# 5. Compute bootstrap CI for each condition
boot_ctrl <- bootstrapCI(pe, condition = 1:20, seed = 42)
boot_stim <- bootstrapCI(pe, condition = 21:40, seed = 42)

# 6. Report significant clusters
for (i in seq_along(cluster_result$clusters)) {
  cat(sprintf("Cluster %d: stat = %.2f, p = %.4f\n",
              i, cluster_result$cluster_stats[i],
              cluster_result$cluster_p[i]))
}
```

## References

- Maris, E. & Oostenveld, R. (2007). “Nonparametric statistical testing
  of EEG- and MEG-data.” *Journal of Neuroscience Methods*, 164(1),
  177-190.
- Luck, S.J. (2014). *An Introduction to the Event-Related Potential
  Technique*. 2nd ed. MIT Press.

## Session Info

``` r

sessionInfo()
```
