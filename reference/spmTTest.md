# Statistical Parametric Mapping (SPM1D) for Biomechanics

Functions for statistical analysis of continuous waveform data using
Statistical Parametric Mapping (SPM) methodology adapted from
neuroimaging. These methods test hypotheses over entire waveforms rather
than discrete points. SPM t-test for waveform comparison

## Usage

``` r
spmTTest(x, group1 = NULL, group2 = NULL, alpha = 0.05, two_tailed = TRUE)
```

## Arguments

- x:

  A PhysioExperiment object or matrix (time x observations).

- group1:

  Indices for first group (for two-sample test).

- group2:

  Indices for second group. If NULL, performs one-sample test.

- alpha:

  Significance level (default: 0.05).

- two_tailed:

  Logical; if TRUE, performs two-tailed test.

## Value

A list of class "spm_result" containing:

- t:

  T-statistic at each time point

- threshold:

  Critical threshold from RFT

- clusters:

  Significant clusters (start, end, extent, p-value)

- p_values:

  Pointwise p-values

- alpha:

  Significance level used

## Details

Performs a t-test at each time point and computes an SPM t-statistic map
with Random Field Theory (RFT) correction for multiple comparisons.

SPM analyzes continuous biomechanical waveforms (e.g., joint angles,
moments) by computing t-statistics at each time point and using Random
Field Theory to control family-wise error rate across the entire
waveform.

## References

Pataky TC (2012). One-dimensional statistical parametric mapping in
Python. Computer Methods in Biomechanics and Biomedical Engineering.

## Examples

``` r
# Create example gait data (100 time points x 20 subjects)
set.seed(123)
# Group 1: normal gait
g1 <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
# Group 2: altered gait (effect at 40-60% of cycle)
g2 <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
g2[40:60, ] <- g2[40:60, ] + 1.5

data <- cbind(g1, g2)
pe <- PhysioExperiment(
  assays = list(values = data),
  samplingRate = 100
)

# Two-sample SPM t-test
result <- spmTTest(pe, group1 = 1:10, group2 = 11:20)
print(result)
#> SPM Analysis Result
#> ==================
#> Test type: two-sample
#> Time points: 100
#> Alpha: 0.050
#> Threshold: 4.375
#> FWHM (smoothness): 2.24
#> Resel count: 44.24
#> 
#> Significant clusters: 5
#>   Cluster 1: [44-45] extent=2, p=0.0001
#>   Cluster 2: [50-50] extent=1, p=0.0113
#>   Cluster 3: [52-52] extent=1, p=0.0113
#>   Cluster 4: [54-54] extent=1, p=0.0113
#>   Cluster 5: [58-58] extent=1, p=0.0113
```
