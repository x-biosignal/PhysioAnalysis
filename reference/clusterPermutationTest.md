# Cluster-based permutation test

Performs cluster-based permutation testing for multiple comparison
correction. Identifies clusters of significant effects and computes
cluster-level p-values.

## Usage

``` r
clusterPermutationTest(
  x,
  condition1 = NULL,
  condition2 = NULL,
  n_permutations = 1000L,
  cluster_threshold = 0.05,
  tail = 0L,
  seed = NULL,
  n_cores = NULL,
  method = c("cluster", "tfce"),
  tfce_E = 0.5,
  tfce_H = 2,
  tfce_dh = NULL,
  adjacency = NULL
)
```

## Arguments

- x:

  An epoched PhysioExperiment object (4D data).

- condition1:

  Indices for first condition.

- condition2:

  Indices for second condition. If NULL, tests against zero.

- n_permutations:

  Number of permutations (default: 1000).

- cluster_threshold:

  Initial threshold for cluster formation (p-value).

- tail:

  Test type: 0 = two-tailed, 1 = right tail, -1 = left tail.

- seed:

  Random seed for reproducibility.

- n_cores:

  Number of cores for parallel processing. Default NULL uses sequential
  processing. Set to parallel::detectCores() - 1 for maximum speed.

- method:

  `"cluster"` (default) for cluster-mass inference, or `"tfce"` for
  Threshold-Free Cluster Enhancement
  ([`tfce()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tfce.md))
  with a per-point max-statistic FWER correction.

- tfce_E, tfce_H:

  Extent and height exponents for `method = "tfce"` (defaults: 0.5 and
  2).

- tfce_dh:

  Threshold step for TFCE (`NULL` integrates exactly over the distinct
  statistic levels).

- adjacency:

  Optional channel-adjacency matrix for TFCE; `NULL` uses
  neighbouring-channel adjacency.

## Value

A list containing:

- clusters:

  List of identified clusters with their statistics

- cluster_p:

  P-values for each cluster

- t_obs:

  Observed t-values matrix

- cluster_mask:

  Logical matrix indicating cluster membership

## References

Maris, E. & Oostenveld, R. (2007). "Nonparametric statistical testing of
EEG- and MEG-data." *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

## See also

[`tTestEpochs()`](https://x-biosignal.github.io/PhysioAnalysis/reference/tTestEpochs.md)
for pointwise t-tests,
[`correctPValues()`](https://x-biosignal.github.io/PhysioAnalysis/reference/correctPValues.md)
for traditional multiple comparison correction,
[`findSignificantWindows()`](https://x-biosignal.github.io/PhysioAnalysis/reference/findSignificantWindows.md)
for identifying significant time windows.

## Examples

``` r
# Create example epoched data
set.seed(123)
epochs <- array(rnorm(50 * 4 * 20 * 1), dim = c(50, 4, 20, 1))
# Add effect to condition 2
epochs[20:30, 1:2, 11:20, 1] <- epochs[20:30, 1:2, 11:20, 1] + 1
pe <- PhysioExperiment(
  assays = list(epoched = epochs),
  samplingRate = 100
)
# Cluster permutation test
result <- clusterPermutationTest(pe, 1:10, 11:20, n_permutations = 100)

# With parallel processing (faster for large datasets)
# result <- clusterPermutationTest(pe, 1:10, 11:20, n_cores = 4)
```
