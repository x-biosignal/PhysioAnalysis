# Network-Based Statistic (NBS)

Identifies connected subnetworks of edges that differ between two groups
(or two within-subject conditions) of connectivity matrices, with
family-wise error control from a permutation null on the largest
suprathreshold component (Zalesky, Fornito & Bullmore 2010). Each edge
is first tested with a mass-univariate t-test; edges whose statistic
exceeds `thresh` form a suprathreshold graph whose connected components
are the candidate subnetworks. Their size (edge count) or mass (summed
statistic exceedance, Smith 2009) is compared against the permutation
distribution of the largest component to obtain FWER-corrected component
p-values.

## Usage

``` r
networkBasedStatistic(
  mats_group1,
  mats_group2,
  thresh = 3,
  n_perm = 1000L,
  tail = c("both", "right", "left"),
  paired = FALSE,
  component = c("size", "mass"),
  directed = FALSE,
  alpha = 0.05,
  seed = NULL
)
```

## Arguments

- mats_group1:

  A list of n x n connectivity matrices or an n x n x N array (group 1,
  or condition 1 when `paired`).

- mats_group2:

  The corresponding matrices for group 2 (or condition 2 for a paired
  design).

- thresh:

  Primary edge-statistic (t) threshold (default: 3).

- n_perm:

  Number of permutations (default: 1000).

- tail:

  `"both"`, `"right"` (group1 \> group2), or `"left"`.

- paired:

  Logical; `TRUE` for a within-subject (paired) design, in which case
  the two inputs are paired condition matrices (default: `FALSE`).

- component:

  Component measure: `"size"` (edge count) or `"mass"` (summed statistic
  exceedance).

- directed:

  Logical; use all off-diagonal (directed) edges rather than the upper
  triangle (default: `FALSE`).

- alpha:

  Significance level for the returned adjacency mask (default: 0.05).

- seed:

  Optional RNG seed for reproducible permutations.

## Value

A list with `components` (a data.frame of component `size`, `mass`, and
`p_value`), `component_edges` (a list of node-pair matrices, one per
component), the `adjacency` mask of edges in significant components, the
`edge_stats` matrix, the `suprathreshold` mask, the permutation
`null_distribution`, and the settings used.

## References

Zalesky, A., Fornito, A., & Bullmore, E. T. (2010). Network-based
statistic: identifying differences in brain networks. NeuroImage, 53(4),
1197-1207.

## See also

[`plotNBSnetwork()`](https://x-biosignal.github.io/PhysioAnalysis/reference/plotNBSnetwork.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g1 <- replicate(20, matrix(rnorm(400), 20), simplify = FALSE)
g2 <- replicate(20, matrix(rnorm(400), 20), simplify = FALSE)
nbs <- networkBasedStatistic(g1, g2, thresh = 3, n_perm = 500)
nbs$components
} # }
```
