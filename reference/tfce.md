# Threshold-Free Cluster Enhancement (TFCE)

Enhances a statistic map by integrating, at every point, the extent of
the suprathreshold cluster containing that point raised to the power `E`
times the threshold height raised to the power `H`, over all thresholds
(Smith & Nichols 2009). This rewards signal that is either spatially
broad or locally high without committing to a single cluster-forming
threshold. Works on a time-by-channel (or time-by-frequency) matrix, or
a 1-D vector, with time (and channel) adjacency.

## Usage

``` r
tfce(stat_map, E = 0.5, H = 2, dh = NULL, adjacency = NULL, tail = 0)
```

## Arguments

- stat_map:

  A numeric matrix (time x channel) or a numeric vector.

- E:

  Extent exponent (default: 0.5).

- H:

  Height exponent (default: 2).

- dh:

  Threshold step. `NULL` (default) integrates exactly over the distinct
  statistic levels; a numeric value uses a regular grid.

- adjacency:

  Optional channel-adjacency matrix (columns of `stat_map`). `NULL` uses
  neighbouring-column (4-connectivity) adjacency.

- tail:

  `1` for positive effects, `-1` for negative, or `0` (default) for a
  signed two-sided map (positive TFCE minus negative TFCE).

## Value

A TFCE map with the same shape as `stat_map`.

## References

Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
enhancement: addressing problems of smoothing, threshold dependence and
localisation in cluster inference. NeuroImage, 44(1), 83-98.

## See also

[`clusterPermutationTest()`](https://x-biosignal.github.io/PhysioAnalysis/reference/clusterPermutationTest.md)

## Examples

``` r
m <- matrix(0, 100, 1); m[20:80, 1] <- 2; m[90:93, 1] <- 6
tf <- tfce(m)
c(broad = tf[50, 1], narrow = tf[91, 1])
#>     broad    narrow 
#>  20.82733 144.00000 
```
