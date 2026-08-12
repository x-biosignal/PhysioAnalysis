# Threshold network by density

Thresholds a connectivity matrix to achieve a target network density.

## Usage

``` r
thresholdNetwork(connectivity, density = 0.2, absolute = TRUE)
```

## Arguments

- connectivity:

  A connectivity matrix.

- density:

  Target density (proportion of edges to keep, 0-1).

- absolute:

  If TRUE, uses absolute values for ranking.

## Value

A thresholded adjacency matrix.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 4), nrow = 500, ncol = 4)),
  samplingRate = 100
)
conn <- correlationMatrix(pe)
# Keep top 30% of connections
adj <- thresholdNetwork(conn, density = 0.3)
```
