# Create adjacency matrix from connectivity

Converts a connectivity matrix into an adjacency matrix for network
analysis.

## Usage

``` r
adjacencyMatrix(connectivity, threshold = NULL, absolute = FALSE)
```

## Arguments

- connectivity:

  A connectivity matrix or result from connectivityMatrix().

- threshold:

  Threshold value for binarization. If NULL, keeps weighted edges.

- absolute:

  If TRUE, uses absolute values before thresholding.

## Value

A square adjacency matrix.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 4), nrow = 500, ncol = 4)),
  samplingRate = 100
)
conn <- correlationMatrix(pe)
adj <- adjacencyMatrix(conn, threshold = 0.3)
```
