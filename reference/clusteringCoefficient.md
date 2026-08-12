# Compute clustering coefficient

Calculates the local clustering coefficient for each node.

## Usage

``` r
clusteringCoefficient(
  adjacency,
  weighted = FALSE,
  use_cpp = TRUE,
  n_cores = 1L
)
```

## Arguments

- adjacency:

  An adjacency matrix.

- weighted:

  If TRUE, uses weighted clustering coefficient.

- use_cpp:

  If TRUE (default), use the compiled C++ backend, falling back to R if
  unavailable.

- n_cores:

  Number of OpenMP threads for the C++ backend (default 1).

## Value

A numeric vector of clustering coefficients.

## Examples

``` r
# Fully connected triangle
adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
clusteringCoefficient(adj)  # Should be 1 for all nodes
#> [1] 1 1 1
```
