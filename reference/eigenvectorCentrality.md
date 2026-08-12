# Compute eigenvector centrality

Calculates eigenvector centrality for each node.

## Usage

``` r
eigenvectorCentrality(adjacency, max_iter = 100, tol = 1e-06, use_cpp = TRUE)
```

## Arguments

- adjacency:

  An adjacency matrix.

- max_iter:

  Maximum iterations for power method.

- tol:

  Convergence tolerance.

- use_cpp:

  If TRUE (default), use the compiled C++ backend, falling back to R if
  unavailable.

## Value

A numeric vector of eigenvector centrality values.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
eigenvectorCentrality(adj)
#> [1] 0.5773503 0.5773503 0.5773503
```
