# Compute betweenness centrality

Calculates betweenness centrality for each node.

## Usage

``` r
betweennessCentrality(
  adjacency,
  normalized = TRUE,
  use_cpp = TRUE,
  n_cores = 1L
)
```

## Arguments

- adjacency:

  An adjacency matrix.

- normalized:

  If TRUE, normalizes by (n-1)(n-2)/2.

- use_cpp:

  If TRUE (default), use the compiled C++ backend, falling back to R if
  unavailable.

- n_cores:

  Number of OpenMP threads for the C++ backend (default 1).

## Value

A numeric vector of betweenness centrality values.

## Examples

``` r
# Star network - center node should have highest betweenness
adj <- matrix(0, 4, 4)
adj[1, 2:4] <- 1
adj[2:4, 1] <- 1
betweennessCentrality(adj)
#> [1] 1 0 0 0
```
