# Compute shortest path length

Calculates the shortest path length between all node pairs using
Floyd-Warshall.

## Usage

``` r
pathLength(adjacency, weighted = FALSE, use_cpp = TRUE)
```

## Arguments

- adjacency:

  An adjacency matrix.

- weighted:

  If TRUE, uses edge weights as distances.

- use_cpp:

  If TRUE (default), use the compiled C++ backend, falling back to R if
  unavailable.

## Value

A matrix of shortest path lengths.

## Examples

``` r
adj <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
pathLength(adj)
#>      [,1] [,2] [,3]
#> [1,]    0    1    2
#> [2,]    1    0    1
#> [3,]    2    1    0
```
