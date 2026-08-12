# Compute global efficiency

Calculates the global efficiency of a network.

## Usage

``` r
globalEfficiency(adjacency, weighted = FALSE, use_cpp = TRUE)
```

## Arguments

- adjacency:

  An adjacency matrix.

- weighted:

  If TRUE, uses weighted paths.

- use_cpp:

  If TRUE (default), use the compiled C++ backend, falling back to R if
  unavailable.

## Value

Global efficiency value (0-1).

## Examples

``` r
# Fully connected - maximum efficiency
adj <- matrix(1, 4, 4)
diag(adj) <- 0
globalEfficiency(adj)
#> [1] 1
```
