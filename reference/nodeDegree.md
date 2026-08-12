# Compute node degree

Calculates the degree (number of connections) for each node.

## Usage

``` r
nodeDegree(adjacency, weighted = FALSE, use_cpp = TRUE)
```

## Arguments

- adjacency:

  An adjacency matrix (binary or weighted).

- weighted:

  If TRUE, returns weighted degree (strength).

- use_cpp:

  If TRUE (default), use the compiled C++ backend, falling back to R if
  unavailable.

## Value

A numeric vector of node degrees.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3)
nodeDegree(adj)
#> [1] 2 1 1
```
