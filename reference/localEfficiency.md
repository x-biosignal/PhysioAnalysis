# Compute local efficiency

Calculates the local efficiency for each node.

## Usage

``` r
localEfficiency(adjacency, weighted = FALSE)
```

## Arguments

- adjacency:

  An adjacency matrix.

- weighted:

  If TRUE, uses weighted paths.

## Value

A numeric vector of local efficiency values.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
localEfficiency(adj)
#> [1] 1 1 1
```
