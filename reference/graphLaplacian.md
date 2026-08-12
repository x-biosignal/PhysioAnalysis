# Compute graph Laplacian

Calculates the graph Laplacian matrix.

## Usage

``` r
graphLaplacian(adjacency, normalized = FALSE)
```

## Arguments

- adjacency:

  An adjacency matrix.

- normalized:

  If TRUE, returns normalized Laplacian.

## Value

The Laplacian matrix.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
L <- graphLaplacian(adj)
```
