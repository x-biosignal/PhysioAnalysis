# Binarize network

Converts a weighted adjacency matrix to a binary network.

## Usage

``` r
binarizeNetwork(adjacency, threshold = 0)
```

## Arguments

- adjacency:

  A weighted adjacency matrix.

- threshold:

  Threshold for binarization. Default 0 (any non-zero edge).

## Value

A binary adjacency matrix (0s and 1s).

## Examples

``` r
set.seed(123)
adj <- matrix(c(0, 0.5, 0.3, 0.5, 0, 0.8, 0.3, 0.8, 0), 3, 3)
bin <- binarizeNetwork(adj)
```
