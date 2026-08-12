# Spectral decomposition of graph Laplacian

Computes eigenvalues and eigenvectors of the graph Laplacian.

## Usage

``` r
spectralDecomposition(adjacency, normalized = TRUE, n_components = NULL)
```

## Arguments

- adjacency:

  An adjacency matrix.

- normalized:

  If TRUE, uses normalized Laplacian.

- n_components:

  Number of components to return. If NULL, returns all.

## Value

A list with eigenvalues and eigenvectors.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
spec <- spectralDecomposition(adj)
```
