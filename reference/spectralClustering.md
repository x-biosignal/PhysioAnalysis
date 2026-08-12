# Spectral clustering of network nodes

Performs spectral clustering on the network.

## Usage

``` r
spectralClustering(adjacency, n_clusters = 2, normalized = TRUE)
```

## Arguments

- adjacency:

  An adjacency matrix.

- n_clusters:

  Number of clusters.

- normalized:

  If TRUE, uses normalized Laplacian.

## Value

A list with cluster assignments and spectral embedding.

## Examples

``` r
# Create block-structured network
adj <- matrix(0, 6, 6)
adj[1:3, 1:3] <- 1
adj[4:6, 4:6] <- 1
adj[3, 4] <- adj[4, 3] <- 0.5
diag(adj) <- 0
clusters <- spectralClustering(adj, n_clusters = 2)
```
