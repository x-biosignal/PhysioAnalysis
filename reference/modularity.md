# Compute network modularity

Calculates the modularity of a network given a community partition.

## Usage

``` r
modularity(adjacency, communities)
```

## Arguments

- adjacency:

  An adjacency matrix.

- communities:

  A vector of community assignments for each node.

## Value

Modularity value (-0.5 to 1).

## Examples

``` r
# Two clear communities
adj <- matrix(0, 6, 6)
adj[1:3, 1:3] <- 1
adj[4:6, 4:6] <- 1
diag(adj) <- 0
communities <- c(1, 1, 1, 2, 2, 2)
modularity(adj, communities)
#> [1] 0.5
```
