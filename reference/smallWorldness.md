# Compute small-worldness

Calculates the small-world coefficient (sigma) of a network.

## Usage

``` r
smallWorldness(adjacency, n_rand = 100, n_cores = NULL)
```

## Arguments

- adjacency:

  An adjacency matrix.

- n_rand:

  Number of random networks for comparison.

- n_cores:

  Number of cores for parallel processing. Default NULL uses sequential
  processing.

## Value

A list with small-worldness metrics.

## Examples

``` r
# Create a small-world-like network
set.seed(123)
n <- 20
adj <- matrix(0, n, n)
for (i in 1:n) {
  adj[i, ((i) %% n) + 1] <- 1
  adj[i, ((i + 1) %% n) + 1] <- 1
}
adj <- adj + t(adj)
adj[adj > 0] <- 1
# Add some random long-range connections
adj[sample(which(adj == 0 & row(adj) < col(adj)), 5)] <- 1
adj <- adj + t(adj)
adj[adj > 0] <- 1
diag(adj) <- 0
sw <- smallWorldness(adj, n_rand = 10)
```
