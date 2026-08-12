# Plot network metrics

Creates a bar plot of network metrics for each node.

## Usage

``` r
plotNetworkMetrics(
  adjacency,
  metrics = c("degree", "clustering", "betweenness"),
  node_names = NULL,
  title = "Network Metrics"
)
```

## Arguments

- adjacency:

  An adjacency matrix.

- metrics:

  Vector of metrics to compute: "degree", "clustering", "betweenness",
  "eigenvector", "local_efficiency".

- node_names:

  Optional vector of node names.

- title:

  Plot title.

## Value

A ggplot object.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0), 4, 4)
plotNetworkMetrics(adj, node_names = c("Fz", "Cz", "Pz", "Oz"))
```
