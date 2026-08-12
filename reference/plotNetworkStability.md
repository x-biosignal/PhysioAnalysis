# Plot network stability over time

Visualizes the temporal stability of network topology.

## Usage

``` r
plotNetworkStability(stability, title = "Network Stability")
```

## Arguments

- stability:

  Result from temporalStability().

- title:

  Plot title.

## Value

A ggplot object.

## Examples

``` r
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)),
  samplingRate = 100
)
dyn_conn <- slidingWindowConnectivity(pe, window_size = 200, step = 50)
stab <- temporalStability(dyn_conn)
plotNetworkStability(stab)
```
