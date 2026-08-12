# Plot dynamic connectivity

Creates a visualization of time-varying connectivity.

## Usage

``` r
plotDynamicConnectivity(
  dyn_conn,
  node_pair = "mean",
  title = "Dynamic Connectivity"
)
```

## Arguments

- dyn_conn:

  Result from slidingWindowConnectivity().

- node_pair:

  Vector of two node indices to plot, or "mean" for average.

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
plotDynamicConnectivity(dyn_conn, node_pair = c(1, 2))
```
